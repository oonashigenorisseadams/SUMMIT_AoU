version 1.0

# merge_ld_panels.wdl
#
# Gather stage for aou_ld_panels_v9. Merges the 22 per-chromosome filtered
# filesets produced by filter_ld_panels.wdl into one genome-wide fileset
# per panel, via plink2 --pmerge-list.
#
# Why this stage exists: SUMMIT's randomized genome-wide LD-score estimator
# opens one persistent PgenReader per analysis and accumulates
#   W_kv = sum_j sqrt(a_jk) X_j z_jv
# as a sum over every variant in the panel, with a genome-wide null
# subtraction (sum_j a_jk)/d. It takes one --geno prefix resolving to exactly
# one complete trio, so it cannot consume 22 per-chromosome filesets.
# (The deterministic windowed estimator, --ld-wind-kb, masks to an exact
# base-pair window and never crosses chromosomes, so it would not need this
# stage -- but it is a different estimand with a different standardization
# convention (ddof=0 vs ddof=1) and is not the h2/rg path.)
#
# Inputs are PANEL-MAJOR: pgen_by_panel[i] is the 22 chromosomes of panels[i],
# in chromosome order. This shape works two ways:
#
#   standalone  -- supply panel-major arrays in the inputs JSON
#   chained     -- inside filter_ld_panels.wdl, after the chromosome scatter:
#
#                    scatter (i in range(length(panels))) {
#                      call merge_panel {
#                        input:
#                          panel  = panels[i],
#                          chroms = chroms,
#                          pgen   = transpose(filter.out_pgen)[i],
#                          pvar   = transpose(filter.out_pvar)[i],
#                          psam   = transpose(filter.out_psam)[i],
#                          plink2 = plink2
#                      }
#                    }
#
#                  transpose() is only valid if every filter shard emits its
#                  panel outputs in the SAME order. Construct those output
#                  arrays from the panel array, not from glob().

workflow ld_panel_merge {

  input {
    Array[String] panels
    Array[String] chroms = ["1","2","3","4","5","6","7","8","9","10","11",
                            "12","13","14","15","16","17","18","19","20","21","22"]

    Array[Array[File]] pgen_by_panel
    Array[Array[File]] pvar_by_panel
    Array[Array[File]] psam_by_panel

    File plink2
    String docker_image = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
  }

  scatter (i in range(length(panels))) {
    call merge_panel {
      input:
        panel        = panels[i],
        chroms       = chroms,
        pgen         = pgen_by_panel[i],
        pvar         = pvar_by_panel[i],
        psam         = psam_by_panel[i],
        plink2       = plink2,
        docker_image = docker_image
    }
  }

  output {
    Array[File] merged_pgen     = merge_panel.merged_pgen
    Array[File] merged_pvar     = merge_panel.merged_pvar
    Array[File] merged_psam     = merge_panel.merged_psam
    Array[File] merged_snplist  = merge_panel.snplist
    Array[File] merge_manifests = merge_panel.manifest
  }
}


task merge_panel {

  input {
    String panel
    Array[String] chroms
    Array[File] pgen
    Array[File] pvar
    Array[File] psam

    File plink2
    String docker_image

    Int cpu    = 8
    Int mem_gb = 32

    # Derived so the two can't drift apart. plink2 is given the container's
    # memory less a 4 GB margin for the JVM-less shell, awk/sort passes, and
    # allocator slack. If a large panel OOMs, raise mem_gb -- do not reach for
    # --sort-vars, which costs more memory, not less.
    Int plink_mem_mb = (mem_gb - 4) * 1024
  }

  # Inputs land once, output is roughly the same size again, plus slack.
  # The 3x multiplier is also buying throughput: PD-SSD scales at
  # 240 MiB/s + 0.48 MiB/s per GiB provisioned.
  #
  # Measured genome-wide input sizes (production run, 2026-09-09):
  #   afr 84.4 GB -> ~303 GB disk   (largest panel; NOT eur)
  #   amr 57.0 GB -> ~221 GB
  #   eur 45.0 GB -> ~185 GB        (capped at N=75,000)
  #   eas  7.1 GB, sas 4.1 GB, mid 1.6 GB
  Float in_gb   = size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB")
  Int   disk_gb = ceil(in_gb * 3 + 50)

  command <<<
    set -euo pipefail

    PLINK="~{plink2}"
    chmod 750 "$PLINK"
    mkdir -p out

    # ---------------------------------------------------------------------
    # Build the merge list.
    #
    # Order comes straight from the WDL array, which preserves scatter order.
    # Do NOT sort these paths -- lexicographic order puts chr10 before chr2
    # and breaks --pmerge-list's chromosome-contiguity requirement.
    # ---------------------------------------------------------------------
    printf '%s\n' ~{sep=" " pgen} | sed 's/\.pgen$//' > merge_list.txt

    printf '%s\n' ~{sep=" " chroms} > expected_chroms.txt

    n_in=$(awk 'END{print NR}' merge_list.txt)
    n_exp=$(awk 'END{print NR}' expected_chroms.txt)
    echo "[~{panel}] merge_list entries: ${n_in} (expected ${n_exp})"
    test "${n_in}" -eq "${n_exp}"

    # Audit: entry k must be the fileset for expected chromosome k.
    #
    # filter_ld_panels.wdl writes --out out/<panel>.chr<C>, so a localized
    # prefix looks like
    #   /mnt/disks/cromwell_root/.../anc_afr.chr7
    # Strip the directory, then everything up to and including the last
    # ".chr", leaving "chr7". Panel names may themselves contain dots or
    # digits, which is why this anchors on the final ".chr" rather than
    # splitting on the first dot.
    # ---------------------------------------------------------------------
    sed 's|.*/||; s|^.*\.chr|chr|' merge_list.txt > basenames.txt
    paste expected_chroms.txt basenames.txt > order_audit.tsv
    cat order_audit.tsv
    awk '$2 != "chr" $1 { print "ORDER MISMATCH: expected chr" $1 ", got " $2 > "/dev/stderr"; exit 1 }' \
      order_audit.tsv

    # ---------------------------------------------------------------------
    # Every trio member must be present. Only the .pgen paths were used to
    # build the list; .pvar and .psam are passed in purely so Cromwell
    # localizes them alongside.
    # ---------------------------------------------------------------------
    while read -r p; do
      for ext in pgen pvar psam; do
        test -s "${p}.${ext}" || { echo "missing ${p}.${ext}" >&2; exit 1; }
      done
    done < merge_list.txt

    # ---------------------------------------------------------------------
    # Sample sets must be identical across all 22 chromosomes. They come from
    # one keep file so they should be; a mismatch silently turns a
    # concatenation into a real merge with missing calls.
    #
    # grep -c exits 1 on zero matches and trips pipefail -- use awk | wc -l.
    # ---------------------------------------------------------------------
    ref_hash=""
    n_samples=0
    while read -r p; do
      awk '!/^#/' "${p}.psam" | LC_ALL=C sort > ids.txt
      h=$(md5sum ids.txt | cut -d' ' -f1)
      if [ -z "${ref_hash}" ]; then
        ref_hash="${h}"
        n_samples=$(awk 'END{print NR}' ids.txt)
      fi
      if [ "${h}" != "${ref_hash}" ]; then
        echo "SAMPLE SET MISMATCH at ${p}" >&2
        exit 1
      fi
    done < merge_list.txt
    echo "[~{panel}] sample set consistent across all inputs: n=${n_samples}"

    # Expected genome-wide M, for the post-merge check.
    n_var_in=0
    while read -r p; do
      n=$(awk '!/^#/' "${p}.pvar" | wc -l)
      n_var_in=$(( n_var_in + n ))
    done < merge_list.txt
    echo "[~{panel}] summed input variants: ${n_var_in}"

    # ---------------------------------------------------------------------
    # Merge.
    #
    # Sample sets are identical and variant sets are disjoint (one per
    # chromosome), so this is a concatenation -- no --merge-mode needed.
    # If plink2 complains about variant ordering, add --sort-vars; it should
    # not be necessary given the ordering audit above, and it costs memory.
    #
    # Output stays PGEN: SUMMIT reads .pgen/.pvar/.psam natively. Do not add
    # --make-bed. The PVAR must remain plain text -- SUMMIT does not read
    # .pvar.zst.
    # ---------------------------------------------------------------------
    t0=$(date +%s)
    "$PLINK" --pmerge-list merge_list.txt pfile \
      --make-pgen \
      --threads ~{cpu} --memory ~{plink_mem_mb} \
      --out out/~{panel}
    t1=$(date +%s)

    # ---------------------------------------------------------------------
    # Post-merge checks.
    #
    # SUMMIT requires unique PVAR variant IDs and unique PSAM sample IDs, and
    # chromosome-contiguous variant order. The first two are asserted here;
    # the third follows from the merge-list order audit above.
    # ---------------------------------------------------------------------
    n_var_out=$(awk '!/^#/' out/~{panel}.pvar | wc -l)
    n_samp_out=$(awk '!/^#/' out/~{panel}.psam | wc -l)

    test "${n_var_out}" -eq "${n_var_in}"
    test "${n_samp_out}" -eq "${n_samples}"

    # pvar columns: #CHROM POS ID REF ALT
    n_dup=$(awk '!/^#/{print $3}' out/~{panel}.pvar | LC_ALL=C sort | uniq -d | wc -l)
    test "${n_dup}" -eq 0

    # Canonical variant list, in panel order. This is what any downstream
    # .annot file has to match: a thin annotation needs exactly one row per
    # genotype variant, and a full annotation is aligned by SNP ID and then
    # required to agree on chromosome and base-pair coordinate.
    awk '!/^#/{print $3}' out/~{panel}.pvar > out/~{panel}.snplist

    bytes_out=$(du -cb out/~{panel}.pgen out/~{panel}.pvar out/~{panel}.psam | tail -1 | cut -f1)

    {
      printf 'panel\tn_chrom\tn_samples\tn_variants\tn_dup_ids\tbytes_out\tseconds\n'
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "~{panel}" "${n_in}" "${n_samp_out}" "${n_var_out}" "${n_dup}" "${bytes_out}" "$(( t1 - t0 ))"
    } | tee out/~{panel}.merge_manifest.tsv
  >>>

  output {
    File merged_pgen = "out/~{panel}.pgen"
    File merged_pvar = "out/~{panel}.pvar"
    File merged_psam = "out/~{panel}.psam"
    File snplist     = "out/~{panel}.snplist"
    File manifest    = "out/~{panel}.merge_manifest.tsv"
  }

  runtime {
    docker: docker_image
    memory: mem_gb + " GB"
    cpu: cpu
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: 0
  }
}
