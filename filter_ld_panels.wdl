version 1.0

## aou_ld_panels_v9 — plink2 filter for SUMMIT reference panels
##
## First run: chroms = ["22"], keep_files = [anc_eur.keep, anc_mid.keep].

workflow ld_panel_filter {
    input {
        # Parallel arrays, same length and order. index i = one chromosome.
        Array[String] chroms
        Array[File]   pgen_files
        Array[File]   pvar_files
        Array[File]   psam_files

        # Panel definition. Output names are derived from basename(keep, ".keep"),
        # so anc_eur.keep -> anc_eur.chr22.{pgen,pvar,psam}
        Array[File]   keep_files

        File plink2

        # QC — per panel, not intersected
        Float  maf  = 0.01
        Float  geno = 0.05
        String hwe  = "1e-6 keep-fewhet"

        # Runtime
        String docker = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
        Int    cpu             = 8
        Int    memory_gb       = 32
        String disk_type       = "SSD"   # see note in the reply — NOT HDD
        Int    disk_multiplier = 3
        Int    disk_pad_gb     = 50
        Int    preemptible     = 0
    }

    scatter (i in range(length(chroms))) {
        call filter_chrom {
            input:
                chrom           = chroms[i],
                pgen            = pgen_files[i],
                pvar            = pvar_files[i],
                psam            = psam_files[i],
                keeps           = keep_files,
                plink2          = plink2,
                maf             = maf,
                geno            = geno,
                hwe             = hwe,
                docker          = docker,
                cpu             = cpu,
                memory_gb       = memory_gb,
                disk_type       = disk_type,
                disk_multiplier = disk_multiplier,
                disk_pad_gb     = disk_pad_gb,
                preemptible     = preemptible
        }
    }

    output {
        Array[Array[File]] panel_pgen = filter_chrom.out_pgen
        Array[Array[File]] panel_pvar = filter_chrom.out_pvar
        Array[Array[File]] panel_psam = filter_chrom.out_psam
        Array[Array[File]] plink_logs = filter_chrom.out_logs
        Array[File]        manifests  = filter_chrom.manifest
    }
}

task filter_chrom {
    input {
        String      chrom
        File        pgen
        File        pvar
        File        psam
        Array[File] keeps
        File        plink2
        Float       maf
        Float       geno
        String      hwe
        String      docker
        Int         cpu
        Int         memory_gb
        String      disk_type
        Int         disk_multiplier
        Int         disk_pad_gb
        Int         preemptible
    }

    Int input_gb     = ceil(size(pgen, "GB") + size(pvar, "GB") + size(psam, "GB"))
    Int disk_gb      = input_gb * disk_multiplier + disk_pad_gb
    Int plink_mem_mb = memory_gb * 1024 - 4096

    command <<<
        set -euo pipefail

        echo "=== environment ==="
        echo "nproc:      $(nproc)"
        echo "mem:        $(awk '/MemTotal/{print $2/1048576" GB"}' /proc/meminfo)"
        echo "cwd device: $(df -h . | tail -1)"
        echo "rotational: $(cat /sys/block/sd*/queue/rotational 2>/dev/null | tr '\n' ' ' || echo unknown)"
        echo "start:      $(date -u +%FT%TZ)"

        cp "~{plink2}" ./plink2
        chmod 750 ./plink2
        ./plink2 --version

        # Prefix reconstruction
        # co-locates and --pfile resolves all three.
        PGEN_PATH="~{pgen}"
        PFX="${PGEN_PATH%.pgen}"
        ls -lL "${PFX}.pgen" "${PFX}.pvar" "${PFX}.psam"

        mkdir -p out

        MANIFEST="manifest.chr~{chrom}.tsv"
        printf 'chrom\tpanel\tn_keep\tn_samples_out\tn_variants_out\tn_dup_ids\tbytes_out\tseconds\n' > "$MANIFEST"

        # Loop panels. set -e above means panel 3 failing aborts the task
        # instead of leaving five good outputs and rc=0.
        while read -r K; do
            PANEL="$(basename "$K" .keep)"
            OUT="out/${PANEL}.chr~{chrom}"

            echo "=== panel ${PANEL} start $(date -u +%FT%TZ) ==="
            T0=$(date +%s)

            ./plink2 \
                --pfile "$PFX" \
                --keep "$K" \
                --snps-only --max-alleles 2 \
                --geno ~{geno} \
                --maf ~{maf} \
                --hwe ~{hwe} \
                --make-pgen \
                --threads ~{cpu} \
                --memory ~{plink_mem_mb} \
                --out "$OUT"

            T1=$(date +%s)
            echo "=== panel ${PANEL} done  $(date -u +%FT%TZ) ==="

            # awk/wc rather than grep -c: grep exits 1 on zero matches and
            # would trip pipefail.
            NKEEP=$(awk '!/^#/ && NF' "$K" | wc -l)
            NSAMP=$(awk '!/^#/' "${OUT}.psam" | wc -l)
            NVAR=$(awk '!/^#/'  "${OUT}.pvar" | wc -l)

            # ID format is chrom:pos:ref. Uniqueness is an assumption, not a
            # guarantee — check it here rather than discovering it at --pmerge-list.
            NDUP=$(awk '!/^#/{print $3}' "${OUT}.pvar" | sort | uniq -d | wc -l)

            BYTES=$(du -cb "${OUT}".pgen "${OUT}".pvar "${OUT}".psam | tail -1 | cut -f1)

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "~{chrom}" "$PANEL" "$NKEEP" "$NSAMP" "$NVAR" "$NDUP" "$BYTES" "$((T1-T0))" \
                >> "$MANIFEST"
        done < ~{write_lines(keeps)}

        echo "=== manifest ==="
        cat "$MANIFEST"
        echo "end: $(date -u +%FT%TZ)"
    >>>

    runtime {
        docker:      docker
        cpu:         cpu
        memory:      "~{memory_gb} GB"
        disks:       "local-disk ~{disk_gb} ~{disk_type}"
        preemptible: preemptible
    }

    output {
        Array[File] out_pgen = glob("out/*.pgen")
        Array[File] out_pvar = glob("out/*.pvar")
        Array[File] out_psam = glob("out/*.psam")
        Array[File] out_logs = glob("out/*.log")
        File        manifest = "manifest.chr~{chrom}.tsv"
    }
}
