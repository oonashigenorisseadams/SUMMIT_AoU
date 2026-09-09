version 1.0

task probe_cdr {
  input {
    File plink2_bin
    File pgen
    File pvar
    File psam
  }

command <<<
    set -uo pipefail

    echo "=== staged files ==="
    ls -lh ~{pgen} ~{pvar} ~{psam}

    echo "=== plink2 ==="
    chmod 750 ~{plink2_bin}
    ~{plink2_bin} --version

    echo "=== can plink2 open the trio? ==="
    ~{plink2_bin} --pfile ~{sub(pgen, "\.pgen$", "")} --validate
  >>>

  output {
    File stdout_log = stdout()
    File stderr_log = stderr()
  }

  runtime {
    docker: "ubuntu:22.04"
    cpu: "2"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
  }
}

workflow probe_cdr_access {
  input {
    File plink2_bin
    File pgen
    File pvar
    File psam
  }
  call probe_cdr {
    input: plink2_bin = plink2_bin, pgen = pgen, pvar = pvar, psam = psam
  }
  output {
    File out = probe_cdr.stdout_log
    File err = probe_cdr.stderr_log
  }
}
