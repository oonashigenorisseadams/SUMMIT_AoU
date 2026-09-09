version 1.0

## Measures Cromwell's concurrent-job-limit on the Workbench-managed engine.
## No inputs, no localization, 1 vCPU, 10 GB boot-only disk. 24 shards x 4 min
## is on the order of 2 vCPU-hours total and no egress.
##
## Read the result by sorting stamps on start_epoch: shards that start within
## seconds of each other are one batch. The batch size is the limit.
## If all 24 start together, raise n_shards to 150 and re-run before assuming
## 132 jobs will actually run in parallel.
##
## Turn call caching OFF before running, and before any re-run.

workflow concurrency_probe {
    input {
        Int    n_shards      = 24
        Int    sleep_seconds = 240
        String docker = "ubuntu@sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc"
    }

    scatter (i in range(n_shards)) {
        call ping {
            input:
                idx           = i,
                sleep_seconds = sleep_seconds,
                docker        = docker
        }
    }

    call collate { input: stamps = ping.stamp, docker = docker }

    output {
        File timeline = collate.timeline
    }
}

task ping {
    input {
        Int    idx
        Int    sleep_seconds
        String docker
    }

    command <<<
        set -euo pipefail
        S=$(date -u +%s)
        sleep ~{sleep_seconds}
        E=$(date -u +%s)
        printf '%s\t%s\t%s\t%s\n' "~{idx}" "$(hostname)" "$S" "$E" > stamp.tsv
    >>>

    runtime {
        docker:      docker
        cpu:         1
        memory:      "1 GB"
        disks:       "local-disk 10 HDD"
        preemptible: 0
    }

    output {
        File stamp = "stamp.tsv"
    }
}

task collate {
    input {
        Array[File] stamps
        String      docker
    }

    command <<<
        set -euo pipefail
        printf 'shard\thost\tstart_epoch\tend_epoch\n' > timeline.tsv
        while read -r F; do
            cat "$F" >> timeline.tsv
        done < ~{write_lines(stamps)}

        # sorted by start time, with seconds since the first shard started
        { head -1 timeline.tsv; tail -n +2 timeline.tsv | sort -k3,3n; } \
            | awk 'NR==1{print $0"\toffset_s"; next} NR==2{t0=$3} {print $0"\t"($3-t0)}'
    >>>

    runtime {
        docker:      docker
        cpu:         1
        memory:      "1 GB"
        disks:       "local-disk 10 HDD"
        preemptible: 0
    }

    output {
        File timeline = "timeline.tsv"
    }
}
