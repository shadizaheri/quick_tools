version 1.0

workflow InsertSizeMetrics {
    input {
        File cram_file
        File cram_index
        File reference_fasta

        # User-defined runtime settings (from Terra input JSON)
        Int? disk_size_override
        Int? memory_override
        Int? cpu_override
        Int? preemptible_override
    }

    call SamtoolsInsertSize {
        input:
            cram=cram_file,
            reference_fasta=reference_fasta,
            disk_size=disk_size_override,
            memory=memory_override,
            cpu=cpu_override,
            preemptible=preemptible_override
    }

    call PicardInsertSize {
        input:
            cram=cram_file,
            cram_index=cram_index,
            reference_fasta=reference_fasta,
            disk_size=disk_size_override,
            memory=memory_override,
            cpu=cpu_override,
            preemptible=preemptible_override
    }

    output {
        File samtools_output = SamtoolsInsertSize.stats_output
        File picard_metrics = PicardInsertSize.metrics_output
        File picard_histogram = PicardInsertSize.histogram_output
    }
}

task SamtoolsInsertSize {
    input {
        File cram
        File reference_fasta
        Int? disk_size
        Int? memory
        Int? cpu
        Int? preemptible
    }

    command {
        samtools stats -@ ~{if defined(cpu) then cpu else 4} --reference ~{reference_fasta} ~{cram} | grep "insert size" > insert_size_stats.txt
    }

    output {
        File stats_output = "insert_size_stats.txt"
    }

    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        memory: "~{if defined(memory) then memory else 8}G"
        cpu: "~{if defined(cpu) then cpu else 4}"
        disks: "local-disk ~{if defined(disk_size) then disk_size else 100} HDD"
        preemptible: "~{if defined(preemptible) then preemptible else 0}"
    }
}

task PicardInsertSize {
    input {
        File cram
        File cram_index
        File reference_fasta
        Int? disk_size
        Int? memory
        Int? cpu
        Int? preemptible
    }

    command {
        gatk CollectInsertSizeMetrics \
            -I ~{cram} \
            -O insert_size_metrics.txt \
            -H insert_size_histogram.pdf \
            -M 0.5
    }

    output {
        File metrics_output = "insert_size_metrics.txt"
        File histogram_output = "insert_size_histogram.pdf"
    }

    runtime {
        docker: "broadinstitute/gatk:4.6.1.0"
        memory: "~{if defined(memory) then memory else 8}G"
        cpu: "~{if defined(cpu) then cpu else 4}"
        disks: "local-disk ~{if defined(disk_size) then disk_size else 100} HDD"
        preemptible: "~{if defined(preemptible) then preemptible else 0}"
    }
}
