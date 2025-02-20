version 1.0
workflow InsertSizeMetrics {
    input {
        File cram_file
        File cram_index
        File reference_fasta
        File reference_fasta_index
    }

    # Get CRAM file size to estimate disk space
    call GetFileSize {
        input: input_file = cram_file
    }

    call SamtoolsInsertSize {
        input:
            cram=cram_file,
            reference_fasta=reference_fasta,
            disk_size=GetFileSize.disk_space_needed
    }

    call PicardInsertSize {
        input:
            cram=cram_file,
            cram_index=cram_index,
            reference_fasta=reference_fasta,
            disk_size=GetFileSize.disk_space_needed
    }

    output {
        File samtools_output = SamtoolsInsertSize.stats_output
        File picard_metrics = PicardInsertSize.metrics_output
        File picard_histogram = PicardInsertSize.histogram_output
    }
}

task GetFileSize {
    input {
        File input_file
    }

    command <<<
        # Get the file size in GB
        FILE_SIZE_GB=$(ls -l ~{input_file} | awk '{print $5 / (1024*1024*1024)}')
        
        # Ensure at least 50GB disk or 3x file size
        DISK_NEEDED=$(echo "$FILE_SIZE_GB * 3" | bc)
        if (( $(echo "$DISK_NEEDED < 50" | bc -l) )); then
            DISK_NEEDED=50
        fi
        echo ${DISK_NEEDED%.*} > disk_size.txt
    >>>

    output {
        Int disk_space_needed = read_int("disk_size.txt")
    }

    runtime {
        docker: "ubuntu:latest"
        memory: "2G"
        cpu: 1
    }
}

task SamtoolsInsertSize {
    input {
        File cram
        File reference_fasta
        Int? disk_size  # Optional, can be overridden in Terra
    }

    command {
        samtools stats -@ 4 --reference ~{reference_fasta} ~{cram} | grep "insert size" > insert_size_stats.txt
    }

    output {
        File stats_output = "insert_size_stats.txt"
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.17--h6899075_1"
        memory: "8G"
        cpu: 4
        disks: "local-disk ~{if defined(disk_size) then disk_size else 100} HDD"  # Fixed syntax
    }
}

task PicardInsertSize {
    input {
        File cram
        File cram_index
        File reference_fasta
        Int? disk_size  # Optional, can be overridden in Terra
    }

    command {
        picard CollectInsertSizeMetrics \
            I=~{cram} \
            O=insert_size_metrics.txt \
            H=insert_size_histogram.pdf \
            M=0.5
    }

    output {
        File metrics_output = "insert_size_metrics.txt"
        File histogram_output = "insert_size_histogram.pdf"
    }

    runtime {
        docker: "broadinstitute/picard:latest"
        memory: "8G"
        cpu: 4
        disks: "local-disk ~{if defined(disk_size) then disk_size else 100} HDD"  # Fixed syntax
    }
}
