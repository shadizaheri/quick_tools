version 1.0

workflow fastqc_workflow {
    input {
        File fastq1                     # Required FASTQ file for single-end reads
        File? fastq2                    # Optional second FASTQ file for paired-end reads
        Int disk_space = 10             # Default disk space in GB
        Int memory = 4                  # Default memory in GB
        Int num_threads = 2             # Default number of CPU threads
    }

    call fastqc {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            disk_space = disk_space,
            memory = memory,
            num_threads = num_threads
    }

    output {
        File fastq1_report = fastqc.fastq1_report
        File? fastq2_report = fastqc.fastq2_report
        Array[File] zips = fastqc.zips
    }
}

task fastqc {
    input {
        File fastq1                  # Required FASTQ file
        File? fastq2                 # Optional FASTQ file for paired-end reads
        Int disk_space               # Disk space in GB
        Int memory                   # Memory in GB
        Int num_threads              # Number of CPU threads
    }

    command <<<
        set -euo pipefail

        # Create output directory
        mkdir -p fastqc_output

        # Run FastQC on the first FASTQ file
        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq1}

        # Optionally run FastQC on the second FASTQ file
        if [ -n "~{fastq2}" ]; then
            fastqc \
                --outdir fastqc_output \
                --threads ${num_threads} \
                ~{fastq2}
        fi
    >>>

    output {
        File fastq1_report = glob("fastqc_output/*${basename(fastq1)}*_fastqc.html")[0]
        File? fastq2_report = select_first([glob("fastqc_output/*${basename(fastq2)}*_fastqc.html")[0], ""])
        Array[File] zips = glob("fastqc_output/*.zip")
    }

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Shadi Zaheri"
        description: "Run FastQC to perform quality control on FASTQ files, with optional paired-end processing."
    }
}
