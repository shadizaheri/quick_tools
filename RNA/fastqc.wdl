version 1.0
workflow fastqc_workflow {
    input {
        File fastq1                     # First FASTQ file (required)
        File fastq2                     # Second FASTQ file (required for paired-end reads)
        Int disk_space = 50             # Default disk space in GB
        Int memory = 8                  # Default memory in GB
        Int num_threads = 4             # Default number of CPU threads
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
        File report1 = fastqc.report1         # HTML report for fastq1
        File report2 = fastqc.report2         # HTML report for fastq2
        #File zip1 = fastqc.zip1               
        #File zip2 = fastqc.zip2              
    }
}

task fastqc {
    input {
        File fastq1                     # First FASTQ file
        File fastq2                     # Second FASTQ file
        Int disk_space                  # Disk space in GB
        Int memory                      # Memory in GB
        Int num_threads                 # Number of CPU threads
    }

    command <<<
        set -euo pipefail

        # Create output directory
        mkdir -p fastqc_output

        echo "Running FastQC with the following inputs:"
        echo "FASTQ1: ${fastq1}"
        echo "FASTQ2: ${fastq2}"
        echo "Threads: ${num_threads}"

        # Run FastQC for each file
        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq1}

        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq2}
    >>>

    output {
        # Match FastQC's default output naming
        File report1 = glob("fastqc_output/*_fastqc.html")[0]  # First FASTQ HTML report
        File report2 = glob("fastqc_output/*_fastqc.html")[1]  # Second FASTQ HTML report
        #File zip1 = glob("fastqc_output/*_fastqc.zip")[0]      # First FASTQ ZIP file
        #File zip2 = glob("fastqc_output/*_fastqc.zip")[1]      # Second FASTQ ZIP file
    }

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Shadi Zaheri"
        description: "Run FastQC to perform quality control on FASTQ files"
    }
}
