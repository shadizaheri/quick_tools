version 1.0

workflow fastqc_workflow {
    input {
        File fastq1                     # First FASTQ file (required)
        File fastq2                     # Second FASTQ file (required for paired-end reads)
        Int disk_space = 50             # Default disk space in GB
        Int memory = 8                  # Default memory in GB
        Int num_threads = 4             # Default number of CPU threads
    }

    # Run FastQC for fastq1
    call run_fastqc {
        input:
            fastq = fastq1,
            disk_space = disk_space,
            memory = memory,
            num_threads = num_threads
    }

    # Run FastQC for fastq2
    call run_fastqc as run_fastqc_2 {
        input:
            fastq = fastq2,
            disk_space = disk_space,
            memory = memory,
            num_threads = num_threads
    }

    # Outputs
    output {
        File fastqc1_html = run_fastqc.fastqc_html
        File fastqc2_html = run_fastqc_2.fastqc_html
    }
}

task run_fastqc {
    # Inputs
    input {
        File fastq
        Int memory                      # Memory in GB
        Int disk_space                  # Disk space in GB
        Int num_threads                 # Number of CPU threads
    }

    # Command to run FastQC
    command {
        fastqc ${fastq} --threads ${num_threads}
    }

    # Runtime settings
    runtime {
        memory: "${memory} GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: num_threads
        docker: "biocontainers/fastqc:v0.11.9_cv8"  # Use FastQC Docker container
    }

    # Outputs
    output {
        File fastqc_zip = "${basename(fastq)}_fastqc.zip"
        File fastqc_html = "${basename(fastq)}_fastqc.html"
    }
    meta {
      author: "Shadi Zaheri"
      description: "Run FastQC to perform quality control on FASTQ files"
    }
}
