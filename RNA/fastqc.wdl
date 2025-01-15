version 1.0

workflow fastqc_workflow {
    input {
        File fastq1                     
        File? fastq2                    
        String output_prefix            
        Int disk_space = 10             # Default disk space in GB
        Int memory = 4                  # Default memory in GB
        Int num_threads = 2             # Default number of CPU threads
    }

    call fastqc {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            output_prefix = output_prefix,
            disk_space = disk_space,
            memory = memory,
            num_threads = num_threads
    }

    output {
        File fastq1_report = fastqc.fastq1_report
        File? fastq2_report = fastqc.fastq2_report
    }
}

task fastqc {
    input {
        File fastq1                  
        File? fastq2                 # Optional: Second FASTQ file for paired-end reads
        String output_prefix         
        Int disk_space               # GB
        Int memory                   # GB
        Int num_threads              # Number of CPU threads
    }

    command {
        set -euo pipefail

        # Create output directory
        mkdir -p fastqc_output

        # Run FastQC on fastq1
        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq1}
        mv fastqc_output/*.html fastqc_output/${output_prefix}_fastq1.html
        mv fastqc_output/*.zip fastqc_output/${output_prefix}_fastq1.zip

        # Optionally run FastQC on fastq2 if provided
        if [ -n "${fastq2}" ]; then
            fastqc \
                --outdir fastqc_output \
                --threads ${num_threads} \
                ${fastq2}
            mv fastqc_output/*.html fastqc_output/${output_prefix}_fastq2.html
            mv fastqc_output/*.zip fastqc_output/${output_prefix}_fastq2.zip
        fi
    }

    output {
        File fastq1_report = glob("fastqc_output/*_fastq1.html")[0]
        File? fastq2_report = select_first([glob("fastqc_output/*_fastq2.html")[0], ""])
    }

    runtime {
        docker: "biocontainers/fastqc:v0.11.9_cv8"          
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
    }

    meta {
        author: "Shadi Zaheri"
        description: "Run FastQC to perform quality control on FASTQ files, with separate outputs for fastq1 and fastq2."
    }
}
