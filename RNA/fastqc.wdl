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
        Array[File] reports = fastqc.reports
        Array[File] zips = fastqc.zips
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

        # Run FastQC
        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq1} ${"-" + fastq2}

        # Rename outputs with the specified prefix
        for report in fastqc_output/*.html; do
            mv $report fastqc_output/${output_prefix}_$(basename $report)
        done
        for zip in fastqc_output/*.zip; do
            mv $zip fastqc_output/${output_prefix}_$(basename $zip)
        done
    }

    output {
        Array[File] reports = glob("fastqc_output/*.html")  
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
        description: "Run FastQC to perform quality control on FASTQ files"
    }
}

