version 1.0

workflow fastqc_workflow {
    input {
        File fastq1                     # First FASTQ file (required)
        File fastq2                     # Second FASTQ file (required for paired-end reads)
        String output_prefix1           # Prefix for output of fastq1
        String output_prefix2           # Prefix for output of fastq2
        Int disk_space = 50             # Default disk space in GB (adjustable)
        Int memory = 8                  # Default memory in GB (adjustable)
        Int num_threads = 4             # Default number of CPU threads (adjustable)
    }

    call fastqc {
        input:
            fastq1 = fastq1,
            fastq2 = fastq2,
            output_prefix1 = output_prefix1,
            output_prefix2 = output_prefix2,
            disk_space = disk_space,
            memory = memory,
            num_threads = num_threads
    }

    output {
        File report1 = fastqc.report1         
        File report2 = fastqc.report2         
        #File zip1 = fastqc.zip1               # ZIP output for fastq1
        #File zip2 = fastqc.zip2               # ZIP output for fastq2
    }
}

task fastqc {
    input {
        File fastq1                     # First FASTQ file
        File fastq2                     # Second FASTQ file
        String output_prefix1           
        String output_prefix2           
        Int disk_space                  # Disk space in GB
        Int memory                      # Memory in GB
        Int num_threads                 # Number of CPU threads
    }

    command <<<
        set -euo pipefail

        mkdir -p fastqc_output

        echo "Running FastQC with the following inputs:"
        echo "FASTQ1: ${fastq1}"
        echo "FASTQ2: ${fastq2}"
        echo "Threads: ${num_threads}"
        echo "Output directory: fastqc_output"

        # Run FastQC for each file separately
        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq1}

        fastqc \
            --outdir fastqc_output \
            --threads ${num_threads} \
            ${fastq2}

        for report in fastqc_output/$(basename ${fastq1%.*})_fastqc.html; do
            mv $report fastqc_output/${output_prefix1}.html
        done

        #for zip in fastqc_output/$(basename ${fastq1%.*})_fastqc.zip; do
        #    mv $zip fastqc_output/${output_prefix1}.zip
        #done

        for report in fastqc_output/$(basename ${fastq2%.*})_fastqc.html; do
            mv $report fastqc_output/${output_prefix2}.html
        done

        #for zip in fastqc_output/$(basename ${fastq2%.*})_fastqc.zip; do
        #    mv $zip fastqc_output/${output_prefix2}.zip
        #done
    >>>

    output {
        File report1 = "fastqc_output/${output_prefix1}.html"  
        File report2 = "fastqc_output/${output_prefix2}.html"  
        #File zip1 = "fastqc_output/${output_prefix1}.zip"      
        #File zip2 = "fastqc_output/${output_prefix2}.zip"      
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
