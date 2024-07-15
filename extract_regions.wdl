version 1.0

workflow ExtractRegions {
    input {
        File input_fasta
        File regions_file
        Int extract_memory = 2  # Default memory in GB for the Extract task
        Int extract_cpu = 1     # Default number of CPUs for the Extract task
        String disk_size = "50G"   # Default disk size for the task
        String disk_type = "pd-ssd" # Default disk type for the task
    }

    meta {
        author: "Shadi Zaheri"
        email: "szaheri@broadinstitute.org"
        version: "1.0"
        description: "Workflow to extract specific regions or sequences from a FASTA or FASTQ file based on a list of regions specified in a text file."
    }

    call Extract {
        input:
            input_fasta = input_fasta,
            regions_file = regions_file,
            memory = extract_memory,
            cpu = extract_cpu,
            disk_size = disk_size,
            disk_type = disk_type
    }

    output {
        File extracted_fasta = Extract.extracted_fasta
    }
}

task Extract {
    input {
        File input_fasta
        File regions_file
        Int memory
        Int cpu
        String disk_size
        String disk_type
    }

    command {
        # Extract specific regions from the input FASTA file
        samtools faidx --region-file ${regions_file} ${input_fasta} > extracted_regions.fasta
    }

    output {
        File extracted_fasta = "extracted_regions.fasta"
    }

    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        memory: "${memory}G"
        cpu: "${cpu}"
        disks: "local-disk ${disk_size} ${disk_type}"
    }
}
