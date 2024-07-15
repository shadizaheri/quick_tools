version 1.0

workflow DownsampleAndIndexBam {
    input {
        File input_bam
        String fraction
        Int downsample_memory = 2  # Default memory in GB for the Downsample task
        Int downsample_cpu = 1     # Default number of CPUs for the Downsample task
        Int index_memory = 2       # Default memory in GB for the Index task
        Int index_cpu = 1          # Default number of CPUs for the Index task
        String disk_size = "50G"   # Default disk size for both tasks
        String disk_type = "pd-ssd" # Default disk type for both tasks
    }

    meta {
        author: "Shadi Zaheri"
        email: "szaheri@broadinstitute.org"
        version: "1.0"
        description: "Workflow to downsample a BAM file based on reads and index the resulting BAM file."
    }

    call Downsample {
        input:
            input_bam = input_bam,
            fraction = fraction,
            memory = downsample_memory,
            cpu = downsample_cpu,
            disk_size = disk_size,
            disk_type = disk_type
    }

    call IndexBam {
        input:
            bam = Downsample.downsampled_bam,
            memory = index_memory,
            cpu = index_cpu,
            disk_size = disk_size,
            disk_type = disk_type
    }

    output {
        File downsampled_bam = Downsample.downsampled_bam
        File downsampled_bam_index = IndexBam.bam_index
    }
}

task Downsample {
    input {
        File input_bam
        String fraction
        Int memory
        Int cpu
        String disk_size
        String disk_type
    }

    command {
        # Downsample the BAM file based on the specified fraction of reads
        samtools view -s ${fraction} -b ${input_bam} > downsampled.bam
    }

    output {
        File downsampled_bam = "downsampled.bam"
    }

    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        memory: "${memory}G"
        cpu: "${cpu}"
        disks: "local-disk ${disk_size} ${disk_type}"
    }
}

task IndexBam {
    input {
        File bam
        Int memory
        Int cpu
        String disk_size
        String disk_type
    }

    command {
        # Index the downsampled BAM file
        samtools index ${bam}
    }

    output {
        File bam_index = "${bam}.bai"
    }

    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        memory: "${memory}G"
        cpu: "${cpu}"
        disks: "local-disk ${disk_size} ${disk_type}"
    }
}
