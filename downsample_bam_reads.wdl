version 1.0

workflow DownsampleAndIndexBam {
    input {
        File input_bam
        Float fraction
        String sample_name
        Int downsample_memory = 2  # Default memory in GB
        Int downsample_cpu = 1      # Default number of CPUs
        Int index_memory = 2        # Default memory for indexing
        Int index_cpu = 1           # Default CPUs for indexing
        String disk_size = "50G"    # Default disk size
        String disk_type = "pd-ssd" # Default disk type
    }

    meta {
        author: "Shadi Zaheri"
        email: "szaheri@broadinstitute.org"
        version: "1.1"
        description: "Workflow to downsample a BAM file and index the resulting BAM file."
    }

    call Downsample {
        input:
            input_bam = input_bam,
            fraction = fraction,
            sample_name = sample_name,
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
        Float fraction
        String sample_name
        Int memory
        Int cpu
        String disk_size
        String disk_type
    }

    command {
        samtools view -@ ${cpu} -s ${fraction} -b ${input_bam} > ${sample_name}_downsampled.bam
    }

    output {
        File downsampled_bam = "${sample_name}_downsampled.bam"
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.16--h9aed4be_0"
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
        samtools index -@ ${cpu} ${bam} ${bam}.bai
    }

    output {
        File bam_index = "${bam}.bai"
    }

    runtime {
        docker: "quay.io/biocontainers/samtools:1.16--h9aed4be_0"
        memory: "${memory}G"
        cpu: "${cpu}"
        disks: "local-disk ${disk_size} ${disk_type}"
    }
}