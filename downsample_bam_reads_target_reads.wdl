version 1.0

workflow DownsampleAndIndexBam {
    input {
        File input_bam
        Int target_reads
        String sample_name
        Int? downsample_memory = 2  # GB
        Int? downsample_cpu = 1      
        Int? index_memory = 2        
        Int? index_cpu = 1           
        String? disk_space = "local-disk 50G"  # Unified disk setting
    }

    meta {
        author: "Shadi Zaheri"
        email: "szaheri@broadinstitute.org"
        version: "1.2"
        description: "Workflow to downsample a BAM file to an exact number of reads and index the resulting BAM file."
    }

    call CountReads { input: bam = input_bam }

    call CalculateFraction {
        input:
            total_reads = CountReads.total_reads,
            target_reads = target_reads
    }

    call Downsample {
        input:
            input_bam = input_bam,
            fraction = CalculateFraction.fraction,
            sample_name = sample_name,
            memory = downsample_memory,
            cpu = downsample_cpu,
            disk_space = disk_space
    }

    call IndexBam {
        input:
            bam = Downsample.downsampled_bam,
            memory = index_memory,
            cpu = index_cpu,
            disk_space = disk_space
    }

    output {
        File downsampled_bam = Downsample.downsampled_bam
        File downsampled_bam_index = IndexBam.bam_index
    }
}

task CountReads {
    input {
        File bam
    }
    command {
        samtools view -c ${bam} > total_reads.txt
    }
    output {
        Int total_reads = read_int("total_reads.txt")
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
    }
}

task CalculateFraction {
    input {
        Int total_reads
        Int target_reads
    }
    command {
        echo "scale=6; ${target_reads} / ${total_reads}" | bc > fraction.txt
    }
    output {
        Float fraction = read_float("fraction.txt")
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

task Downsample {
    input {
        File input_bam
        Float fraction
        String sample_name
        Int? memory
        Int? cpu
        String? disk_space
    }

    command {
        samtools view -@ ${cpu} -s ${fraction} -b ${input_bam} > ${sample_name}_downsampled.bam
    }

    output {
        File downsampled_bam = "${sample_name}_downsampled.bam"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
        memory: select_first([memory, "8G"])
        cpu: select_first([cpu, 4])
        disks: select_first([disk_space, "local-disk 50G"])
    }
}

task IndexBam {
    input {
        File bam
        Int? memory
        Int? cpu
        String? disk_space
    }

    command {
        samtools index -@ ${cpu} ${bam} ${bam}.bai
    }

    output {
        File bam_index = "${bam}.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
        memory: select_first([memory, "8G"])
        cpu: select_first([cpu, 4])
        disks: select_first([disk_space, "local-disk 50G"])
    }
}