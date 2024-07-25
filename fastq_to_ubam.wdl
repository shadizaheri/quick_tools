version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

workflow FastqToSamWorkflow {
    input {
        String sample_id
        File fastq_file
        String read_group_name
        String sample_name
        String library_name
        String platform
        RuntimeAttr? runtime_attrs
    }

    call FastqToSam {
        input:
            fastq_file = fastq_file,
            read_group_name = read_group_name,
            sample_name = sample_name,
            library_name = library_name,
            platform = platform,
            runtime_attrs = runtime_attrs
    }

    output {
        File ubam_file = FastqToSam.ubam_file
    }
}

task FastqToSam {
    input {
        File fastq_file
        String read_group_name
        String sample_name
        String library_name
        String platform
        RuntimeAttr? runtime_attrs
    }
    Int disk_size = 32

    command {
        java -jar /usr/picard/picard.jar FastqToSam \
            FASTQ=${fastq_file} \
            OUTPUT=${sample_name}.unaligned.bam \
            READ_GROUP_NAME=${read_group_name} \
            SAMPLE_NAME=${sample_name} \
            LIBRARY_NAME=${library_name} \
            PLATFORM=${platform}
    }

    output {
        File ubam_file = "${sample_name}.unaligned.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:                "us.gcr.io/broad-dsp-lrma/picard:lrfp-clr"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attrs, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}