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

struct DataTypeParameters {
    Int num_shards
    String map_preset
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
            runtime_attrs = select_first([runtime_attrs, default_runtime_attrs()])
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
        RuntimeAttr runtime_attrs
    }

    command {
        java -jar picard.jar FastqToSam \
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

    runtime {
        docker: runtime_attrs.docker
        memory: "${runtime_attrs.mem_gb} GB"
        cpu: runtime_attrs.cpu_cores
        disks: "local-disk ${runtime_attrs.disk_gb} HDD"
        bootDiskSizeGb: runtime_attrs.boot_disk_gb
        preemptible: runtime_attrs.preemptible_tries
        maxRetries: runtime_attrs.max_retries
    }
}

# Function to provide default runtime attributes
function default_runtime_attrs() {
    return {
        mem_gb: 4.0,
        cpu_cores: 2,
        disk_gb: 50,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 2,
        docker: "broadinstitute/gatk:latest"
    }
}
