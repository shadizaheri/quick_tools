version 1.0

workflow DeepVariantRunner {
  input {
    File bam
    File bai

    File ref_fasta
    File ref_fasta_fai

    String model_type
    Int threads = 4
    Int memory = 16
    String zones = "us-central1-a"

    String? haploid_contigs
    File? par_regions_bed

    Object? runtime_attr_override
  }

  call DV {
    input:
      bam = bam,
      bai = bai,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      model_type = model_type,
      threads = threads,
      memory = memory,
      zones = zones,
      haploid_contigs = haploid_contigs,
      par_regions_bed = par_regions_bed,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File VCF = DV.VCF
    File VCF_tbi = DV.VCF_tbi
    File gVCF = DV.gVCF
    File gVCF_tbi = DV.gVCF_tbi
    File visual_report_html = DV.visual_report_html
    File resouce_monitor_log = DV.resouce_monitor_log
    File? gpu_monitor_log = DV.gpu_monitor_log
    File output_dir_structure = DV.output_dir_structure
  }
}

task DV {

    parameter_meta {
        model_type: "which DV pre-trained model to use. Must be one of [PACBIO, ONT_R104] (or anything later supported after DV's 1.5.0 release)."
    }
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        String? haploid_contigs
        File? par_regions_bed

        String model_type

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    command <<<
        set -euxo pipefail

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        export MONITOR_MOUNT_POINT="/cromwell_root/"
        bash /opt/vm_local_monitoring_script.sh &> resources.log &
        job_id=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=~{model_type} \
            --ref=~{ref_fasta} \
            ~{true='--haploid_contigs ' false='' defined(haploid_contigs)} ~{select_first([haploid_contigs, ""])} \
            ~{true='--par_regions_bed ' false='' defined(par_regions_bed)} ~{select_first([par_regions_bed, ""])} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}" || cat resources.log
        if ps -p "${job_id}" > /dev/null; then kill "${job_id}"; fi

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"
    >>>

    output {

        File resouce_monitor_log = "resources.log"
        File? gpu_monitor_log = "gpu.usages.log"
        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"
    }

    #########################
    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 20
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-deepvariant:1.6.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        noAddress: true

        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones: zones
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
