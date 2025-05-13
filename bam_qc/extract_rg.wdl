version 1.0

workflow SplitBAMByReadGroupOnly {
  input {
    String sample_name
    File input_bam
    File input_bam_index
    String disk_space
  }

  call ExtractReadGroups {
    input: bam = input_bam
  }

  scatter (rg_id in ExtractReadGroups.read_group_ids) {
    call ExtractRG_BAM {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        rg_id = rg_id,
        sample_name = sample_name
    }
  }

  output {
    Array[File] per_rg_bams = ExtractRG_BAM.output_bam
    Array[File] per_rg_bam_indexes = ExtractRG_BAM.output_bam_index
  }
}

task ExtractReadGroups {
  input {
    File bam
  }

  command {
    samtools view -H ~{bam} | grep '^@RG' | sed 's/.*ID:\([^\t]*\).*/\1/' > rg_ids.txt
  }

  output {
    Array[String] read_group_ids = read_lines("rg_ids.txt")
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/samtools:1.18"
    cpu: 1
    memory: "1G"
    disks: "local-disk 10 HDD"
  }
}

task ExtractRG_BAM {
  input {
    File input_bam
    File input_bam_index
    String rg_id
    String sample_name
  }

  command {
    set -euo pipefail
    samtools view -@ 2 -b -r ~{rg_id} ~{input_bam} > ~{sample_name}.~{rg_id}.bam
    samtools index ~{sample_name}.~{rg_id}.bam
  }

  output {
    File output_bam = "~{sample_name}.~{rg_id}.bam"
    File output_bam_index = "~{sample_name}.~{rg_id}.bam.bai"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-gcr-public/samtools:1.18"
    cpu: 2
    memory: "4G"
    disks: "local-disk 20 HDD"
  }
}
