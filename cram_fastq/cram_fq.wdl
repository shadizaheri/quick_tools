version 1.0

workflow cram_to_fastq_workflow {
  input {
    File cram_file
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    String? cram_to_bam_memory = "4G"
    Int? cram_to_bam_cpu = 2
    String? cram_to_bam_disk = "local-disk 20 HDD"
    Int? cram_to_bam_preemptible = 1

    String? bam_to_ubam_memory = "4G"
    Int? bam_to_ubam_cpu = 2
    String? bam_to_ubam_disk = "local-disk 20 HDD"
    Int? bam_to_ubam_preemptible = 1

    String? ubam_to_fastq_memory = "4G"
    Int? ubam_to_fastq_cpu = 2
    String? ubam_to_fastq_disk = "local-disk 20 HDD"
    Int? ubam_to_fastq_preemptible = 1

    String prefix
  }

  call CramToBam {
    input:
      cram = cram_file,
      reference = reference_fasta,
      memory = select_first([cram_to_bam_memory, "4G"]),
      cpu = select_first([cram_to_bam_cpu, 2]),
      disk = select_first([cram_to_bam_disk, "local-disk 20 HDD"]),
      preemptible = select_first([cram_to_bam_preemptible, 1])
  }

  call BamToUbam {
    input:
      bam = CramToBam.bam,
      reference = reference_fasta,
      memory = select_first([bam_to_ubam_memory, "4G"]),
      cpu = select_first([bam_to_ubam_cpu, 2]),
      disk = select_first([bam_to_ubam_disk, "local-disk 20 HDD"]),
      preemptible = select_first([bam_to_ubam_preemptible, 1])
  }

  call UbamToFastq {
    input:
      ubam = BamToUbam.ubam,
      prefix = prefix,
      memory = select_first([ubam_to_fastq_memory, "4G"]),
      cpu = select_first([ubam_to_fastq_cpu, 2]),
      disk = select_first([ubam_to_fastq_disk, "local-disk 20 HDD"]),
      preemptible = select_first([ubam_to_fastq_preemptible, 1])
  }

  output {
    File read1_fastq = UbamToFastq.read1_fastq
    File read2_fastq = UbamToFastq.read2_fastq
  }
}

task CramToBam {
  input {
    File cram
    File reference
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command {
    samtools view -b -T ~{reference} ~{cram} > aligned.bam
  }
  output {
    File bam = "aligned.bam"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }
}

task BamToUbam {
  input {
    File bam
    File reference
    String memory
    Int cpu
    String disk
    Int preemptible
  }
  command {
    gatk RevertSam \
      -I ~{bam} \
      -O unmapped.bam \
      --SANITIZE true \
      --REMOVE_ALIGNMENT_INFORMATION true \
      --RESTORE_ORIGINAL_QUALITIES true
  }
  output {
    File ubam = "unmapped.bam"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }
}

task UbamToFastq {
  input {
    File ubam
    String memory
    Int cpu
    String disk
    Int preemptible
    String prefix
  }
  command {
    gatk SamToFastq \
      -I ~{ubam} \
      -F ~{prefix}_read1.fastq \
      -F2 ~{prefix}_read2.fastq
  }

  output {
    File read1_fastq = "~{prefix}_read1.fastq"
    File read2_fastq = "~{prefix}_read2.fastq"
  }
  runtime {
    docker: "us.gcr.io/broad-gatk/gatk:latest"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }
}
