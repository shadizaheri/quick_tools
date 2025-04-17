version 1.0

workflow BamMd5Workflow {
  input {
    File input_bam
    String md5_memory = "2G"
    Int md5_cpu = 1
    String md5_disk = "local-disk 20 HDD"
    Int md5_preemptible = 1
  }

  call ComputeMd5 {
    input:
      bam_file = input_bam,
      memory = md5_memory,
      cpu = md5_cpu,
      disk = md5_disk,
      preemptible = md5_preemptible
  }

  output {
    File md5_output = ComputeMd5.md5_file
  }
}

task ComputeMd5 {
  input {
    File bam_file
    String memory
    Int cpu
    String disk
    Int preemptible
  }

  command {
    md5sum ~{bam_file} > md5.txt
  }

  output {
    File md5_file = "md5.txt"
  }

  runtime {
    docker: "ubuntu:20.04"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

  meta {
    author: "Shadi Zaheri"
  }
}

