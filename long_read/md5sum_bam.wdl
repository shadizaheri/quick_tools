version 1.0 

task ComputeMd5 {
  input {
    File bam_file
    Int cpu_cores = 1
    String memory_gb = "2 GB"
    String disk_space = "10 GB"
  }

  command {
    md5sum ~{bam_file} > md5.txt
  }

  output {
    File md5_file = "md5.txt"
  }

  runtime {
    cpu: cpu_cores
    memory: memory_gb
    disks: "local-disk " + disk_space
    docker: "ubuntu:20.04"
  }
}

workflow BamMd5Workflow {
  input {
    File input_bam
    Int cpu_cores = 1
    String memory_gb = "2 GB"
    String disk_space = "10 GB"
  }

  call ComputeMd5 {
    input:
      bam_file = input_bam,
      cpu_cores = cpu_cores,
      memory_gb = memory_gb,
      disk_space = disk_space
  }

  output {
    File md5_output = ComputeMd5.md5_file
  }
  meta {
    author: "Shadi Zaheri"
  }
}

