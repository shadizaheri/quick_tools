version 1.0

task count_clips {
  input {
    File bam_file
    Int cpu = 2
    String memory = "4G"
    String disks = "local-disk 50 HDD"
  }

  command {
    set -e

    soft=$(samtools view ${bam_file} | grep -cE '[0-9]+S')

    hard=$(samtools view ${bam_file} | grep -cE '[0-9]+H')

    echo -e "file\tsoft_clipped_reads\thard_clipped_reads" > counts.tsv
    echo -e "${bam_file}\t\$soft\t\$hard" >> counts.tsv
  }

  output {
    File counts = "counts.tsv"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    cpu: cpu
    memory: memory
    disks: disks
  }
}

workflow count_clips_workflow {
  input {
    Array[File] bam_files
    Int cpu = 2
    String memory = "4G"
    String disks = "local-disk 100 SSD"
  }

  scatter (bam in bam_files) {
    call count_clips {
      input:
        bam_file = bam,
        cpu       = cpu,
        memory    = memory,
        disks     = disks
    }
  }

  output {
    Array[File] counts_reports = count_clips.counts
  }
}
