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

    # Count soft‑clipped reads
    soft=$(samtools view ${bam_file} | grep -cE '[0-9]+S')

    # Count hard‑clipped reads
    hard=$(samtools view ${bam_file} | grep -cE '[0-9]+H')

    # Write header line
    echo -e "file\tsoft_clipped_reads\thard_clipped_reads" > counts.tsv
    # Use shell variables for soft and hard; WDL variable for bam_file
    echo -e "${bam_file}\t\$soft\t\$hard" >> counts.tsv
  }

  output {
    File counts = "counts.tsv"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/gtex_v8_star_2.7.10a_custom"
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
