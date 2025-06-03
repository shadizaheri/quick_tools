version 1.0

workflow SamtoolsStatsWorkflow {
  input {
    File input_bam_or_cram
    File reference_fasta
    String sample_name

    # Adjustable runtime settings
    Int cpu = 1
    String memory = "2G"
    String disk = "10 HDD"
  }

  call SamtoolsStatsTask {
    input:
      input_bam_or_cram = input_bam_or_cram,
      reference_fasta = reference_fasta,
      sample_name = sample_name,
      cpu = cpu,
      memory = memory,
      disk = disk
  }

  output {
    File stats_output = SamtoolsStatsTask.stats_output
  }
}

task SamtoolsStatsTask {
  input {
    File input_bam_or_cram
    File reference_fasta
    String sample_name

    Int cpu
    String memory
    String disk
  }

  command {
    set -e

    # Optional CRAM reference setup
    if [[ "${input_bam_or_cram}" == *.cram ]]; then
      export SAMTOOLS_REF="${reference_fasta}"
    fi

    samtools stats ${input_bam_or_cram} > ${sample_name}.stats
  }

  output {
    File stats_output = "${sample_name}.stats"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    cpu: cpu
    memory: memory
    disks: "local-disk ${disk}"
  }
}

