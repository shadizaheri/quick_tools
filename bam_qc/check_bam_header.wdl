version 1.0

workflow ExtractReadGroupHeaderAndQNames {
  input {
    String sample_name
    File cram_file
    String disk_space  # e.g., "10 HDD"
    Int num_qnames_to_sample = 10000
  }

  call ExtractHeader {
    input:
      sample_name = sample_name,
      cram_file = cram_file,
      disk_space = disk_space
  }

  call SampleReadNames {
    input:
      sample_name = sample_name,
      cram_file = cram_file,
      disk_space = disk_space,
      num_reads = num_qnames_to_sample
  }

  output {
    File header_txt = ExtractHeader.header_txt
    File qnames_txt = SampleReadNames.qnames_txt
  }
}

task ExtractHeader {
  input {
    String sample_name
    File cram_file
    String disk_space
  }

  command {
    set -euo pipefail
    echo "Extracting header for ~{sample_name}"
    samtools view -H ~{cram_file} | grep -E '^@RG|^@PG' > "~{sample_name}_header.txt"
  }

  output {
    File header_txt = "~{sample_name}_header.txt"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    cpu: 1
    memory: "1G"
    disks: "local-disk ~{disk_space}"
  }
}

task SampleReadNames {
  input {
    String sample_name
    File cram_file
    String disk_space
    Int num_reads
  }

  command {
    set -euo pipefail
    echo "Sampling ~{num_reads} read names from ~{sample_name}"
    samtools view ~{cram_file} | head -n ~{num_reads} | cut -f1 > "~{sample_name}_qnames.txt"
  }

  output {
    File qnames_txt = "~{sample_name}_qnames.txt"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    cpu: 1
    memory: "1G"
    disks: "local-disk ~{disk_space}"
  }
}
