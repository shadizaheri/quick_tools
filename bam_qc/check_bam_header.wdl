version 1.0

workflow ExtractReadGroupHeader {
  input {
    String sample_name
    File cram_file
    String disk_space  # e.g., "10 HDD"
  }

  call ExtractHeader {
    input:
      sample_name = sample_name,
      cram_file = cram_file,
      disk_space = disk_space
  }

  output {
    File header_txt = ExtractHeader.header_txt
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
    gsutil cat ~{cram_file} | samtools view -H - | grep -E '^@RG|^@PG' > "~{sample_name}_header.txt"
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

