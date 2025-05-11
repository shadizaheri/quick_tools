version 1.0

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

task count_clips {
  input {
    File bam_file
    Int cpu
    String memory
    String disks
  }

  command <<<
    set -e

    FILE=~{bam_file}
    if [[ ! -f "$FILE" ]]; then
      echo "Error: File not found: $FILE" >&2
      exit 2
    fi

    count_clips() {
      local label=$1
      local flags=$2
      local soft=$(samtools view $flags "$FILE" \
        | cut -f6 | grep -cE '[0-9]+S')
      local hard=$(samtools view $flags "$FILE" \
        | cut -f6 | grep -cE '[0-9]+H')
      echo -e "$label\t$soft\t$hard"
    }

    # Header line
    echo -e "category\tsoft_clipped_reads\thard_clipped_reads" > counts.tsv

    # Primary (neither secondary nor supplementary)
    count_clips "Primary"   "-F 0x100 -F 0x800" >> counts.tsv

    # Secondary (0x100 but not 0x800)
    count_clips "Secondary" "-f 0x100 -F 0x800" >> counts.tsv

    # Supplementary (0x800)
    count_clips "Supplementary" "-f 0x800" >> counts.tsv
  >>>

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