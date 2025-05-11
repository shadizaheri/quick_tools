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
    # Ensure input exists
    FILE=\"${bam_file}\"
    if [[ ! -f \"$FILE\" ]]; then
      echo "Error: File not found: $FILE" >&2
      exit 2
    fi

    # Function to count clipped reads
    count_clips() {
      local label=$1
      local flags=$2

    # count soft‑clips in CIGAR
    local soft=$(samtools view $flags "${FILE}" \
      | cut -f6 \
      | grep -cE '[0-9]+S')

    # count hard‑clips in CIGAR
  local hard=$(samtools view $flags "${FILE}" \
    | cut -f6 \
    | grep -cE '[0-9]+H')

  echo "${label} — Soft-clips: $soft; Hard-clips: $hard"
}

# Run counts and write to TSV
{
  count_clips "Primary"     "-F 0x100 -F 0x800"
  count_clips "Secondary"   "-f 0x100 -F 0x800"
  count_clips "Supplementary" "-f 0x800"
} > counts.tsv

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