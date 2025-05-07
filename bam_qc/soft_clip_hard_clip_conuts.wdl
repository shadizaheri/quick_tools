version 1.0

task count_clips {
  input {
    File bam           # BAM file
    File bam_index     # .bai
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb   = 10
  }

  command <<<
#!/usr/bin/env bash
set -euo pipefail

# Bind WDL path
bam_path="${bam}"

# Build samtools command (no -T)
samtools_cmd=( samtools view -@ ${cpu} )

# Single‑pass count via awk
printf 'type\tcount\n' > counts.tsv
"${samtools_cmd[@]}" "${bam_path}" | \
  awk '
    BEGIN { soft=hard=sec=0 }
    {
      if ($6 ~ /[0-9]+S/) soft++
      if ($6 ~ /[0-9]+H/) hard++
      if (($2 & 0x100) || ($2 & 0x800)) sec++
    }
    END {
      print "soft_clipped\t" soft
      print "hard_clipped\t" hard
      print "secondary_supplementary\t" sec
    }
  ' >> counts.tsv

# Extract individual counts
awk 'NR==2{print $2 > "soft_clipped.txt"}
     NR==3{print $2 > "hard_clipped.txt"}
     NR==4{print $2 > "secondary_supplementary.txt"}' counts.tsv

echo "Counts written to counts.tsv"
>>>

  output {
    File   counts_tsv             = "counts.tsv"
    Int    soft_clipped            = read_int("soft_clipped.txt")
    Int    hard_clipped            = read_int("hard_clipped.txt")
    Int    secondary_supplementary = read_int("secondary_supplementary.txt")
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/gtex_v8_star_2.7.10a_custom"
    cpu:    cpu
    memory: "${memory_gb} GB"
    disks:  "local-disk ${disk_gb} SSD"
  }
}

workflow clip_counter {
  input {
    File bam_file
    File bam_index
    Int cpu        = 1
    Int memory_gb  = 4
    Int disk_gb    = 10
  }

  call count_clips {
    input:
      bam        = bam_file,
      bam_index  = bam_index,
      cpu        = cpu,
      memory_gb  = memory_gb,
      disk_gb    = disk_gb
  }

  output {
    File counts_tsv             = count_clips.counts_tsv
    Int  soft_clipped           = count_clips.soft_clipped
    Int  hard_clipped           = count_clips.hard_clipped
    Int  secondary_supplementary= count_clips.secondary_supplementary
  }
}
