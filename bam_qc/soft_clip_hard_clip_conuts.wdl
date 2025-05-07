version 1.0

task count_clips {
  input {
    File bam             # CRAM or BAM file
    File bam_index       # corresponding index (.crai or .bai)
    Int memory_gb        # Runtime memory in GB
    Int disk_gb          # Runtime disk size in GB
  }

  command <<<
#!/usr/bin/env bash
set -euo pipefail

# Count soft‑clipped reads
soft=$(samtools view "${bam}" | grep -cE '[0-9]+S')

# Count hard‑clipped reads
hard=$(samtools view "${bam}" | grep -cE '[0-9]+H')

# Count secondary or supplementary alignments (flag 0x100 or 0x800)
sec_sup=$(samtools view "${bam}" | awk '{ if ((and($2,0x100)) || (and($2,0x800))) c++ } END { print c+0 }')

# Write combined TSV
cat <<EOF > counts.tsv
"type\tcount"
"soft_clipped\t${soft}"
"hard_clipped\t${hard}"
"secondary_supplementary\t${sec_sup}"
EOF

# Write individual count files
printf "%d" ${soft}    > soft_clipped.txt
printf "%d" ${hard}    > hard_clipped.txt
printf "%d" ${sec_sup} > secondary_supplementary.txt

echo "Counts written to counts.tsv"
>>>

  output {
    File counts_tsv             = "counts.tsv"
    Int soft_clipped            = read_int("soft_clipped.txt")
    Int hard_clipped            = read_int("hard_clipped.txt")
    Int secondary_supplementary = read_int("secondary_supplementary.txt")
  }

  runtime {
    docker: "quay.io/biocontainers/samtools:1.17--h2e538c0_1"
    cpu: 1
    memory: "${memory_gb} GB"
    disks: "local-disk ${disk_gb} SSD"
  }
}

workflow clip_counter {
  input {
    File bam_file
    File bam_index
    Int memory_gb = 4
    Int disk_gb   = 10
  }

  call count_clips {
    input:
      bam       = bam_file,
      bam_index = bam_index,
      memory_gb = memory_gb,
      disk_gb   = disk_gb
  }

  output {
    File counts_tsv             = count_clips.counts_tsv
    Int soft_clipped            = count_clips.soft_clipped
    Int hard_clipped            = count_clips.hard_clipped
    Int secondary_supplementary = count_clips.secondary_supplementary
  }
}
