version 1.0

task count_clips {
  input {
    File bam_or_cram
    File reference_fasta
    File reference_fasta_index
    Int cpu = 2
    String memory = "4G"
    String disks = "local-disk 50 HDD"
  }

  command {
    set -e

    # point samtools at your FASTA (the .fai will be next to it)
    REF_OPTS="-T ${reference_fasta}"

    # Primary reads (neither secondary 0x100 nor supplementary 0x800)
    primary_soft=$(samtools view $REF_OPTS -F 0x100 -F 0x800 ${bam_or_cram} \
                  | cut -f6 | grep -cE '[0-9]+S')
    primary_hard=$(samtools view $REF_OPTS -F 0x100 -F 0x800 ${bam_or_cram} \
                  | cut -f6 | grep -cE '[0-9]+H')

    # Secondary reads (0x100 but not 0x800)
    secondary_soft=$(samtools view $REF_OPTS -f 0x100 -F 0x800 ${bam_or_cram} \
                    | cut -f6 | grep -cE '[0-9]+S')
    secondary_hard=$(samtools view $REF_OPTS -f 0x100 -F 0x800 ${bam_or_cram} \
                    | cut -f6 | grep -cE '[0-9]+H')

    # Supplementary reads (0x800)
    supp_soft=$(samtools view $REF_OPTS -f 0x800 ${bam_or_cram} \
               | cut -f6 | grep -cE '[0-9]+S')
    supp_hard=$(samtools view $REF_OPTS -f 0x800 ${bam_or_cram} \
               | cut -f6 | grep -cE '[0-9]+H')

    # write header + one line of counts
    printf "file\tprimary_soft\tprimary_hard\tsecondary_soft\tsecondary_hard\tsupplementary_soft\tsupplementary_hard\n" \
      > counts.tsv

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${bam_or_cram}" \
      "$primary_soft" "$primary_hard" \
      "$secondary_soft" "$secondary_hard" \
      "$supp_soft"   "$supp_hard"   \
      >> counts.tsv
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
    File reference_fasta
    File reference_fasta_index
    Int cpu = 2
    String memory = "4G"
    String disks = "local-disk 100 SSD"
  }

  scatter (bam in bam_files) {
    call count_clips {
      input:
        bam_or_cram           = bam,
        reference_fasta       = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        cpu                   = cpu,
        memory                = memory,
        disks                 = disks
    }
  }

  output {
    Array[File] counts_reports = count_clips.counts
  }
}
