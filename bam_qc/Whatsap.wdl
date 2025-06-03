version 1.0

workflow WhatshapPhasing {
  input {
    File vcf_gz
    File vcf_gz_tbi
    File bam
    File bai
    File reference_fasta
    File reference_fasta_fai
    Boolean output_haplotagged_bam = true
    Boolean distrust_genotypes = true
    Boolean include_indels = true
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "whatshap/whatshap:latest"
  }

  call WhatshapPhase {
    input:
      vcf_gz = vcf_gz,
      vcf_gz_tbi = vcf_gz_tbi,
      bam = bam,
      bai = bai,
      reference_fasta = reference_fasta,
      reference_fasta_fai = reference_fasta_fai,
      output_haplotagged_bam = output_haplotagged_bam,
      distrust_genotypes = distrust_genotypes,
      include_indels = include_indels,
      memory_gb = memory_gb,
      cpu_cores = cpu_cores,
      docker_image = docker_image
  }

  output {
    File phased_vcf = WhatshapPhase.phased_vcf
    File phased_vcf_index = WhatshapPhase.phased_vcf_index
    File? haplotagged_bam = WhatshapPhase.haplotagged_bam
    File? haplotagged_bai = WhatshapPhase.haplotagged_bai
  }
}

task WhatshapPhase {
  input {
    File vcf_gz
    File vcf_gz_tbi
    File bam
    File bai
    File reference_fasta
    File reference_fasta_fai

    Boolean output_haplotagged_bam
    Boolean distrust_genotypes
    Boolean include_indels

    Int memory_gb
    Int cpu_cores
    String docker_image
  }

  command <<<
    set -euxo pipefail

    # Build optional flags
    EXTRA_FLAGS=""
    if [[ ~{distrust_genotypes} == true ]]; then
      EXTRA_FLAGS="$EXTRA_FLAGS --distrust-genotypes"
    fi
    if [[ ~{include_indels} == true ]]; then
      EXTRA_FLAGS="$EXTRA_FLAGS --indels"
    fi
    if [[ ~{output_haplotagged_bam} == true ]]; then
      EXTRA_FLAGS="$EXTRA_FLAGS --output-haplotagged-bam haplotagged.bam"
    fi

    # Run whatshap
    whatshap phase \
      --reference ~{reference_fasta} \
      $EXTRA_FLAGS \
      --output phased.vcf.gz \
      ~{vcf_gz} \
      ~{bam}

    tabix phased.vcf.gz
  >>>

  output {
    File phased_vcf = "phased.vcf.gz"
    File phased_vcf_index = "phased.vcf.gz.tbi"
    File? haplotagged_bam = "haplotagged.bam"
    File? haplotagged_bai = "haplotagged.bam.bai"
  }

  runtime {
    memory: "~{memory_gb} GiB"
    cpu: cpu_cores
    docker: docker_image
    disks: "local-disk 50 SSD"
  }
}

