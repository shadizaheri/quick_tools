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
    String output_bam_name = "haplotagged.bam"
    Int memory_gb = 8
    Int cpu_cores = 2
    String docker_image = "us.gcr.io/broad-dsp-lrma/lr-whatshap:2.2"
  }

  call whatshap_phase {
    input:
      vcf = vcf_gz,
      vcf_index = vcf_gz_tbi,
      aligned_bams = [bam],
      aligned_bam_indices = [bai],
      reference = reference_fasta,
      reference_index = reference_fasta_fai
  }

  call whatshap_haplotag {
    input:
      phased_vcf = whatshap_phase.phased_vcf,
      phased_vcf_index = whatshap_phase.phased_vcf_index,
      aligned_bam = bam,
      aligned_bam_index = bai,
      reference = reference_fasta,
      reference_index = reference_fasta_fai,
      output_bam_name = output_bam_name,
  }

  output {
    File phased_vcf = whatshap_phase.phased_vcf
    File phased_vcf_index = whatshap_phase.phased_vcf_index
    File? haplotagged_bam = whatshap_haplotag.haplotagged_bam
    File? haplotagged_bai = whatshap_haplotag.haplotagged_bam_index
  }
}
task whatshap_phase {
	input {
		File vcf
		File vcf_index
		String? chromosome

		Array[File] aligned_bams
		Array[File] aligned_bam_indices

		File reference
		File reference_index
	}

	String vcf_basename = basename(vcf, ".vcf.gz")
	Int disk_size = ceil((size(vcf, "GB") + size(reference, "GB") + size(aligned_bams, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap --version

		whatshap phase \
			--indels \
			--reference ~{reference} \
			~{"--chromosome " + chromosome} \
			--output ~{vcf_basename}.phased.vcf.gz \
			~{vcf} \
			~{sep=' ' aligned_bams}

		tabix ~{vcf_basename}.phased.vcf.gz
	>>>

	output {
		File phased_vcf = "~{vcf_basename}.phased.vcf.gz"
		File phased_vcf_index = "~{vcf_basename}.phased.vcf.gz.tbi"
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-lrma/lr-whatshap:2.2"
		cpu: 2
		memory: "8 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 0
		maxRetries: 0
	}
}

task whatshap_haplotag {
	input {
		File phased_vcf
		File phased_vcf_index

		File aligned_bam
		File aligned_bam_index

		File reference
		File reference_index

		String? output_bam_name
	}

	String output_bam = select_first([output_bam_name, "~{basename(aligned_bam, '.bam')}.haplotagged.bam"])
	Int threads = 4
	Int disk_size = ceil((size(phased_vcf, "GB") + size(aligned_bam, "GB") + size(reference, "GB")) * 2 + 20)

	command <<<
		set -euo pipefail

		whatshap --version

		whatshap haplotag \
			--tag-supplementary \
			--output-threads ~{threads} \
			--reference ~{reference} \
			--output ~{output_bam} \
			~{phased_vcf} \
			~{aligned_bam}

		samtools --version

		samtools index \
			-@ ~{threads - 1} \
			~{output_bam}
	>>>

	output {
		File haplotagged_bam = "~{output_bam}"
		File haplotagged_bam_index = "~{output_bam}.bai"
	}

	runtime {
		docker: "us.gcr.io/broad-dsp-lrma/lr-whatshap:2.2"
		cpu: threads
		memory: "4 GB"
		disk: disk_size + " GB"
		disks: "local-disk " + disk_size + " HDD"
		preemptible: 0
		maxRetries: 0
	}
}

