version 1.0

workflow RNAContaminationWorkflow {

  # Input Parameters
  input {
    File fastq_file     # RNA-Seq reads (FASTQ file)
    File reference_fasta # Reference genome (FASTA)
    File sites_vcf      # Known SNP sites (VCF)
    File ref_flat       # Reference flat file for RNA-SeQC
    File kraken2_db     # Kraken2 database
    File bam_file       # Aligned RNA BAM file

    # User-adjustable resource parameters with default values
    Int fastqc_memory_gb = 2  # Default: 2 GB for FastQC
    Int fastqc_disk_gb = 10    # Default: 10 GB for FastQC
    Int kraken2_memory_gb = 8  # Default: 8 GB for Kraken2
    Int kraken2_disk_gb = 50    # Default: 50 GB for Kraken2
    Int rnaseqc_memory_gb = 4  # Default: 4 GB for RNA-SeQC
    Int rnaseqc_disk_gb = 20    # Default: 20 GB for RNA-SeQC
    Int verifybamid2_memory_gb = 8  # Default: 8 GB for VerifyBamID2
    Int verifybamid2_disk_gb = 30    # Default: 30 GB for VerifyBamID2
  }

  # Run FastQC
  call FastQC {
    input:
      fastq_file = fastq_file,
      memory_gb = fastqc_memory_gb,
      disk_gb = fastqc_disk_gb
  }

  # Run Kraken2
  call Kraken2 {
    input:
      fastq_file = fastq_file,
      kraken2_db = kraken2_db,
      memory_gb = kraken2_memory_gb,
      disk_gb = kraken2_disk_gb
  }

  # Run Picard RNA-SeQC
  call RNASeQC {
    input:
      bam_file = bam_file,
      reference_fasta = reference_fasta,
      ref_flat = ref_flat,
      memory_gb = rnaseqc_memory_gb,  # Use Elvis operator to provide default
      disk_gb = rnaseqc_disk_gb      # Use Elvis operator to provide default
  }

  # Run VerifyBamID2
  call VerifyBamID2 {
    input:
      bam_file = bam_file,
      sites_vcf = sites_vcf,
      reference_fasta = reference_fasta,
      memory_gb = verifybamid2_memory_gb,  # Use Elvis operator to provide default
      disk_gb = verifybamid2_disk_gb      # Use Elvis operator to provide default
  }

  # Outputs
  output {
    File fastqc_report = FastQC.fastqc_report
    File kraken2_report = Kraken2.kraken2_report
    File rnaseqc_metrics = RNASeQC.metrics_file
    File verifybamid2_output = VerifyBamID2.output_file
  }
}

### FastQC Task ###
task FastQC {
  input {
    File fastq_file
    Int memory_gb
    Int disk_gb
  }

  command {
    fastqc ${fastq_file} -o ./
  }

  output {
    File fastqc_report = "${basename(fastq_file)}_fastqc.zip"
  }

  runtime {
    docker: "biocontainers/fastqc:v0.11.9_cv7"
    cpu: 1
    memory: "${memory_gb}G"
    disks: "local-disk ${disk_gb} HDD"
  }
}

### Kraken2 Task ###
task Kraken2 {
  input {
    File fastq_file
    File kraken2_db
    Int memory_gb
    Int disk_gb
  }

  command {
    kraken2 --db ${kraken2_db} --report kraken2_report.txt --output kraken2_output.txt ${fastq_file}
  }

  output {
    File kraken2_report = "kraken2_report.txt"
  }

  runtime {
    docker: "biocontainers/kraken2:v2.1.2"
    cpu: 4
    memory: "${memory_gb}G"
    disks: "local-disk ${disk_gb} HDD"
  }
}

### RNA-SeQC Task ###
task RNASeQC {
  input {
    File bam_file
    File reference_fasta
    File ref_flat
    Int memory_gb = 4  # Default value if not provided
    Int disk_gb = 20    # Default value if not provided
  }

  command {
    java -jar /path/to/picard.jar CollectRnaSeqMetrics \
      I=${bam_file} \
      O=rnaseqc_metrics.txt \
      REF_FLAT=${ref_flat} \
      STRAND=NONE \
      R=${reference_fasta}
  }

  output {
    File metrics_file = "rnaseqc_metrics.txt"
  }

  runtime {
    docker: "broadinstitute/picard:2.27.4"
    cpu: 2
    memory: "${memory_gb}G"
    disks: "local-disk ${disk_gb} HDD"
  }
}

### VerifyBamID2 Task ###
task VerifyBamID2 {
  input {
    File bam_file
    File sites_vcf
    File reference_fasta
    Int memory_gb = 8  # Default value if not provided
    Int disk_gb = 30    # Default value if not provided
  }

  command {
    VerifyBamID2 --BAM ${bam_file} --VCF ${sites_vcf} --OUT contamination_results --REF ${reference_fasta}
  }

  output {
    File output_file = "contamination_results.selfSM"
  }

  runtime {
    docker: "broadinstitute/verifybamid2:v1.0.6"
    cpu: 4
    memory: "${memory_gb}G"
    disks: "local-disk ${disk_gb} HDD"
  }
}

