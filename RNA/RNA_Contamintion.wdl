version 1.0
# Author: Shadi Zaheri
# Description: Workflow for RNA-seq contamination analysis


workflow RNAContaminationWorkflow {
  input {
    File   fastq_file
    File   reference_fasta
    File   sites_vcf
    File   ref_flat
    String sample_id
    File   bam_file
    File   bamixchecker_sites_bed

    File   kraken2_db              # <— ADDED!
    Int    fastqc_memory_gb       = 2
    Int    fastqc_disk_gb         = 10
    Int    kraken2_memory_gb      = 8
    Int    kraken2_disk_gb        = 50
    Int    rnaseqc_memory_gb      = 4
    Int    rnaseqc_disk_gb        = 20
    Int    verifybamid2_memory_gb = 8
    Int    verifybamid2_disk_gb   = 30
    Int    bamixchecker_memory_gb = 4
    Int    bamixchecker_disk_gb   = 20
  }

  call FastQC {
    input:
      fastq_file = fastq_file,
      memory_gb  = fastqc_memory_gb,
      disk_gb    = fastqc_disk_gb
  }

  call Kraken2 {
    input:
      fastq_file = fastq_file,
      kraken2_db = kraken2_db,
      memory_gb  = kraken2_memory_gb,
      disk_gb    = kraken2_disk_gb
  }

  
  call RNASeQC {
    input:
      bam_file = bam_file,
      reference_fasta = reference_fasta,
      ref_flat = ref_flat,
      memory_gb = rnaseqc_memory_gb,
      disk_gb = rnaseqc_disk_gb
  }

  call VerifyBamID2 {
    input:
      bam_file = bam_file,
      sites_vcf = sites_vcf,
      reference_fasta = reference_fasta,
      memory_gb = verifybamid2_memory_gb,
      disk_gb = verifybamid2_disk_gb
  }

  call BamixChecker {
    input:
      bam_file = bam_file,
      sites_bed = bamixchecker_sites_bed,
      memory_gb = bamixchecker_memory_gb,
      disk_gb = bamixchecker_disk_gb
  }

  output {
    File fastqc_report = FastQC.fastqc_report
    File kraken2_report = Kraken2.kraken2_report
    File rnaseqc_metrics = RNASeQC.metrics_file
    File verifybamid2_output = VerifyBamID2.output_file
    File bamixchecker_output = BamixChecker.contamination_file
  }
}

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
    File fastqc_report = glob("*_fastqc.zip")[0]
  }

  runtime {
    docker: "biocontainers/fastqc:v0.11.9_cv7"
    cpu: 1
    memory: "${memory_gb}G"
    disks: "local-disk ${disk_gb} HDD"
  }
}

task Kraken2 {
  input {
    File   fastq_file
    File   kraken2_db
    Int    memory_gb
    Int    disk_gb
  }

  command <<<                                
    kraken2 \
      --db ${kraken2_db} \
      --report kraken2_report.txt \
      --output kraken2_output.txt \
      ${fastq_file}
  >>>

  output {
    File kraken2_report = "kraken2_report.txt"
  }

  runtime {
    docker: "biocontainers/kraken2:v2.1.2"
    cpu:    4
    memory: "${memory_gb}G"
    disks:  "local-disk ${disk_gb} HDD"
  }
}


task RNASeQC {
  input {
    File bam_file
    File reference_fasta
    File ref_flat
    Int memory_gb = 4
    Int disk_gb = 20
  }

  command {
    picard CollectRnaSeqMetrics \
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

task VerifyBamID2 {
  input {
    File bam_file
    File sites_vcf
    File reference_fasta
    Int memory_gb = 8
    Int disk_gb = 30
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

task BamixChecker {
  input {
    File bam_file
    File sites_bed
    Int memory_gb
    Int disk_gb
  }

  command {
    bamixchecker --bam ${bam_file} --sites ${sites_bed} --out contamination_bamixchecker.txt
  }

  output {
    File contamination_file = "contamination_bamixchecker.txt"
  }

  runtime {
    docker: "maxplanck/bamixchecker:v1.0.1"
    cpu: 2
    memory: "${memory_gb}G"
    disks: "local-disk ${disk_gb} HDD"
  }
}
