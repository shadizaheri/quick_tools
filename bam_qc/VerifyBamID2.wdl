version 1.0

## VerifyBamID2 Workflow for Contamination Detection
## 
## This workflow runs VerifyBamID2 to detect DNA sample contamination in CRAM files.
## VerifyBamID2 estimates the proportion of reads from a different individual than 
## the supposed source by comparing allele frequencies at known variant sites.
##
## Input Requirements:
## - CRAM file and its index (.crai)
## - Reference genome (FASTA, index, and dictionary)
## - SNP VCF file with common variants and its index
## - Sample identifier for output naming
##
## Outputs:
## - .selfSM file containing contamination metrics
## - Log file with detailed run information

workflow RunVerifyBamID2 {
  input {
    File cram                    # Input CRAM alignment file
    File crai                    # CRAM index file (.crai)
    File ref_fasta              # Reference genome FASTA file
    File ref_fasta_index        # Reference FASTA index (.fai)
    File ref_dict               # Reference dictionary file (.dict)
    File snp_vcf                # VCF file with SNP sites for contamination checking
    File snp_vcf_index          # VCF index file (.tbi or .csi)
    String sample_id            # Sample identifier for output file naming
    String docker_image = "us.gcr.io/broad-dsp-lrma/verifybamid2:v2.0.1"  # Docker container image
  }
  
  # Call the main task to run VerifyBamID2
  call VerifyBamID2Task {
    input:
      cram = cram,
      crai = crai,
      ref_fasta = ref_fasta,
      snp_vcf = snp_vcf,
      snp_vcf_index = snp_vcf_index,
      sample_id = sample_id,
      docker_image = docker_image
  }
  
  # Workflow outputs
  output {
    File contamination_metrics = VerifyBamID2Task.contamination_metrics
    File log_file = VerifyBamID2Task.log_file
  }

  # Workflow metadata
  meta {
    description: "Runs VerifyBamID2 to detect DNA sample contamination in CRAM files"
    version: "1.0"
    author: "Broad Institute"
    email: "support@broadinstitute.org"
    keywords: ["contamination", "quality-control", "genomics", "verifybamid2"]
    
    # External references
    verifybamid2_paper: "https://genome.cshlp.org/content/22/5/1053"
    verifybamid2_github: "https://github.com/Griffan/VerifyBamID"
    
    # Typical runtime for various input sizes
    runtime_notes: "Runtime typically 5-30 minutes depending on CRAM size and SNP density"
  }

  # Parameter descriptions for workflow inputs
  parameter_meta {
    cram: {
      description: "Input CRAM alignment file to check for contamination",
      category: "required",
      help: "Must be properly aligned to the same reference genome used in ref_fasta"
    }
    crai: {
      description: "Index file for the input CRAM file",
      category: "required", 
      help: "Should have same base name as CRAM file with .crai extension"
    }
    ref_fasta: {
      description: "Reference genome FASTA file used for alignment",
      category: "required",
      help: "Must match the reference used to create the CRAM file"
    }
    ref_fasta_index: {
      description: "FASTA index file (.fai) for the reference genome",
      category: "required"
    }
    ref_dict: {
      description: "Sequence dictionary file (.dict) for the reference genome", 
      category: "required",
      help: "Can be created with Picard CreateSequenceDictionary"
    }
    snp_vcf: {
      description: "VCF file containing SNP sites for contamination estimation",
      category: "required",
      help: "Should contain common population variants (e.g., from 1000 Genomes)"
    }
    snp_vcf_index: {
      description: "Index file for the SNP VCF (.tbi or .csi)",
      category: "required"
    }
    sample_id: {
      description: "Identifier for the sample, used in output file naming",
      category: "required",
      help: "Should be unique and descriptive for the sample"
    }
    docker_image: {
      description: "Docker container image containing VerifyBamID2",
      category: "optional",
      help: "Default uses latest Broad Institute VerifyBamID2 image"
    }
  }
}

## Task to execute VerifyBamID2 contamination detection
task VerifyBamID2Task {
  input {
    File cram                    # Input CRAM file
    File crai                    # CRAM index file  
    File ref_fasta              # Reference genome FASTA
    File snp_vcf                # SNP VCF for contamination checking
    File snp_vcf_index          # VCF index file
    String sample_id            # Sample identifier
    String docker_image         # Docker container image
  }
  
  command {
    set -e  # Exit on any error
    
    # Run VerifyBamID2 with verbose output
    # --Verbose: Enable detailed logging
    # --Bam: Input alignment file (works with CRAM too)
    # --Reference: Reference genome FASTA
    # --Site: VCF file with variant sites 
    # --SiteIndex: VCF index file
    # --Output: Output file prefix
    verifybamid2 \
      --Verbose \
      --Bam ${cram} \
      --Reference ${ref_fasta} \
      --Site ${snp_vcf} \
      --SiteIndex ${snp_vcf_index} \
      --Output ${sample_id}
  }
  
  # Task outputs
  output {
    File contamination_metrics = "${sample_id}.selfSM"    # Main contamination results
    File log_file = "${sample_id}.log"                    # Detailed run log
  }
  
  # Runtime configuration
  runtime {
    docker: docker_image        # Container with VerifyBamID2 installed
    memory: "4G"                # Memory allocation
    cpu: 2                      # CPU cores
    disks: "local-disk 50 HDD"  # Disk space allocation
    preemptible: 2              # Use preemptible instances to reduce cost
    bootDiskSizeGb: 15          # Boot disk size
  }

  # Task metadata
  meta {
    description: "Runs VerifyBamID2 to estimate DNA contamination from sequencing data"
    help: "VerifyBamID2 estimates contamination by comparing observed vs expected allele frequencies"
    contamination_metrics_description: "Tab-delimited file with contamination estimates and statistics"
    log_file_description: "Detailed log file with runtime information and potential errors"
    FREEMIX_column: "Estimated contamination fraction (0-1)"
    FREELK1_column: "Log-likelihood of estimated contamination"
    FREELK0_column: "Log-likelihood of no contamination"
    CHIPMIX_column: "Contamination estimate from array data (if available)"
    SNPS_column: "Number of SNPs used for estimation"
    READS_column: "Number of reads examined"
    AVG_DP_column: "Average depth at examined sites"
  }

  # Task parameter descriptions
  parameter_meta {
    cram: {
      description: "Input CRAM alignment file",
      stream: true,
      help: "File will be streamed to reduce storage requirements"
    }
    crai: {
      description: "CRAM index file enabling random access",
      localization_optional: false
    }
    ref_fasta: {
      description: "Reference genome in FASTA format",
      help: "Must match reference used for CRAM alignment"
    }
    snp_vcf: {
      description: "VCF file with variant sites for contamination estimation",
      help: "Common variants work best (MAF > 0.05 recommended)"
    }
    snp_vcf_index: {
      description: "Tabix index for the SNP VCF file"
    }
    sample_id: {
      description: "Sample identifier used for output file naming",
      pattern: "[A-Za-z0-9_-]+",
      help: "Should contain only alphanumeric characters, underscores, and hyphens"
    }
    docker_image: {
      description: "Docker image containing VerifyBamID2 software",
      default: "us.gcr.io/broad-dsp-lrma/verifybamid2:v2.0.1"
    }
  }
}