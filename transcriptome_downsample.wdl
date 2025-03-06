version 1.0

workflow DownsampleTranscriptomeBam {
  parameter_meta {
    transcriptome_bam: "Path to the input transcriptome BAM file."
    downsampled_genomic_bam: "Path to the already downsampled genomic BAM file."
    sample_id: "Sample name used for output file naming."
    memory_gb: "Memory allocated for the task in gigabytes."
    disk_gb: "Disk space allocated for the task in gigabytes."
    cpu: "Number of CPU cores allocated for the task."
    docker: "Docker image with samtools."
  }

  input {
    File transcriptome_bam
    File downsampled_genomic_bam
    String sample_id
    Int memory_gb = 4
    Int disk_gb = 50
    Int cpu = 2
    String docker = "biocontainers/samtools:v1.10.0_cv1"
  }
  
  call ExtractReadNames {
    input:
      downsampled_genomic_bam = downsampled_genomic_bam,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }
  
  call DownsampleTranscriptome {
    input:
      transcriptome_bam = transcriptome_bam,
      read_names_file = ExtractReadNames.read_names_file,
      sample_id = sample_id,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }
  
  output {
    File downsampled_transcriptome_bam = DownsampleTranscriptome.downsampled_transcriptome_bam
    File downsampled_transcriptome_bai = DownsampleTranscriptome.downsampled_transcriptome_bai
  }
}

# Step 1: Extract read names from the downsampled genomic BAM
task ExtractReadNames {
  input {
    File downsampled_genomic_bam
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    samtools view ~{downsampled_genomic_bam} | cut -f1 | sort | uniq > retained_read_names.txt
  }
  output {
    File read_names_file = "retained_read_names.txt"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} SSD"
    cpu: cpu
  }
}

# Step 2: Downsample transcriptome BAM using extracted read names
task DownsampleTranscriptome {
  input {
    File transcriptome_bam
    File read_names_file
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    samtools view -b -N ~{read_names_file} ~{transcriptome_bam} -o tmp.bam
    
    # Ensure BAM is not empty
    if [ ! -s tmp.bam ]; then
        echo "Error: The downsampled transcriptome BAM is empty." >&2
        exit 1
    fi
    
    # Sort the BAM before indexing
    samtools sort -o ~{sample_id}_downsampled_Aligned.toTranscriptome.out.bam tmp.bam
    samtools index ~{sample_id}_downsampled_Aligned.toTranscriptome.out.bam
  }
  output {
    File downsampled_transcriptome_bam = "~{sample_id}_downsampled_Aligned.toTranscriptome.out.bam"
    File downsampled_transcriptome_bai = "~{sample_id}_downsampled_Aligned.toTranscriptome.out.bam.bai"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} SSD"
    cpu: cpu
  }
}