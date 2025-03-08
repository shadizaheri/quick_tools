version 1.0

workflow FixTranscriptomeBAM {
  input {
    File genomic_bam
    File transcriptome_bam
    String sample_id
    Int memory_gb = 8
    Int disk_gb = 50
    Int cpu = 2
    String docker = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
  }

  call ExtractReadNames {
    input:
      genomic_bam = genomic_bam,
      transcriptome_bam = transcriptome_bam,
      sample_id = sample_id,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }

  call FilterTranscriptomeBAM {
    input:
      transcriptome_bam = transcriptome_bam,
      retained_reads = ExtractReadNames.retained_reads,
      sample_id = sample_id,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }

  call ValidateAndSortBAM {
    input:
      filtered_transcriptome_bam = FilterTranscriptomeBAM.filtered_transcriptome_bam,
      sample_id = sample_id,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }

  output {
    File fixed_transcriptome_bam = ValidateAndSortBAM.sorted_filtered_transcriptome_bam
  }
}

### Task 1: Extract Read Names from Both BAMs
task ExtractReadNames {
  input {
    File genomic_bam
    File transcriptome_bam
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    # Extract read names from both BAMs
    samtools view ${genomic_bam} | cut -f1 | sort | uniq > genomic_names.txt
    samtools view ${transcriptome_bam} | cut -f1 | sort | uniq > transcriptome_names.txt

    # Identify retained reads (reads that exist in genomic BAM)
    comm -12 genomic_names.txt transcriptome_names.txt > retained_reads.txt
  }
  output {
    File retained_reads = "retained_reads.txt"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

### Task 2: Filter Transcriptome BAM Based on Retained Reads
task FilterTranscriptomeBAM {
  input {
    File transcriptome_bam
    File retained_reads
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    # Keep only reads present in retained_reads.txt
    samtools view -b -N retained_reads ${transcriptome_bam} -o filtered_transcriptome.bam
  }
  output {
    File filtered_transcriptome_bam = "filtered_transcriptome.bam"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

### Task 3: Validate BAM Order and Sort if Needed
task ValidateAndSortBAM {
  input {
    File filtered_transcriptome_bam
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    # Check if the BAM is already sorted by name
    samtools sort -n -c ${filtered_transcriptome_bam} || \
    samtools sort -n -o sorted_filtered_transcriptome.bam ${filtered_transcriptome_bam}
  }
  output {
    File sorted_filtered_transcriptome_bam = "sorted_filtered_transcriptome.bam"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

