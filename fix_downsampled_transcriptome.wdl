version 1.0

workflow FixTranscriptomeBAM {
  input {
    File transcriptome_bam
    String sample_id
    Int memory_gb = 8
    Int disk_gb = 50
    Int cpu = 2
    String docker = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
  }

  call RemoveSecondaryAlignments {
    input:
      transcriptome_bam = transcriptome_bam,
      sample_id = sample_id,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }

  call CountReads {
    input:
      bam = RemoveSecondaryAlignments.filtered_transcriptome_bam,
      sample_id = sample_id,
      docker = docker
  }

  call ValidateAndSortBAM {
    input:
      filtered_transcriptome_bam = RemoveSecondaryAlignments.filtered_transcriptome_bam,
      sample_id = sample_id,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }

  call CountReads as CountFinalReads {
    input:
      bam = ValidateAndSortBAM.sorted_filtered_transcriptome_bam,
      sample_id = sample_id,
      docker = docker
  }

  output {
    File fixed_transcriptome_bam = ValidateAndSortBAM.sorted_filtered_transcriptome_bam
    Int initial_transcriptome_reads = CountReads.total_reads
    Int final_transcriptome_reads = CountFinalReads.total_reads
  }
}

### Task 1: Remove Secondary Alignments from Transcriptome BAM
task RemoveSecondaryAlignments {
  input {
    File transcriptome_bam
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    samtools view -b -F 256 ${transcriptome_bam} -o filtered_transcriptome.bam
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

### Task 2: Count Reads in BAM File
task CountReads {
  input {
    File bam
    String sample_id
    String docker
  }
  command {
    samtools view -c ${bam} > read_count.txt
  }
  output {
    Int total_reads = read_int("read_count.txt")
  }
  runtime {
    docker: docker
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    cpu: 1
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
    # Always sort BAM by name for RSEM compatibility
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
