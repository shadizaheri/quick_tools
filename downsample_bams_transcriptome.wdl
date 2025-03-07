version 1.0

workflow DownsampleSTARBAMs {
  parameter_meta {
    genome_bam: "Path to the input STAR genome-aligned BAM file."
    transcriptome_bam: "Path to the input STAR transcriptome BAM file."
    genome_bai: "Path to the input genome BAM index file."
    transcriptome_bai: "Path to the input transcriptome BAM index file."
    target_reads: "Target number of reads after downsampling. Default: 75,000,000."
    sample_id: "Sample name used for output file naming."
    memory_gb_count_reads: "Memory allocated for the CountReads task in gigabytes."
    disk_gb_count_reads: "Disk space allocated for the CountReads task in gigabytes."
    memory_gb_downsample_genome: "Memory allocated for DownsampleGenomeBAM task in gigabytes."
    disk_gb_downsample_genome: "Disk space allocated for DownsampleGenomeBAM task in gigabytes."
    memory_gb_downsample_transcriptome: "Memory allocated for DownsampleTranscriptomeBAM task in gigabytes."
    disk_gb_downsample_transcriptome: "Disk space allocated for DownsampleTranscriptomeBAM task in gigabytes."
    memory_gb_count_downsampled: "Memory allocated for the CountDownsampledReads task in gigabytes."
    disk_gb_count_downsampled: "Disk space allocated for the CountDownsampledReads task in gigabytes."
    cpu_count_reads: "Number of CPU cores allocated for CountReads."
    cpu_downsample_genome: "Number of CPU cores allocated for DownsampleGenomeBAM."
    cpu_downsample_transcriptome: "Number of CPU cores allocated for DownsampleTranscriptomeBAM."
    cpu_count_downsampled: "Number of CPU cores allocated for CountDownsampledReads."
    docker: "Docker image used to run the tasks."
  }

  input {
    File genome_bam
    File transcriptome_bam
    File genome_bai
    File transcriptome_bai
    Int target_reads = 75000000
    String sample_id
    Int memory_gb_count_reads = 4
    Int disk_gb_count_reads = 10
    Int memory_gb_downsample_genome = 16
    Int disk_gb_downsample_genome = 100
    Int memory_gb_downsample_transcriptome = 16
    Int disk_gb_downsample_transcriptome = 100
    Int memory_gb_count_downsampled = 4
    Int disk_gb_count_downsampled = 10
    Int cpu_count_reads = 1
    Int cpu_downsample_genome = 4
    Int cpu_downsample_transcriptome = 4
    Int cpu_count_downsampled = 1
    String docker = "quay.io/biocontainers/samtools:1.15--h1170115_0"
  }
  
  call CountReads { 
    input: 
      bam = genome_bam,
      memory_gb = memory_gb_count_reads,
      disk_gb = disk_gb_count_reads,
      cpu = cpu_count_reads,
      docker = docker
  }
  
  call DownsampleGenomeBAM {
    input:
      genome_bam = genome_bam,
      genome_bai = genome_bai,
      total_reads = CountReads.total_reads,
      target_reads = target_reads,
      sample_id = sample_id,
      memory_gb = memory_gb_downsample_genome,
      disk_gb = disk_gb_downsample_genome,
      cpu = cpu_downsample_genome,
      docker = docker
  }

  call DownsampleTranscriptomeBAM {
    input:
      transcriptome_bam = transcriptome_bam,
      transcriptome_bai = transcriptome_bai,
      read_ids = DownsampleGenomeBAM.read_ids,
      sample_id = sample_id,
      memory_gb = memory_gb_downsample_transcriptome,
      disk_gb = disk_gb_downsample_transcriptome,
      cpu = cpu_downsample_transcriptome,
      docker = docker
  }
  
  call CountDownsampledReads {
    input:
      genome_bam = DownsampleGenomeBAM.downsampled_bam,
      transcriptome_bam = DownsampleTranscriptomeBAM.sorted_transcriptome_bam,
      memory_gb = memory_gb_count_downsampled,
      disk_gb = disk_gb_count_downsampled,
      cpu = cpu_count_downsampled,
      docker = docker
  }
  
  output {
    File downsampled_genome_bam = DownsampleGenomeBAM.downsampled_bam
    File downsampled_transcriptome_bam = DownsampleTranscriptomeBAM.sorted_transcriptome_bam
  }
}

task CountReads {
  input {
    File bam
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    samtools view -c ${bam} > total_reads.txt
  }
  output {
    Int total_reads = read_int("total_reads.txt")
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

task DownsampleGenomeBAM {
  input {
    File genome_bam
    File genome_bai
    Int total_reads
    Int target_reads
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }

  command {
    fraction=$(echo "scale=6; ${target_reads}/${total_reads}" | bc)
    samtools view -b -s $fraction -f 1 -F 12 ${genome_bam} -o downsampled_genome.bam
    samtools index downsampled_genome.bam
    samtools view downsampled_genome.bam | cut -f1 | sort | uniq > read_ids.txt
  }

  output {
    File downsampled_bam = "downsampled_genome.bam"
    File read_ids = "read_ids.txt"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

task DownsampleTranscriptomeBAM {
  input {
    File transcriptome_bam
    File transcriptome_bai
    File read_ids
    String sample_id
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }

  command {
    samtools view -H ${transcriptome_bam} > header.sam
    samtools view ${transcriptome_bam} | grep -F -f read_ids.txt > filtered.sam
    cat header.sam filtered.sam | samtools view -bS - > filtered_transcriptome.bam
    samtools sort -n -o sorted_transcriptome.bam filtered_transcriptome.bam
    samtools index sorted_transcriptome.bam
  }

  output {
    File sorted_transcriptome_bam = "sorted_transcriptome.bam"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

task CountDownsampledReads {
  input {
    File genome_bam
    File transcriptome_bam
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }

  command {
    samtools view -c ${genome_bam} > downsampled_total_reads.txt
    samtools view -c ${transcriptome_bam} > downsampled_transcriptome_total_reads.txt
  }

  output {
    Int downsampled_total_reads = read_int("downsampled_total_reads.txt")
    Int downsampled_transcriptome_total_reads = read_int("downsampled_transcriptome_total_reads.txt")
  }

  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}
