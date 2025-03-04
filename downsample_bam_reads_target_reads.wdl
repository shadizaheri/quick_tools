version 1.0

workflow DownsampleBam {
  parameter_meta {
    bam: "Path to the input BAM file."
    bai: "Path to the input BAI index file."
    target_reads: "Target number of reads after downsampling."
    sample_id: "Sample name used for output file naming."
    memory_gb_count_reads: "Memory allocated for the CountReads task in gigabytes."
    disk_gb_count_reads: "Disk space allocated for the CountReads task in gigabytes."
    memory_gb_downsample: "Memory allocated for the Downsample task in gigabytes."
    disk_gb_downsample: "Disk space allocated for the Downsample task in gigabytes."
    cpu_count_reads: "Number of CPU cores allocated for CountReads."
    cpu_downsample: "Number of CPU cores allocated for Downsample."
    docker: "Docker image used to run the tasks."
  }

  input {
    File bam
    File bai
    Int target_reads = 75000000
    String sample_id
    Int memory_gb_count_reads = 4
    Int disk_gb_count_reads = 10
    Int memory_gb_downsample = 8
    Int disk_gb_downsample = 50
    Int cpu_count_reads = 1
    Int cpu_downsample = 2
    String docker = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
  }
  
  call CountReads { 
    input: 
      bam = bam,
      memory_gb = memory_gb_count_reads,
      disk_gb = disk_gb_count_reads,
      cpu = cpu_count_reads,
      docker = docker
  }
  
  call Downsample {
    input:
      bam = bam,
      bai = bai,
      total_reads = CountReads.total_reads,
      target_reads = target_reads,
      sample_id = sample_id,
      memory_gb = memory_gb_downsample,
      disk_gb = disk_gb_downsample,
      cpu = cpu_downsample,
      docker = docker
  }
  
  output {
    File downsampled_bam = Downsample.downsampled_bam
    File downsampled_bai = Downsample.downsampled_bai
  }
}

# Step 1: Count total reads
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
    disks: "~{disk_gb} GiB"
    cpu: cpu
  }
}

# Step 2: Downsample BAM
task Downsample {
  input {
    File bam
    File bai
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
    samtools view -b -s $fraction ${bam} -o ${sample_id}_downsampled.bam
    samtools index ${sample_id}_downsampled.bam
  }
  output {
    File downsampled_bam = "${sample_id}_downsampled.bam"
    File downsampled_bai = "${sample_id}_downsampled.bam.bai"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "~{disk_gb} GiB"
    cpu: cpu
  }
}
