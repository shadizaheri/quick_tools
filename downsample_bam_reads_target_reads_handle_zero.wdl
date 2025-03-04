version 1.0

workflow DownsampleBam {
  parameter_meta {
    bam: "Path to the input BAM file."
    bai: "Path to the input BAI index file."
    target_reads: "Target number of reads after downsampling."
    sample_id: "Sample name used for output file naming."
    preemptible_tries: "Number of preemptible retries allowed for each task."
    memory_gb: "Memory allocated for each task in gigabytes."
    disk_gb: "Disk space allocated for each task in gigabytes."
    cpu: "Number of CPU cores allocated for each task."
    docker: "Docker image used to run the tasks."
  }

  input {
    File bam
    File bai
    Int target_reads = 75000000
    String sample_id
    Int preemptible_tries = 1
    Int memory_gb = 8
    Int disk_gb = 50
    Int cpu = 2
    String docker = "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
  }
  
  call CountReads { 
    input: 
      bam = bam,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
      docker = docker
  }
  
  call Downsample {
    input:
      bam = bam,
      bai = bai,
      total_reads = CountReads.total_reads,
      target_reads = target_reads,
      sample_id = sample_id,
      preemptible_tries = preemptible_tries,
      memory_gb = memory_gb,
      disk_gb = disk_gb,
      cpu = cpu,
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
    Int preemptible_tries
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
    preemptible: preemptible_tries
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
    Int preemptible_tries
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    fraction=$(echo "scale=6; if (${total_reads} < ${target_reads}) print 1 else print ${target_reads}/${total_reads}" | bc)
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
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: preemptible_tries
    cpu: cpu
  }
}
