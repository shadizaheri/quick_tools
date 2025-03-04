version 1.0

workflow DownsampleBam {
  input {
    File bam
    File bai
    Int target_reads = 75000000
    String sample_id
    String count_reads_memory = "4G"
    String count_reads_disk = "10G"
    String downsample_memory = "8G"
    String downsample_disk = "50G"
  }
  
  call CountReads { 
    input: 
      bam = bam,
      memory = count_reads_memory,
      disk = count_reads_disk
  }
  
  call Downsample {
    input:
      bam = bam,
      bai = bai,
      total_reads = CountReads.total_reads,
      target_reads = target_reads,
      sample_id = sample_id,
      memory = downsample_memory,
      disk = downsample_disk
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
    String memory
    String disk
  }
  command {
    samtools view -c ${bam} > total_reads.txt
  }
  output {
    Int total_reads = read_int("total_reads.txt")
  }
  runtime {
    docker: "biocontainers/samtools:v1.10.0_cv1"
    memory: memory
    cpu: 1
    disk: disk
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
    String memory
    String disk
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
    docker: "biocontainers/samtools:v1.10.0_cv1"
    memory: memory
    cpu: 2
    disk: disk
  }
}
