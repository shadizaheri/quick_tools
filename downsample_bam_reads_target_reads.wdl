version 1.0

workflow DownsampleBam {
  parameter_meta {
    bam: "Path to the input genomic BAM file."
    transcriptome_bam: "Path to the input transcriptome BAM file."
    bai: "Path to the input BAI index file."
    target_reads: "Target number of reads after downsampling."
    sample_id: "Sample name used for output file naming."
    memory_gb_count_reads: "Memory allocated for the CountReads task in gigabytes."
    disk_gb_count_reads: "Disk space allocated for the CountReads task in gigabytes."
    memory_gb_downsample: "Memory allocated for the Downsample task in gigabytes."
    disk_gb_downsample: "Disk space allocated for the Downsample task in gigabytes."
    memory_gb_count_downsampled: "Memory allocated for the CountDownsampledReads task in gigabytes."
    disk_gb_count_downsampled: "Disk space allocated for the CountDownsampledReads task in gigabytes."
    cpu_count_reads: "Number of CPU cores allocated for CountReads."
    cpu_downsample: "Number of CPU cores allocated for Downsample."
    cpu_count_downsampled: "Number of CPU cores allocated for CountDownsampledReads."
    docker: "Docker image used to run the tasks."
  }

  input {
    File bam
    File transcriptome_bam
    File bai
    Int target_reads = 75000000
    String sample_id
    Int memory_gb_count_reads = 4
    Int disk_gb_count_reads = 10
    Int memory_gb_downsample = 8
    Int disk_gb_downsample = 50
    Int memory_gb_count_downsampled = 4
    Int disk_gb_count_downsampled = 10
    Int cpu_count_reads = 1
    Int cpu_downsample = 2
    Int cpu_count_downsampled = 1
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
      transcriptome_bam = transcriptome_bam,
      bai = bai,
      total_reads = CountReads.total_reads,
      target_reads = target_reads,
      sample_id = sample_id,
      memory_gb = memory_gb_downsample,
      disk_gb = disk_gb_downsample,
      cpu = cpu_downsample,
      docker = docker
  }
  
  call CountDownsampledReads {
    input:
      bam = Downsample.downsampled_bam,
      transcriptome_bam = Downsample.downsampled_transcriptome_bam,
      memory_gb = memory_gb_count_downsampled,
      disk_gb = disk_gb_count_downsampled,
      cpu = cpu_count_downsampled,
      docker = docker
  }
  
  output {
    File downsampled_bam = Downsample.downsampled_bam
    File downsampled_bai = Downsample.downsampled_bai
    File downsampled_transcriptome_bam = Downsample.downsampled_transcriptome_bam
    File downsampled_transcriptome_bai = Downsample.downsampled_transcriptome_bai
    Int downsampled_total_reads = CountDownsampledReads.downsampled_total_reads
    Int downsampled_transcriptome_total_reads = CountDownsampledReads.downsampled_transcriptome_total_reads
  }
}

# Count total reads
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

# Downsample Both Genomic and Transcriptome BAMs with Paired-End Read Retention
task Downsample {
  input {
    File bam
    File transcriptome_bam
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
    
    # Downsample genomic BAM while retaining paired-end reads
    samtools view -b -s $fraction -f 1 -F 12 ${bam} -o tmp_genomic.bam  # Retain properly paired reads

    # Extract retained read names
    samtools view tmp_genomic.bam | cut -f1 | sort | uniq > retained_reads.txt

    # Downsample transcriptome BAM based on retained reads
    samtools view -b -N retained_reads.txt ${transcriptome_bam} -o tmp_transcriptome.bam

    # Ensure BAMs are not empty
    if [ ! -s tmp_genomic.bam ]; then
        echo "Error: The downsampled genomic BAM is empty." >&2
        exit 1
    fi
    if [ ! -s tmp_transcriptome.bam ]; then
        echo "Error: The downsampled transcriptome BAM is empty." >&2
        exit 1
    fi

    # Sort genomic BAM by coordinate before indexing
    samtools sort -o ${sample_id}_downsampled.bam tmp_genomic.bam
    samtools index ${sample_id}_downsampled.bam

    # Sort transcriptome BAM by name (important for RSEM compatibility)
    samtools sort -n -o ${sample_id}_downsampled_Aligned.toTranscriptome.out.bam tmp_transcriptome.bam
    samtools index ${sample_id}_downsampled_Aligned.toTranscriptome.out.bam
  }
  output {
    File downsampled_bam = "~{sample_id}_downsampled.bam"
    File downsampled_bai = "~{sample_id}_downsampled.bam.bai"
    File downsampled_transcriptome_bam = "~{sample_id}_downsampled_Aligned.toTranscriptome.out.bam"
    File downsampled_transcriptome_bai = "~{sample_id}_downsampled_Aligned.toTranscriptome.out.bam.bai"
  }
  runtime {
    docker: docker
    memory: "~{memory_gb} GiB"
    disks: "local-disk ~{disk_gb} HDD"
    cpu: cpu
  }
}

# Count reads after downsampling
task CountDownsampledReads {
  input {
    File bam
    File transcriptome_bam
    Int memory_gb
    Int disk_gb
    Int cpu
    String docker
  }
  command {
    samtools view -c ${bam} > downsampled_total_reads.txt
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
