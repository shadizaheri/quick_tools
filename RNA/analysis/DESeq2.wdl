version 1.0

workflow deseq2_workflow {
  input {
    File count_matrix
    File sample_metadata
    String output_prefix

    String deseq2_memory = "8G"
    Int deseq2_cpu = 2
    String deseq2_disk = "local-disk 50 HDD"
    Int deseq2_preemptible = 1
  }

  call RunDESeq2 {
    input:
      count_matrix = count_matrix,
      sample_metadata = sample_metadata,
      output_prefix = output_prefix,
      memory = deseq2_memory,
      cpu = deseq2_cpu,
      disk = deseq2_disk,
      preemptible = deseq2_preemptible
  }

  output {
    File results_table = RunDESeq2.results_table
    File log_file = RunDESeq2.log_file
  }
}

task RunDESeq2 {
  input {
    File count_matrix
    File sample_metadata
    String output_prefix

    # Runtime parameters
    String memory
    Int cpu
    String disk
    Int preemptible
  }

  command <<<
    echo '
    args <- commandArgs(trailingOnly = TRUE)
    count_file <- args[1]
    meta_file <- args[2]
    prefix <- args[3]

    library(DESeq2)

    counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t")
    meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t")
    counts <- counts[, rownames(meta)]

    dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
    dds <- DESeq(dds)
    res <- results(dds)
    res_df <- as.data.frame(res)
    write.table(res_df, file = paste0(prefix, "_deseq2_results.tsv"), sep = "\t", quote = FALSE)
    ' > run_deseq2.R

    Rscript run_deseq2.R ~{count_matrix} ~{sample_metadata} ~{output_prefix} > deseq2.log 2>&1
  >>>

  output {
    File results_table = "~{output_prefix}_deseq2_results.tsv"
    File log_file = "deseq2.log"
  }

  runtime {
    docker: "bioconductor/bioconductor_docker:RELEASE_3_17"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

  meta {
    author: "Shadi Zaheri"
    version: "1.0.0"
    description: "Runs DESeq2 on RNA-seq count data using inline R and configurable resources."
  }
}
