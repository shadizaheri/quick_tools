version 1.0

workflow deseq2_workflow {
  # Runs DESeq2 differential expression analysis across samples, using tissue types
  # (or other groupings) as conditions. Outputs DE results, normalized counts,
  # and PCA plots for visualization.

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
    File normalized_counts = RunDESeq2.normalized_counts
    File pca_plot = RunDESeq2.pca_plot
    File log_file = RunDESeq2.log_file
  }
}

task RunDESeq2 {
  # Runs DESeq2 differential expression analysis.
  # Inputs:
  #   - count_matrix: gene expression raw counts table.
  #   - sample_metadata: metadata linking samples to tissue types or conditions.
  #   - output_prefix: prefix for all generated output files.
  #
  # Outputs:
  #   - Differential expression results table.
  #   - Normalized counts table.
  #   - PCA plot visualizing sample clustering.
  #   - Execution log.

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
    library(ggplot2)

    # Load counts
    counts <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t", skip = 2)
    if ("Description" %in% colnames(counts)) {
      counts <- counts[, !(colnames(counts) %in% c("Description"))]
    }

    # Load metadata
    meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t")
    counts <- counts[, rownames(meta)]

    # Create DESeq2 dataset
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
    dds <- DESeq(dds)

    # Differential Expression Results
    res <- results(dds)
    res_df <- as.data.frame(res)
    write.table(res_df, file = paste0(prefix, "_deseq2_results.tsv"), sep = "\t", quote = FALSE)

    # Normalized Counts
    norm_counts <- counts(dds, normalized=TRUE)
    write.table(as.data.frame(norm_counts), file = paste0(prefix, "_normalized_counts.tsv"), sep = "\t", quote = FALSE)
    missing_samples <- setdiff(rownames(meta), colnames(counts))
    if (length(missing_samples) > 0) {
      stop(paste("Samples in metadata but not in counts:", paste(missing_samples, collapse=", ")))
    }

    # PCA Plot
    vsd <- vst(dds)
    pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))

    p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
      geom_point(size=3) +
      xlab(paste0("PC1: ", percentVar[1], "% variance")) +
      ylab(paste0("PC2: ", percentVar[2], "% variance")) +
      theme_minimal()

    ggsave(filename = paste0(prefix, "_pca_plot.pdf"), plot = p)
    ' > run_deseq2.R

    Rscript run_deseq2.R ~{count_matrix} ~{sample_metadata} ~{output_prefix} > deseq2.log 2>&1
  >>>

  output {
    File results_table = "~{output_prefix}_deseq2_results.tsv"
    File normalized_counts = "~{output_prefix}_normalized_counts.tsv"
    File pca_plot = "~{output_prefix}_pca_plot.pdf"
    File log_file = "deseq2.log"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/deseq2:v1"
    memory: memory
    cpu: cpu
    disks: disk
    preemptible: preemptible
  }

  meta {
    author: "Shadi Zaheri"
    version: "1.1.0"
    description: "Runs DESeq2 on RNA-seq count data across conditions (e.g., tissues), outputs differential results, normalized counts, and PCA visualization."
  }
}
