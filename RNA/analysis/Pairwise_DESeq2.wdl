version 1.0

workflow deseq2_workflow {
  input {
    File count_matrix
    File sample_metadata
    String output_prefix

    String reference_tissue  # new input
    String target_tissue     # new input

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
      reference_tissue = reference_tissue,
      target_tissue = target_tissue,
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
  input {
    File count_matrix
    File sample_metadata
    String output_prefix
    String reference_tissue
    String target_tissue

    String memory
    Int cpu
    String disk
    Int preemptible
  }

  command <<<!
    echo '
    args <- commandArgs(trailingOnly = TRUE)
    count_file <- args[1]
    meta_file <- args[2]
    prefix <- args[3]
    reference <- args[4]
    target <- args[5]

    library(DESeq2)
    library(ggplot2)

    counts <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t", skip = 2, check.names = FALSE)
    counts <- counts[, !colnames(counts) %in% c("Name", "Description")]

    meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t")
    meta <- meta[rownames(meta) %in% colnames(counts), , drop = FALSE]
    counts <- counts[, rownames(meta)]

    meta$condition <- factor(meta$condition)
    meta$condition <- relevel(meta$condition, ref = reference)

    dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
    dds <- DESeq(dds)

    res <- results(dds, contrast = c("condition", target, reference))
    write.table(as.data.frame(res), file = paste0(prefix, "_", target, "_vs_", reference, "_deseq2_results.tsv"), sep = "\t", quote = FALSE)

    norm_counts <- counts(dds, normalized = TRUE)
    write.table(as.data.frame(norm_counts), file = paste0(prefix, "_normalized_counts.tsv"), sep = "\t", quote = FALSE)

    tryCatch({
      vsd <- vst(dds)
      pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
      percentVar <- round(100 * attr(pcaData, "percentVar"))

      p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        theme_minimal()

      ggsave(filename = paste0(prefix, "_pca_plot.pdf"), plot = p)
    }, error = function(e) {
      cat("Warning: PCA plot failed but continuing. Error message:\n")
      cat(e$message, "\n")
    })
    ' > run_deseq2.R

    Rscript run_deseq2.R ~{count_matrix} ~{sample_metadata} ~{output_prefix} ~{reference_tissue} ~{target_tissue} > deseq2.log 2>&1
  >>>

  output {
    File results_table = "~{output_prefix}_~{target_tissue}_vs_~{reference_tissue}_deseq2_results.tsv"
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
    version: "1.2.0"
    description: "Runs DESeq2 pairwise tissue comparison with donor adjustment support."
  }
}

