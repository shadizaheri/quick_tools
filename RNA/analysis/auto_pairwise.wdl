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

  call RunAllPairwiseDESeq2 {
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
    Array[File] all_result_tables = RunAllPairwiseDESeq2.all_result_tables
    File normalized_counts = RunAllPairwiseDESeq2.normalized_counts
    File pca_plot = RunAllPairwiseDESeq2.pca_plot
    File log_file = RunAllPairwiseDESeq2.log_file
  }
}

task RunAllPairwiseDESeq2 {
  input {
    File count_matrix
    File sample_metadata
    String output_prefix

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

    library(DESeq2)
    library(ggplot2)
    library(tools)

    counts <- read.table(gzfile(count_file), header = TRUE, row.names = 1, sep = "\t", skip = 2, check.names = FALSE)
    counts <- counts[, !colnames(counts) %in% c("Name", "Description")]

    meta <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t")
    meta <- meta[rownames(meta) %in% colnames(counts), , drop = FALSE]
    counts <- counts[, rownames(meta)]

    meta$condition <- factor(meta$condition)
    dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ condition)
    dds <- DESeq(dds)

    tissues <- levels(meta$condition)
    result_files <- c()

    for (i in 1:(length(tissues)-1)) {
      for (j in (i+1):length(tissues)) {
        ref <- tissues[i]
        tgt <- tissues[j]
        res <- results(dds, contrast = c("condition", tgt, ref))
        fname <- paste0(prefix, "_", tgt, "_vs_", ref, "_deseq2_results.tsv")
        write.table(as.data.frame(res), file = fname, sep = "\t", quote = FALSE)
        result_files <- c(result_files, fname)
      }
    }

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

    Rscript run_deseq2.R ~{count_matrix} ~{sample_metadata} ~{output_prefix} > deseq2.log 2>&1
  >>>

  output {
    Array[File] all_result_tables = glob("*deseq2_results.tsv")
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
    version: "1.3.0"
    description: "Runs all pairwise DESeq2 comparisons across tissues using condition column. Outputs all DE tables, normalized counts, and PCA."
  }
}

