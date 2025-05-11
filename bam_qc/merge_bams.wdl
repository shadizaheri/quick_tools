version 1.0

task merge_crams {
  input {
    File cram_1
    File cram_2
    File reference_fasta
    File reference_fasta_index
    String output_name = "merged.cram"
    Int cpu = 2
    String memory = "8G"
    String disks = "local-disk 100 HDD"
  }

  command <<<
    set -e
    samtools merge -O cram -@ ~{cpu} -R ~{reference_fasta} -o ~{output_name} ~{cram_1} ~{cram_2}
  >>>

  output {
    File merged_cram = output_name
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.3"
    cpu: cpu
    memory: memory
    disks: disks
  }
}

workflow merge_two_crams_workflow {
  input {
    File cram_1
    File cram_2
    File reference_fasta
    File reference_fasta_index
  }

  call merge_crams {
    input:
      cram_1 = cram_1,
      cram_2 = cram_2,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index
  }

  output {
    File final_merged_cram = merge_crams.merged_cram
  }
}
