version 1.0

workflow FastaToBamWorkflow {
    input {
        File fasta
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Int cpu = 2
        String memory = "8G"
        String disks = "local-disk 100 HDD"
        Int? num_threads
        Int? num_sort_threads
    }

    String prefix = basename(fasta, ".fasta")

    call FastaToBam {
        input:
            fasta = fasta,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            prefix = prefix,
            cpu = cpu,
            memory = memory,
            disks = disks,
            num_threads = select_first([num_threads, 16]),
            num_sort_threads = select_first([num_sort_threads, 4])
    }

    output {
        File bam = FastaToBam.bam
        File bai = FastaToBam.bai
    }
}

task FastaToBam {
    input {
        File fasta
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String prefix
        Int cpu = 2
        String memory = "8G"
        String disks = "local-disk 100 HDD"
        Int num_threads = 16
        Int num_sort_threads = 4
    }

    command <<<
        set -euxo pipefail

        # Check inputs
        ls -lh ~{fasta}
        ls -lh ~{ref_fasta}
        ls -lh ~{ref_fasta_index}
        ls -lh ~{ref_dict}

        # Align FASTA to reference using minimap2 and convert to BAM using samtools
        minimap2 -ayYL --MD -eqx -x map-hifi -t~{num_threads} ~{ref_fasta} ~{fasta} | \
        samtools sort -@~{num_sort_threads} --no-PG -o ~{prefix}.bam

        # Verify BAM file
        ls -lh ~{prefix}.bam

        # Index the BAM file
        samtools index ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
        cpu: cpu
        memory: memory
        disks: disks
    }
}
