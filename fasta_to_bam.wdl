version 1.0

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
    }

    command <<<
        set -euxo pipefail

        # Align FASTA to reference using minimap2 and convert to BAM using samtools
        minimap2 -ax map-ont ~{ref_fasta} ~{fasta} | \
        samtools view -Sb - > ~{prefix}.bam

        # Index the BAM file
        samtools index ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/mosdepth:sz_v3272024"
        cpu: cpu
        memory: memory
        disks: disks
    }
}

workflow FastaToBamWorkflow {
    input {
        File fasta
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String prefix = "output"
        Int cpu = 2
        String memory = "8G"
        String disks = "local-disk 100 HDD"
    }

    call FastaToBam {
        input:
            fasta = fasta,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            prefix = prefix,
            cpu = cpu,
            memory = memory,
            disks = disks
    }

    output {
        File bam = FastaToBam.bam
        File bai = FastaToBam.bai
    }
}