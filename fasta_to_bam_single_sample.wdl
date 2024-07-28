version 1.0

version 1.0

workflow FastqToBamWorkflow {
    input {
        File fastq
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        Int cpu = 2
        String memory = "16G"
        String disks = "local-disk 200 HDD"
        Int? num_threads
        Int? num_sort_threads
        String sample_id
        String read_group_id = "rg1"
        String read_group_lb = sample_id
        String read_group_pl = "ILLUMINA"
        String read_group_sm = sample_id
    }

    String prefix = basename(fastq, ".fastq.gz")

    call FastqToBam {
        input:
            fastq = fastq,
            ref_fasta = ref_fasta,
            ref_fasta_index = ref_fasta_index,
            ref_dict = ref_dict,
            prefix = prefix,
            cpu = cpu,
            memory = memory,
            disks = disks,
            num_threads = select_first([num_threads, 16]),
            num_sort_threads = select_first([num_sort_threads, 4]),
            read_group_id = read_group_id,
            read_group_lb = read_group_lb,
            read_group_pl = read_group_pl,
            read_group_sm = read_group_sm
    }

    output {
        File bam = FastqToBam.bam
        File bai = FastqToBam.bai
    }
}


task FastqToBam {
    input {
        File fastq
        File ref_fasta
        File ref_fasta_index
        File ref_dict
        String prefix
        Int cpu = 2
        String memory = "16G"
        String disks = "local-disk 200 HDD"
        Int num_threads = 16
        Int num_sort_threads = 4
        String read_group_id
        String read_group_lb
        String read_group_pl
        String read_group_sm
    }

    command <<<
        set -euxo pipefail

        # Check inputs
        ls -lh ~{fastq}
        ls -lh ~{ref_fasta}
        ls -lh ~{ref_fasta_index}
        ls -lh ~{ref_dict}

        # Align FASTQ to reference using minimap2 and save to temporary file
        minimap2 -ayYL --MD -eqx -x map-hifi -t~{num_threads} ~{ref_fasta} ~{fastq} > ~{prefix}.sam 2> ~{prefix}.minimap2.log

        # Check if minimap2 succeeded
        if [ $? -ne 0 ]; then
            echo "minimap2 failed. Check the log file ~{prefix}.minimap2.log for details." >&2
            cat ~{prefix}.minimap2.log
            exit 1
        fi

        head -n 50 ~{prefix}.sam

        samtools addreplacerg -r "ID:~{read_group_id}" -r "LB:~{read_group_lb}" -r "PL:~{read_group_pl}" -r "SM:~{read_group_sm}" -O BAM -o ~{prefix}.rg.bam ~{prefix}.sam

        samtools sort -@~{num_sort_threads} --no-PG -o ~{prefix}.sorted.bam ~{prefix}.rg.bam

        ls -lh ~{prefix}.sorted.bam

        samtools index ~{prefix}.sorted.bam
    >>>

    output {
        File bam = "~{prefix}.sorted.bam"
        File bai = "~{prefix}.sorted.bam.bai"
        File minimap2_log = "~{prefix}.minimap2.log"
        File sam = "~{prefix}.sam"
        File rg_bam = "~{prefix}.rg.bam"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
        cpu: cpu
        memory: memory
        disks: disks
    }
}
