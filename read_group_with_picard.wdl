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
        String read_group_pl = "illumina"
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
        File bam = FastqToBam.final_bam
        File bai = FastqToBam.final_bai
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

        # Output the first few lines of the SAM file for debugging
        head -n 50 ~{prefix}.sam

        # Convert SAM to BAM and sort
        samtools sort -@~{num_sort_threads} --no-PG -o ~{prefix}.sorted.bam ~{prefix}.sam

        # Add read groups using picard
        picard AddOrReplaceReadGroups \
            I=~{prefix}.sorted.bam \
            O=~{prefix}_with_RG.bam \
            RGID=~{read_group_id} \
            RGLB=~{read_group_lb} \
            RGPL=~{read_group_pl} \
            RGSM=~{read_group_sm}

        # Validate the new BAM file
        samtools view -H ~{prefix}_with_RG.bam | grep '^@RG'

        # Index the BAM file
        samtools index ~{prefix}_with_RG.bam
    >>>

    output {
        File final_bam = "~{prefix}_with_RG.bam"
        File final_bai = "~{prefix}_with_RG.bam.bai"
        File minimap2_log = "~{prefix}.minimap2.log"
        File sam = "~{prefix}.sam"
        File sorted_bam = "~{prefix}.sorted.bam"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
        cpu: cpu
        memory: memory
        disks: disks
    }
}
