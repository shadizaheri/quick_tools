task ExtractAndFilterSampleVCF {
    input {
        File joint_vcf
        String sample_name
    }

    command <<<
        set -euxo pipefail

        bcftools index ~{joint_vcf}

        bcftools view -s ~{sample_name} ~{joint_vcf} -o ~{sample_name}.subset.g.vcf.gz

        bcftools view -s ~{sample_name} -i 'GT="0/0"' ~{joint_vcf} -o ~{sample_name}.homref.g.vcf.gz

        tabix -p vcf ~{sample_name}.subset.g.vcf.gz
        tabix -p vcf ~{sample_name}.homref.g.vcf.gz
    >>>

    output {
        File single_sample_vcf = "${sample_name}.subset.g.vcf.gz"
        File single_sample_vcf_tbi = "${sample_name}.subset.g.vcf.gz.tbi"
        File homref_sample_vcf = "${sample_name}.homref.g.vcf.gz"
        File homref_sample_vcf_tbi = "${sample_name}.homref.g.vcf.gz.tbi"
    }

    runtime {
        cpu: 1
        memory: "64 GiB"
        disks: "local-disk ${disk_size} HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 1
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
    }

    Int disk_size = ceil(2 * size(joint_vcf, "GiB")) + 1
}
