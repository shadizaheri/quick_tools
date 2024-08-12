version 1.0
workflow TransformAndSplitVCFWorkflow {
    input {
        File joint_vcf
        Array[String] raw_sample_names  # Directly read from the Terra data table
    }

    # Transformation step
    call transform_sample_names {
        input:
            raw_sample_names = raw_sample_names
    }

    scatter (samplename in transform_sample_names.formatted_sample_names) {
        call ExtractAndFilterSampleVCF {
            input:
                joint_vcf = joint_vcf,
                sample_name = samplename
        }
    }

    output {
        Array[File] single_sample_vcfs = ExtractAndFilterSampleVCF.single_sample_vcf
        Array[File] single_sample_vcf_tbis = ExtractAndFilterSampleVCF.single_sample_vcf_tbi
        Array[File] homref_sample_vcfs = ExtractAndFilterSampleVCF.homref_sample_vcf
        Array[File] homref_sample_vcf_tbis = ExtractAndFilterSampleVCF.homref_sample_vcf_tbi
    }
}

task transform_sample_names {
    input {
        Array[String] raw_sample_names
    }

    command <<<
        echo ~{sep=" " raw_sample_names} > raw_names.txt
        
        awk '{
            gsub("-1","-1", $0);
            gsub("-2","-2", $0);
            gsub("-3","-3", $0);
            print $0;
        }' raw_names.txt > formatted_names.txt
    >>>

    output {
        Array[String] formatted_sample_names = read_lines("formatted_names.txt")
    }
    
    runtime {
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
    }
}

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
