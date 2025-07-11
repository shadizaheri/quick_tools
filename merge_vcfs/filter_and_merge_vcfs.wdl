version 1.0

workflow FilterAndMergeVCFs {
  input {
    File sample_list_tsv
    String tissue_name
    String gcs_input_prefix  # e.g. gs://bucket/samples
    String gcs_output_prefix # e.g. gs://bucket/outputtedsamples
  }

  Array[String] sample_ids = read_lines(sample_list_tsv)

  scatter (sample_id in sample_ids) {
    call ExtractAndFilterVCF {
      input:
        sample_id = sample_id,
        gcs_input_prefix = gcs_input_prefix
    }

    call FixPloidyAndBgzip {
      input:
        filtered_vcf = ExtractAndFilterVCF.filtered_vcf
    }

    call UploadToGCS {
      input:
        file_to_upload = FixPloidyAndBgzip.bgzipped_vcf,
        sample_id = sample_id,
        gcs_output_prefix = gcs_output_prefix
    }
  }

  call MergeAndIndexVCFs {
    input:
      vcfs_to_merge = FixPloidyAndBgzip.bgzipped_vcf,
      output_name = tissue_name + "_merged.vcf.gz"
  }

  call UploadMergedToGCS {
    input:
      merged_vcf = MergeAndIndexVCFs.merged_vcf,
      merged_vcf_index = MergeAndIndexVCFs.merged_vcf_index,
      gcs_output_prefix = gcs_output_prefix,
      tissue_name = tissue_name
  }

  output {
    File merged_vcf = MergeAndIndexVCFs.merged_vcf
    File merged_vcf_index = MergeAndIndexVCFs.merged_vcf_index
    String merged_vcf_gcs_path = UploadMergedToGCS.merged_vcf_gcs_path
  }
}

task ExtractAndFilterVCF {
  input {
    String sample_id
    String gcs_input_prefix
  }

  command <<<
    bcftools view \
      -i 'FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) > 0.35 && FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) < 0.65' \
      ${gcs_input_prefix}/${sample_id}/${sample_id}.hard-filtered.vcf.gz -Oz -o ${sample_id}_HET.vcf.gz
    tabix -p vcf ${sample_id}_HET.vcf.gz
  >>>

  output {
    File filtered_vcf = "${sample_id}_HET.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    memory: "2G"
    cpu: 1
  }
}

task FixPloidyAndBgzip {
  input {
    File filtered_vcf
  }

  command <<<
    bcftools +fixploidy ${filtered_vcf} -o fixed.vcf
    bgzip -c fixed.vcf > fixed_HET_FP.vcf.gz
    tabix -p vcf fixed_HET_FP.vcf.gz
  >>>

  output {
    File bgzipped_vcf = "fixed_HET_FP.vcf.gz"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    memory: "2G"
    cpu: 1
  }
}

task UploadToGCS {
  input {
    File file_to_upload
    String sample_id
    String gcs_output_prefix
  }

  command <<<
    gsutil cp ${file_to_upload} ${gcs_output_prefix}/${sample_id}_HET_FP.vcf.gz
    gsutil cp ${file_to_upload}.tbi ${gcs_output_prefix}/${sample_id}_HET_FP.vcf.gz.tbi
  >>>

  output {
    File uploaded_vcf = "${file_to_upload}"
  }

  runtime {
    docker: "google/cloud-sdk:slim"
    memory: "1G"
    cpu: 1
  }
}

task MergeAndIndexVCFs {
  input {
    Array[File] vcfs_to_merge
    String output_name
  }

  command <<<
    bcftools merge -Oz -o ${output_name} ${sep=' ' vcfs_to_merge}
    tabix -p vcf ${output_name}
  >>>

  output {
    File merged_vcf = output_name
    File merged_vcf_index = output_name + ".tbi"
  }

  runtime {
    docker: "us.gcr.io/broad-dsp-lrma/lr-align:0.1.28"
    memory: "4G"
    cpu: 2
  }
}

task UploadMergedToGCS {
  input {
    File merged_vcf
    File merged_vcf_index
    String gcs_output_prefix
    String tissue_name
  }

  command <<<
    gsutil cp ${merged_vcf} ${gcs_output_prefix}/${tissue_name}_merged.vcf.gz
    gsutil cp ${merged_vcf_index} ${gcs_output_prefix}/${tissue_name}_merged.vcf.gz.tbi
  >>>

  output {
    String merged_vcf_gcs_path = "${gcs_output_prefix}/${tissue_name}_merged.vcf.gz"
  }

  runtime {
    docker: "google/cloud-sdk:slim"
    memory: "1G"
    cpu: 1
  }
}