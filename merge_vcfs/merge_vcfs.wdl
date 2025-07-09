version 1.0

workflow MergeVCFs {
    meta {
        description: "Merge multiple VCF files using bcftools merge"
        author: "Shadi Zaheri"
        version: "1.0"
    }
    
    input {
        Array[File] vcf_files
        Array[File]? vcf_indices  # corresponding .tbi or .csi index files
        String output_prefix = "merged"
        String merge_strategy = "both"  # "snps", "indels", "both", "all", "none"
        Boolean missing_to_ref = false
        Boolean force_samples = false
        String? regions_file
        String? regions_string
        Int threads = 4
        String bcftools_docker = "quay.io/biocontainers/bcftools:1.17--haef29d1_0"
        String memory = "8 GB"
        String disk_space = "50 GB"
        Int preemptible_tries = 2
    }
    
    # Index VCF files if indices are not provided
    if (!defined(vcf_indices)) {
        scatter (vcf in vcf_files) {
            call IndexVCF {
                input:
                    vcf_file = vcf,
                    docker = bcftools_docker,
                    memory = memory,
                    disk_space = disk_space,
                    preemptible_tries = preemptible_tries
            }
        }
    }
    
    # Use provided indices or generated ones
    Array[File] final_indices = select_first([vcf_indices, IndexVCF.index_file])
    
    call MergeVCF {
        input:
            vcf_files = vcf_files,
            vcf_indices = final_indices,
            output_prefix = output_prefix,
            merge_strategy = merge_strategy,
            missing_to_ref = missing_to_ref,
            force_samples = force_samples,
            regions_file = regions_file,
            regions_string = regions_string,
            threads = threads,
            docker = bcftools_docker,
            memory = memory,
            disk_space = disk_space,
            preemptible_tries = preemptible_tries
    }
    
    output {
        File merged_vcf = MergeVCF.merged_vcf
        File merged_vcf_index = MergeVCF.merged_vcf_index
        File merge_stats = MergeVCF.merge_stats
    }
}

task IndexVCF {
    input {
        File vcf_file
        String docker
        String memory
        String disk_space
        Int preemptible_tries
    }
    
    String base_name = basename(vcf_file, ".vcf.gz")
    
    command <<<
        set -euo pipefail
        
        # Create index
        bcftools index -t ~{vcf_file}
    >>>
    
    runtime {
        docker: docker
        memory: memory
        disk: disk_space
        preemptible: preemptible_tries
    }
    
    output {
        File index_file = "~{vcf_file}.tbi"
    }
}

task MergeVCF {
    input {
        Array[File] vcf_files
        Array[File] vcf_indices
        String output_prefix
        String merge_strategy
        Boolean missing_to_ref
        Boolean force_samples
        File? regions_file
        String? regions_string
        Int threads
        String docker
        String memory
        String disk_space
        Int preemptible_tries
    }
    
    command <<<
        set -euo pipefail
        
        # Create file list for bcftools merge
        printf '%s\n' ~{sep=' ' vcf_files} > vcf_list.txt
        
        # Build bcftools merge command
        MERGE_CMD="bcftools merge"
        MERGE_CMD="$MERGE_CMD --threads ~{threads}"
        MERGE_CMD="$MERGE_CMD --merge ~{merge_strategy}"
        
        # Add optional parameters
        ~{if missing_to_ref then "MERGE_CMD=\"$MERGE_CMD --missing-to-ref\"" else ""}
        ~{if force_samples then "MERGE_CMD=\"$MERGE_CMD --force-samples\"" else ""}
        ~{if defined(regions_file) then "MERGE_CMD=\"$MERGE_CMD --regions-file ~{regions_file}\"" else ""}
        ~{if defined(regions_string) then "MERGE_CMD=\"$MERGE_CMD --regions '~{regions_string}'\"" else ""}
        
        # Output options
        MERGE_CMD="$MERGE_CMD --output-type z"
        MERGE_CMD="$MERGE_CMD --output ~{output_prefix}.vcf.gz"
        
        # Add input files
        MERGE_CMD="$MERGE_CMD --file-list vcf_list.txt"
        
        # Run merge command
        echo "Running: $MERGE_CMD"
        eval $MERGE_CMD
        
        # Index the merged VCF
        bcftools index -t ~{output_prefix}.vcf.gz
        
        # Generate merge statistics
        bcftools stats ~{output_prefix}.vcf.gz > ~{output_prefix}_stats.txt
        
        # Log file information
        echo "Merged VCF info:"
        bcftools query -l ~{output_prefix}.vcf.gz | wc -l | awk '{print "Number of samples: " $1}'
        bcftools view -H ~{output_prefix}.vcf.gz | wc -l | awk '{print "Number of variants: " $1}'
    >>>
    
    runtime {
        docker: docker
        memory: memory
        disk: disk_space
        cpu: threads
        preemptible: preemptible_tries
    }
    
    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_vcf_index = "~{output_prefix}.vcf.gz.tbi"
        File merge_stats = "~{output_prefix}_stats.txt"
    }
}