version 1.0

## VerifyBamID2 Workflow for Short Read CRAM Files
## This workflow runs VerifyBamID2 contamination estimation on CRAM files

workflow VerifyBamID2Workflow {
    input {
        File input_cram
        File input_crai
        File reference_fasta
        File reference_fasta_index
        File? reference_dict
        File svd_prefix_bed
        File svd_prefix_mu
        File svd_prefix_UD
        String sample_name
        
        # Optional parameters
        Float? contamination_threshold = 0.02
        Int? max_depth = 200
        Int? min_mapq = 20
        Int? min_baseq = 20
        String? extra_args = ""
        
        # Runtime parameters
        String docker_image = "griffan/verifybamid2:2.0.1"
        Int cpu = 4
        Int memory_gb = 8
        Int disk_size_gb = 100
        Int preemptible_tries = 2
    }

    call VerifyBamID2 {
        input:
            input_cram = input_cram,
            input_crai = input_crai,
            reference_fasta = reference_fasta,
            reference_fasta_index = reference_fasta_index,
            reference_dict = reference_dict,
            svd_prefix_bed = svd_prefix_bed,
            svd_prefix_mu = svd_prefix_mu,
            svd_prefix_UD = svd_prefix_UD,
            sample_name = sample_name,
            contamination_threshold = contamination_threshold,
            max_depth = max_depth,
            min_mapq = min_mapq,
            min_baseq = min_baseq,
            extra_args = extra_args,
            docker_image = docker_image,
            cpu = cpu,
            memory_gb = memory_gb,
            disk_size_gb = disk_size_gb,
            preemptible_tries = preemptible_tries
    }

    output {
        File selfSM = VerifyBamID2.selfSM
        File ancestry = VerifyBamID2.ancestry
        File depthSM = VerifyBamID2.depthSM
        Float contamination = VerifyBamID2.contamination
        Float contamination_error = VerifyBamID2.contamination_error
        Boolean is_contaminated = VerifyBamID2.is_contaminated
        File log_file = VerifyBamID2.log_file
    }

    meta {
        description: "Run VerifyBamID2 contamination estimation on CRAM files"
        author: "Your Name"
        email: "your.email@example.com"
    }

    parameter_meta {
        input_cram: "Input CRAM file"
        input_crai: "Input CRAM index file (.crai)"
        reference_fasta: "Reference genome FASTA file"
        reference_fasta_index: "Reference genome FASTA index (.fai)"
        reference_dict: "Reference genome dictionary file (.dict)"
        svd_prefix_bed: "VerifyBamID2 SVD resource BED file"
        svd_prefix_mu: "VerifyBamID2 SVD resource MU file"
        svd_prefix_UD: "VerifyBamID2 SVD resource UD file"
        sample_name: "Sample name for output files"
        contamination_threshold: "Contamination threshold for flagging samples"
        max_depth: "Maximum depth for pileup"
        min_mapq: "Minimum mapping quality"
        min_baseq: "Minimum base quality"
        docker_image: "Docker image for VerifyBamID2"
        cpu: "Number of CPUs"
        memory_gb: "Memory in GB"
        disk_size_gb: "Disk size in GB"
    }
}

task VerifyBamID2 {
    input {
        File input_cram
        File input_crai
        File reference_fasta
        File reference_fasta_index
        File? reference_dict
        File svd_prefix_bed
        File svd_prefix_mu
        File svd_prefix_UD
        String sample_name
        
        Float? contamination_threshold
        Int? max_depth
        Int? min_mapq
        Int? min_baseq
        String? extra_args
        
        String docker_image
        Int cpu
        Int memory_gb
        Int disk_size_gb
        Int preemptible_tries
    }

    String base_name = basename(input_cram, ".cram")
    String output_prefix = sample_name + ".verifybamid2"
    
    command <<<
        set -euo pipefail
        
        # Create output directory
        mkdir -p outputs
        
        # Link input files with proper naming
        ln -s ~{input_cram} input.cram
        ln -s ~{input_crai} input.cram.crai
        ln -s ~{reference_fasta} reference.fa
        ln -s ~{reference_fasta_index} reference.fa.fai
        
        # Link dictionary if provided
        ~{if defined(reference_dict) then "ln -s " + reference_dict + " reference.dict" else ""}
        
        # Get the base name for SVD files (remove extension)
        SVD_PREFIX=$(basename ~{svd_prefix_bed} .bed)
        
        # Link SVD resource files
        ln -s ~{svd_prefix_bed} ${SVD_PREFIX}.bed
        ln -s ~{svd_prefix_mu} ${SVD_PREFIX}.mu
        ln -s ~{svd_prefix_UD} ${SVD_PREFIX}.UD
        
        # Run VerifyBamID2
        VerifyBamID \
            --Verbose \
            --NumPC 4 \
            --Output outputs/~{output_prefix} \
            --BamFile input.cram \
            --Reference reference.fa \
            --UDPath . \
            --MeanPath . \
            --BedPath . \
            ~{"--MaxDepth " + max_depth} \
            ~{"--MinMapQ " + min_mapq} \
            ~{"--MinQ " + min_baseq} \
            ~{extra_args} \
            --SVDPrefix ${SVD_PREFIX} \
            2>&1 | tee verifybamid2.log
        
        # Extract contamination values
        CONTAMINATION=$(awk 'NR==2 {print $7}' outputs/~{output_prefix}.selfSM)
        CONTAMINATION_ERROR=$(awk 'NR==2 {print $8}' outputs/~{output_prefix}.selfSM)
        
        # Check if contamination exceeds threshold
        THRESHOLD=~{default="0.02" contamination_threshold}
        IS_CONTAMINATED=$(python3 -c "print('true' if float('${CONTAMINATION}') > float('${THRESHOLD}') else 'false')")
        
        # Write results to files for output
        echo ${CONTAMINATION} > contamination.txt
        echo ${CONTAMINATION_ERROR} > contamination_error.txt
        echo ${IS_CONTAMINATED} > is_contaminated.txt
        
        # Log summary
        echo "Sample: ~{sample_name}" >> verifybamid2.log
        echo "Contamination: ${CONTAMINATION}" >> verifybamid2.log
        echo "Contamination Error: ${CONTAMINATION_ERROR}" >> verifybamid2.log
        echo "Is Contaminated (>${THRESHOLD}): ${IS_CONTAMINATED}" >> verifybamid2.log
    >>>

    runtime {
        docker: docker_image
        cpu: cpu
        memory: memory_gb + " GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: preemptible_tries
        bootDiskSizeGb: 15
    }

    output {
        File selfSM = "outputs/" + output_prefix + ".selfSM"
        File ancestry = "outputs/" + output_prefix + ".ancestry"
        File depthSM = "outputs/" + output_prefix + ".depthSM"
        Float contamination = read_float("contamination.txt")
        Float contamination_error = read_float("contamination_error.txt")
        Boolean is_contaminated = read_boolean("is_contaminated.txt")
        File log_file = "verifybamid2.log"
    }

    parameter_meta {
        input_cram: "Input CRAM file"
        input_crai: "CRAM index file"
        reference_fasta: "Reference genome FASTA"
        reference_fasta_index: "Reference FASTA index"
        svd_prefix_bed: "SVD BED file for VerifyBamID2"
        svd_prefix_mu: "SVD MU file for VerifyBamID2"
        svd_prefix_UD: "SVD UD file for VerifyBamID2"
        sample_name: "Sample identifier"
        contamination_threshold: "Threshold for contamination flagging"
        max_depth: "Maximum read depth for analysis"
        min_mapq: "Minimum mapping quality filter"
        min_baseq: "Minimum base quality filter"
        docker_image: "Container image for VerifyBamID2"
        cpu: "CPU cores"
        memory_gb: "Memory allocation in GB"
        disk_size_gb: "Disk space in GB"
    }
}

# Workflow for processing multiple samples
workflow VerifyBamID2BatchWorkflow {
    input {
        Array[File] input_crams
        Array[File] input_crais
        Array[String] sample_names
        File reference_fasta
        File reference_fasta_index
        File? reference_dict
        File svd_prefix_bed
        File svd_prefix_mu
        File svd_prefix_UD
        
        Float? contamination_threshold = 0.02
        Int? max_depth = 200
        Int? min_mapq = 20
        Int? min_baseq = 20
        String? extra_args = ""
        
        String docker_image = "griffan/verifybamid2:2.0.1"
        Int cpu = 4
        Int memory_gb = 8
        Int disk_size_gb = 100
        Int preemptible_tries = 2
    }

    # Process each sample
    scatter (i in range(length(input_crams))) {
        call VerifyBamID2 {
            input:
                input_cram = input_crams[i],
                input_crai = input_crais[i],
                reference_fasta = reference_fasta,
                reference_fasta_index = reference_fasta_index,
                reference_dict = reference_dict,
                svd_prefix_bed = svd_prefix_bed,
                svd_prefix_mu = svd_prefix_mu,
                svd_prefix_UD = svd_prefix_UD,
                sample_name = sample_names[i],
                contamination_threshold = contamination_threshold,
                max_depth = max_depth,
                min_mapq = min_mapq,
                min_baseq = min_baseq,
                extra_args = extra_args,
                docker_image = docker_image,
                cpu = cpu,
                memory_gb = memory_gb,
                disk_size_gb = disk_size_gb,
                preemptible_tries = preemptible_tries
        }
    }

    # Aggregate results
    call CollectResults {
        input:
            selfSM_files = VerifyBamID2.selfSM,
            sample_names = sample_names,
            contamination_values = VerifyBamID2.contamination,
            contamination_errors = VerifyBamID2.contamination_error,
            is_contaminated_flags = VerifyBamID2.is_contaminated,
            contamination_threshold = contamination_threshold
    }

    output {
        Array[File] selfSM_files = VerifyBamID2.selfSM
        Array[File] ancestry_files = VerifyBamID2.ancestry
        Array[File] depthSM_files = VerifyBamID2.depthSM
        Array[Float] contamination_values = VerifyBamID2.contamination
        Array[Boolean] is_contaminated_flags = VerifyBamID2.is_contaminated
        File summary_report = CollectResults.summary_report
        File contaminated_samples = CollectResults.contaminated_samples
    }
}

task CollectResults {
    input {
        Array[File] selfSM_files
        Array[String] sample_names
        Array[Float] contamination_values
        Array[Float] contamination_errors
        Array[Boolean] is_contaminated_flags
        Float? contamination_threshold
    }

    command <<<
        set -euo pipefail
        
        # Create summary report
        echo -e "Sample\tContamination\tContamination_Error\tIs_Contaminated\tThreshold" > summary_report.txt
        
        # Process arrays
        SAMPLES=(~{sep=" " sample_names})
        CONTAMINATION=(~{sep=" " contamination_values})
        CONTAMINATION_ERROR=(~{sep=" " contamination_errors})
        IS_CONTAMINATED=(~{sep=" " is_contaminated_flags})
        THRESHOLD=~{default="0.02" contamination_threshold}
        
        # Create contaminated samples list
        echo "Contaminated Samples (threshold > ${THRESHOLD}):" > contaminated_samples.txt
        echo "Sample\tContamination" >> contaminated_samples.txt
        
        for i in "${!SAMPLES[@]}"; do
            echo -e "${SAMPLES[$i]}\t${CONTAMINATION[$i]}\t${CONTAMINATION_ERROR[$i]}\t${IS_CONTAMINATED[$i]}\t${THRESHOLD}" >> summary_report.txt
            
            if [ "${IS_CONTAMINATED[$i]}" = "true" ]; then
                echo -e "${SAMPLES[$i]}\t${CONTAMINATION[$i]}" >> contaminated_samples.txt
            fi
        done
        
        # Count contaminated samples
        CONTAMINATED_COUNT=$(grep -c "true" summary_report.txt || echo "0")
        TOTAL_COUNT=${#SAMPLES[@]}
        
        echo "" >> contaminated_samples.txt
        echo "Summary: ${CONTAMINATED_COUNT}/${TOTAL_COUNT} samples contaminated" >> contaminated_samples.txt
    >>>

    runtime {
        docker: "ubuntu:20.04"
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 HDD"
    }

    output {
        File summary_report = "summary_report.txt"
        File contaminated_samples = "contaminated_samples.txt"
    }
}
