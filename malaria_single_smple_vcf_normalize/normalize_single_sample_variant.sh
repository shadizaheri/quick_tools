#!/bin/bash

# Define paths
bcftools_path="bcftools"
bgzip_path="bgzip"
reference_fa="PlasmoDB-61_Pfalciparum3D7_Genome.fasta"
input_vcf="PG0015-C.recalibrated.pass.merged.vcf.gz"

# Output file names
vcf_no_hapcomp="PG0015-C.recalibrated.pass.merged.no_hapcomp.vcf.gz"
vcf_no_hapdom="PG0015-C.recalibrated.pass.merged.no_hapdom.vcf.gz"
vcf_normalized="PG0015-C.recalibrated.pass.merged.norm.vcf.gz"

# Step 1: Remove the HAPCOMP field
echo "Removing HAPCOMP field..."
$bcftools_path annotate -x 'INFO/HAPCOMP' $input_vcf | $bgzip_path -c > $vcf_no_hapcomp

# Step 2: Remove the HAPDOM field
echo "Removing HAPDOM field..."
$bcftools_path annotate -x 'INFO/HAPDOM' $vcf_no_hapcomp | $bgzip_path -c > $vcf_no_hapdom

# Step 3: Index the modified VCF
echo "Indexing the VCF without HAPCOMP and HAPDOM..."
$bcftools_path index $vcf_no_hapdom

# Step 4: Normalize the VCF
echo "Normalizing the VCF..."
$bcftools_path norm -m -any --atom-overlaps . -f $reference_fa $vcf_no_hapdom | $bgzip_path -c > $vcf_normalized

echo "Normalization completed. Output saved to $vcf_normalized."


