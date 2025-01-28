#!/bin/bash

# Script: run.sh
# Description: Pipeline for processing and phasing ONT amplicon data.
#works for phasing when two variants are provided, or QC when one or none is provided.
#===============================================================================
# Configuration and Setup
#===============================================================================
RUNID=$1
basedir="/EBSDataDrive/ONT/Runs/${RUNID}"
workdir="${basedir}/result"

# Reference and tool paths
REFERENCE="/EFSGaiaDataDrive/ref/ONT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CLAIR3_PATH="/EBSDataDrive/ONT/Clair3"
HAPCUT2_PATH="/EBSDataDrive/ONT/HapCUT2-1.3.4/build"
SCRIPT_PATH="/EBSDataDrive/ONT/script"

#===============================================================================
# Functions
#===============================================================================
prepare_vcf() {
    local episode=$1
    local barcode=$2
    local variant1=$3
    local variant2=$4

    cp "${SCRIPT_PATH}/dummy.vcf" "${workdir}/${barcode}/${episode}.vcf"
    sed -i -e "s/sample/${episode}/g" "${workdir}/${barcode}/${episode}.vcf"

    # Process first variant
    variant1=$(echo "${variant1}" | tr ':>\t ' '-' | sed 's/--*/-/g')
    IFS='-' read -r -a array1 <<< "${variant1}"
    sed -i -e "s/chrA/${array1[0]}/g" "${workdir}/${barcode}/${episode}.vcf"
    sed -i -e "s/POS1/${array1[1]}/g" "${workdir}/${barcode}/${episode}.vcf"
    sed -i -e "s/REF1/${array1[2]}/g" "${workdir}/${barcode}/${episode}.vcf"
    sed -i -e "s/ALT1/${array1[3]}/g" "${workdir}/${barcode}/${episode}.vcf"
    unset array1

    if [[ "$variant2" =~ chr ]]; then
        # Process second variant
        variant2=$(echo "${variant2}" | tr ':>\t ' '-' | sed 's/--*/-/g')
        IFS='-' read -r -a array2 <<< "${variant2}"
        sed -i -e "s/chrB/${array2[0]}/g" "${workdir}/${barcode}/${episode}.vcf"
        sed -i -e "s/POS2/${array2[1]}/g" "${workdir}/${barcode}/${episode}.vcf"
        sed -i -e "s/REF2/${array2[2]}/g" "${workdir}/${barcode}/${episode}.vcf"
        sed -i -e "s/ALT2/${array2[3]}/g" "${workdir}/${barcode}/${episode}.vcf"
        unset array2

        # Sort VCF file
        grep -v "^#" "${workdir}/${barcode}/${episode}.vcf" | sort -k1,1 -k2,2n > "${workdir}/${barcode}/${episode}_sorted.vcf"
        grep "^#" "${workdir}/${barcode}/${episode}.vcf" > "${workdir}/${barcode}/${episode}_header.vcf"
        cat "${workdir}/${barcode}/${episode}_header.vcf" "${workdir}/${barcode}/${episode}_sorted.vcf" > "${workdir}/${barcode}/${episode}.vcf"
        rm "${workdir}/${barcode}/${episode}_header.vcf" "${workdir}/${barcode}/${episode}_sorted.vcf"
    else
        # Remove second variant line for single variant case
        sed -i '/chrB/d' "${workdir}/${barcode}/${episode}.vcf"
    fi
}

merge_bam_files() {
    local barcode=$1
    local episode=$2
    local coordinate=$3

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Merging bam files for ${barcode}, sample ${episode}..."
    samtools merge -@ 20 -u - "${basedir}/bam_pass/${barcode}"/*bam | \
        samtools sort -@ 20 -o "${workdir}/${barcode}/temp.bam" && \
        samtools index "${workdir}/${barcode}/temp.bam" && \
        samtools view -@ 20 -b "${workdir}/${barcode}/temp.bam" "${coordinate}" > "${workdir}/${barcode}/${episode}.bam" && \
        samtools index "${workdir}/${barcode}/${episode}.bam" && \
        rm "${workdir}/${barcode}/temp"*
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Bam files merged!"
}

run_clair3() {
    local barcode=$1
    local episode=$2

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Running variant calling with Clair3 for ${barcode}, sample ${episode}..."
    "${CLAIR3_PATH}/run_clair3.sh" \
        --bam_fn="${workdir}/${barcode}/${episode}.bam" \
        --bed_fn="${workdir}/${barcode}/${episode}_coordinate.bed" \
        --ref_fn="${REFERENCE}" \
        --threads=6 \
        --platform=ont \
        --model_path="${CLAIR3_PATH}/models/r1041_e82_400bps_sup_v500" \
        --sample_name="${episode}" \
        --use_whatshap_for_final_output_phasing enable \
        --enable_phasing enable \
        --remove_intermediate_dir enable \
        --var_pct_full=1 \
        --ref_pct_full=1 \
        --var_pct_phasing=1 \
        --output="${workdir}/${barcode}/variant_calling_output" \
        --min_coverage=20 >/dev/null 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Clair3 analysis finished!"
}

run_hapcut2() {
    local barcode=$1
    local episode=$2

    echo "$(date '+%Y-%m-%d %H:%M:%S') - Running Hapcut2 for ${barcode}, ${episode}..."
    "${HAPCUT2_PATH}/extractHAIRS" \
        --ont 1 \
        --bam "${workdir}/${barcode}/${episode}.bam" \
        --VCF "${workdir}/${barcode}/${episode}.vcf" \
        --out "${workdir}/${barcode}/fragment_${episode}" \
        --indels 1 \
        --ref "${REFERENCE}" > "${workdir}/${barcode}/HapCUT2.log" 2>&1

    "${HAPCUT2_PATH}/HAPCUT2" \
        --fragments "${workdir}/${barcode}/fragment_${episode}" \
        --VCF "${workdir}/${barcode}/${episode}.vcf" \
        --output "${workdir}/${barcode}/hap2cut_${episode}" >> "${workdir}/${barcode}/HapCUT2.log" 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') - HapCUT2 analysis finished!"
}

#===============================================================================
# Prepare Input File
#===============================================================================
# Remove header from CSV file and create info file
tail -n +2 "${basedir}/sample_sheet.csv" > "${basedir}/temp.info"

# Process the Barcode column to add leading zeros and 'barcode' prefix
# Also ensure clean file ending with newline
awk -F',' 'BEGIN {OFS=","} {
    # Remove any trailing whitespace or special characters
    gsub(/[[:space:]]*$/, "", $0)
    # Modify the barcode in second column
    if ($2 ~ /^[0-9]+$/) {
        $2 = sprintf("barcode%02d", $2)
    }
    print $0
}' "${basedir}/temp.info" > "${basedir}/${RUNID}.info"

# Ensure there's exactly one newline at the end of the file
sed -i -e '$a\' "${basedir}/${RUNID}.info"

# Clean up temporary file
rm "${basedir}/temp.info"

#===============================================================================
# Main Processing Loop
#===============================================================================
while IFS=, read -r Batch Barcode Episode Coordinate Variant1 Variant2 EpisodeWES remainder
do 
    {
        echo "$Barcode"
        mkdir -p "${workdir}/${Barcode}"

        # Clean up Coordinate by removing trailing spaces
        Coordinate=$(echo "${Coordinate}" | sed 's/[[:space:]]*$//')

        # VCF File Preparation
        if [[ ! "$Variant1" =~ chr ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - No variant is provided, continue with basecalling and QC..."
        elif [[ ! "$Variant2" =~ chr ]]; then
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Only one variant is provided, continue with basecalling and QC..."
            prepare_vcf "$Episode" "$Barcode" "$Variant1" ""
        else
            echo "$(date '+%Y-%m-%d %H:%M:%S') - Two variants are provided, continue with basecalling, QC and phasing..."
            prepare_vcf "$Episode" "$Barcode" "$Variant1" "$Variant2"
        fi

        # BAM File Processing
        merge_bam_files "$Barcode" "$Episode" "$Coordinate"

        # Making coordinate bed file
        echo -e "${Coordinate}" | tr -d ' ' | tr ':' '\t' | sed 's/-/\t/g' > "${workdir}/${Barcode}/${Episode}_coordinate.bed"

        # Variant Calling with Clair3
        current_dir=$PWD
        cd "${workdir}/${Barcode}" || exit
        source activate ONT
        run_clair3 "$Barcode" "$Episode"
        conda deactivate

        # Copy output files
        cp "${workdir}/${Barcode}/variant_calling_output/merge_output.vcf.gz" "${Episode}.wf_snp.vcf.gz"
        cp "${workdir}/${Barcode}/variant_calling_output/merge_output.vcf.gz.tbi" "${Episode}.wf_snp.vcf.gz.tbi"

        # HapCUT2 Phasing (only for two variants)
        if [[ "$Variant1" =~ chr ]] && [[ "$Variant2" =~ chr ]]; then
            cd "$current_dir" || exit
            run_hapcut2 "$Barcode" "$Episode"
        fi

        # Final Analysis and Cleanup
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Analysing the reads and writing the results for ${Barcode}, ${Episode}..."
        
        # Choose appropriate Python script based on variant presence
        if [[ ! "$Variant1" =~ chr ]] || [[ ! "$Variant2" =~ chr ]]; then
            python "${SCRIPT_PATH}/basecalling_QC_amplicon.py" \
                "${workdir}/${Barcode}/${Episode}.bam" \
                "${workdir}/${Barcode}/${Episode}_coordinate.bed"
        else
            python "${SCRIPT_PATH}/basecalling_phasing_amplicon.py" \
                "${workdir}/${Barcode}/${Episode}.bam" \
                "${workdir}/${Barcode}/${Episode}.vcf"
        fi
        echo "$(date '+%Y-%m-%d %H:%M:%S') - Main analysis finished!"

        # Cleanup temporary files
        rm -rf "${workdir}/${Barcode}/work"
        if [[ "$Variant1" =~ chr ]] && [[ "$Variant2" =~ chr ]] ; then
            rm "${workdir}/${Barcode}/${Episode}.vcf"
            rm "${workdir}/${Barcode}/hap2cut_${Episode}"
            rm "${workdir}/${Barcode}/fragment_${Episode}"
            rm "${workdir}/${Barcode}/clean-span-hq.bam"*
        fi
    }
done < "${basedir}/${RUNID}.info"