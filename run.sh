#!/bin/bash
RUNID=$1
basedir=/EBSDataDrive/ONT/Runs/${RUNID}
workdir=${basedir}/result
# Remove header from CSV file
tail -n +2 ${basedir}/sample_sheet.csv >  ${basedir}/${RUNID}.info

# Read the CSV file line by line
while IFS=, read -r Batch Barcode Episode Family Initial Ampliconsize Coordinate Gene Variant1 Variant2 Distance EpisodeWES remainder
do 
    {
        echo $Barcode
	mkdir -p ${workdir}/${Barcode}
        # Create dummy VCF file for this barcode/episode
        cp /EBSDataDrive/ONT/script/dummy.vcf ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/sample/${Episode}/g" ${workdir}/${Barcode}/${Episode}.vcf
        
        # Convert : > and space to - in Variant1 and Variant2
        Variant1=$(echo ${Variant1} | tr ':> ' '-')
        Variant2=$(echo ${Variant2} | tr ':> ' '-')
        
        IFS='-' read -r -a array1 <<< ${Variant1}
        sed -i -e "s/chrA/${array1[0]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/POS1/${array1[1]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/REF1/${array1[2]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/ALT1/${array1[3]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        unset array1

        IFS='-' read -r -a array2 <<< ${Variant2}
        sed -i -e "s/chrB/${array2[0]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/POS2/${array2[1]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/REF2/${array2[2]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        sed -i -e "s/ALT2/${array2[3]}/g" ${workdir}/${Barcode}/${Episode}.vcf
        unset array2

        # Sort the VCF file after replacements
        grep -v "^#" ${workdir}/${Barcode}/${Episode}.vcf | sort -k1,1 -k2,2n > ${workdir}/${Barcode}/${Episode}_sorted.vcf
        grep "^#" ${workdir}/${Barcode}/${Episode}.vcf > ${workdir}/${Barcode}/${Episode}_header.vcf
        cat ${workdir}/${Barcode}/${Episode}_header.vcf ${workdir}/${Barcode}/${Episode}_sorted.vcf > ${workdir}/${Barcode}/${Episode}.vcf
        rm ${workdir}/${Barcode}/${Episode}_header.vcf ${workdir}/${Barcode}/${Episode}_sorted.vcf

        # Merge BAM files
        samtools merge -@ 20 -u - ${basedir}/bam_pass/${Barcode}/*bam | samtools sort -@ 20 -o ${workdir}/${Barcode}/temp.bam && 
        samtools index ${workdir}/${Barcode}/temp.bam && 
        samtools view -@ 20 -b ${workdir}/${Barcode}/temp.bam ${Coordinate} > ${workdir}/${Barcode}/${Episode}.bam && 
        samtools index ${workdir}/${Barcode}/${Episode}.bam && 
        rm ${workdir}/${Barcode}/temp*
        
        # Create coordinate BED file for this barcode/episode
        echo -e "${Coordinate}" | tr ':' '\t' | sed 's/-/\t/g' > ${workdir}/${Barcode}/${Episode}_coordinate.bed

    python /EBSDataDrive/ONT/script/phasing_automation_basecalling.py ${workdir}/${Barcode}/${Episode}.bam ${workdir}/${Barcode}/${Episode}.vcf
	mv ${workdir}/${Barcode}/variant_calling_output/${Episode}.wf_snp.vcf.gz* ${workdir}/${Barcode}/
	rm -rf ${workdir}/${Barcode}/work
	rm ${workdir}/${Barcode}/clean-span-hq.bam*
        rm ${workdir}/${Barcode}/${Episode}.bam*
	rm ${workdir}/${Barcode}/${Episode}.vcf
	rm -rf ${workdir}/${Barcode}/.nextflow
	rm ${workdir}/${Barcode}/.nextflow.log
    }
done <  ${basedir}/${RUNID}.info

