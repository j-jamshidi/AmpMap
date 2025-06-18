run=$1
BASEDIR="/EBSDataDrive/ONT/script"
RUNDIR="/EBSDataDrive/ONT/Runs/${run}"
WORKDIR="/EBSDataDrive/ONT/Runs/${run}/result"

function get_presign() {

	DATID=$(bs -c POWH list datasets --input-biosample $sid --not-type "illumina.fastq.v1.8" --sort-by AppSession.DateCreated --terse | tail -n1)
	BAM=$(bs dataset -c POWH content --id=$DATID --extension=bam --terse)
	BAI=$(bs dataset -c POWH content --id=$DATID --extension=bam.bai --terse)


	AWS_RSA256_link=$(bs -c POWH file link -i "$BAM")
	wget_out=$(wget --save-headers --max-redirect=0 -O - "$AWS_RSA256_link" 2>&1)
	BAMPre=$(echo "$wget_out" | grep -i "Location" | tail -n 1 | awk '{print $2}')

	AWS_RSA256_link=$(bs -c POWH file link -i "$BAI")
	wget_out=$(wget --save-headers --max-redirect=0 -O - "$AWS_RSA256_link" 2>&1)
	BAIPre=$(echo "$wget_out" | grep -i "Location" | tail -n 1 | awk '{print $2}')

	echo $BAMPre $BAIPre

}

while IFS=, read -r Batch Barcode Episode Coordinate Variant1 Variant2 EpisodeWES remainder; do
 {
	OUTXML=${WORKDIR}/$Barcode/${Episode}.xml
	if [ $EpisodeWES == NA ]; 
		then cp ${BASEDIR}/solo_LR.xml ${OUTXML}; 
	else
	        cat /EBSDataDrive/software/sample_ran.txt /EBSDataDrive/software/sample_ran_CRE_BS.txt > /EBSDataDrive/software/sample_SR.txt	
		wesrun=$(cat /EBSDataDrive/software/sample_SR.txt | grep $EpisodeWES | cut -f 6 | tr '[:lower:]' '[:upper:]');
		sid=$(cat /EBSDataDrive/software/sample_SR.txt | grep $EpisodeWES | cut -f 1 | tr '[:lower:]' '[:upper:]'); 
		if [ -z "$wesrun" ]; then cp ${BASEDIR}/solo_LR.xml ${OUTXML}; 
		else 
			cp ${BASEDIR}/solo_LR_SR.xml ${OUTXML};
			read -r BAMpresign BAIpresign <<< $(get_presign);
			bamurl=$(echo $BAMpresign | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
			baiurl=$(echo $BAIpresign | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
			sed -i -e "s/C1_SR_index_URL/$baiurl/g" ${OUTXML};
 			sed -i -e "s/C1_SR_BAM_URL/$bamurl/g" ${OUTXML};
		fi; 
	fi;
	sed -i -e "s/C1_IID/$Episode/g" ${OUTXML};
	raw2=$Variant2
	if [[ "$Variant1" =~ chr ]] || [[ "$Variant2" =~ chr ]]; then
		if [[ ! "$Variant2" =~ chr ]]; then
			Variant2=$Variant1
		fi
		chr1=$(echo $Variant1 | cut -d ':' -f 1)
		pos1=$(echo $Variant1 | cut -d ':' -f 2 | cut -d ' ' -f 1)
		chr2=$(echo $Variant2 | cut -d ':' -f 1)
		pos2=$(echo $Variant2 | cut -d ':' -f 2 | cut -d ' ' -f 1)
		if [ "$chr1" != "$chr2" ]; then
			echo "Variants are on different chromosomes."
			exit 1
		fi
		if [ $pos1 -lt $pos2 ]; then
			start=$((pos1 - 100))
			end=$((pos2 + 100))
		else
			start=$((pos2 - 100))
			end=$((pos1 + 100))
		fi
		Coordinate="${chr1}:${start}-${end}"
	fi
	sed -i -e "s/CHRSTARTEND/$Coordinate/g" ${OUTXML};
	Variant2=$raw2
	if [[ "$Variant2" =~ chr ]]; then
		lrbam=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}_phased.bam --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
		lrbai=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}_phased.bam.bai --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
		phasevcf=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}_Phased.vcf --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
		wfvcf=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}.wf_snp.vcf.gz --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
	else
		lrbam=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}_QC.bam --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
		lrbai=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}_QC.bam.bai --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
		phasevcf=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}.vcf --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
		wfvcf=$(aws s3 presign s3://nswhp-gaia-poc-pl/ONT/$run/$Barcode/${Episode}.wf_snp.vcf.gz --expire=604800 | sed -r 's/\//\\\//g' | sed -r 's/&/\\&amp;/g');
	fi
	sed -i -e "s/WF_SNP_VCF/$wfvcf/g" ${OUTXML};
	sed -i -e "s/PHASED_VCF/$phasevcf/g" ${OUTXML};
	sed -i -e "s/C1_LR_BAM_URL/$lrbam/g" ${OUTXML};
	sed -i -e "s/C1_LR_index_URL/$lrbai/g" ${OUTXML};
    }
done < ${RUNDIR}/${run}.info

