#!/bin/bash

## USAGE -------------------------------------------------------------------

set -e -o pipefail
if [ $# -lt 2 ]; then
    echo "About:   Annotate a VCF file"
    echo "Usage:   annotate_biome_vcfs.sh <full path to input VCF> <output prefix with fullpath>"
    exit 1
fi

## MODULES -------------------------------------------------------------------

module purge
module load vep/96
module load CPAN/5.28.1 perl5
module load mysql
module load bcftools java picard tabix

## SETUP -------------------------------------------------------------------

cachedir=/hpc/packages/minerva-centos7/vep/96/Cache
lofteedir=/hpc/packages/minerva-centos7/vep/96/loftee/
lofteeres=/sc/arion/projects/psychgen2/resources/loftee_resources
INPUT=/sc/hydra/projects/rg_psychgen/scratch/SINAI_Freeze_Two.GL.pVCF.PASS.onTarget.biallelic.chrX.vcf.gz
PERL5LIB=${PERL5LIB}:${lofteedir}
vcf="$1"
out="$2"
dir=`dirname ${out}`
ref="/sc/arion/projects/psychgen2/resources/GRCh38.primary_assembly.genome.fa"
dbnsfp="/sc/arion/projects/psychgen2/resources/dbNSFP4.0a/dbNSFPv4.0a_custombuild.gz"
snpeff="/hpc/packages/minerva-common/snpeff/4.3/snpEff"
dbnsfp_fields="$(echo MPC_score,gnomAD_exomes{,_AFR,_AMR,_ASJ,_EAS,_FIN,_NFE,_SAS}_{AC,AN,AF},gnomAD_genomes{,_AFR,_AMR,_ASJ,_EAS,_FIN,_NFE}_{AC,AN,AF},Ensembl_{geneid,transcriptid})"
tmp=${dir}/TMP_`basename ${out}`
I1=${vcf}.MakeSitesOnlyVcf
O1=${tmp}.MakeSitesOnlyVcf.gnomAD.MPC
I2=${tmp}.MakeSitesOnlyVcf.gnomAD.MPC

## STEP0: convert to sites only -------------------------------------------------------------------

if [[ ! -f ${vcf}.MakeSitesOnlyVcf.vcf.gz ]]
then
    
    java -jar ${PICARD} MakeSitesOnlyVcf INPUT=${vcf}.vcf.gz OUTPUT=${vcf}.MakeSitesOnlyVcf.vcf

    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${tmp}.MakeSitesOnlyVcf.vcf.success
    else
	touch ${tmp}.MakeSitesOnlyVcf.vcf.fail
	exit 1
    fi
    
    cat ${vcf}.MakeSitesOnlyVcf.vcf | bgzip > ${vcf}.MakeSitesOnlyVcf.vcf.gz && tabix ${vcf}.MakeSitesOnlyVcf.vcf.gz
    
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${tmp}.MakeSitesOnlyVcf.vcf.gz.success
    else
	touch ${tmp}.MakeSitesOnlyVcf.vcf.gz.fail
	exit 1
    fi

fi

##STEP1: gnomad and mpc scores -------------------------------------------------------------------

if [[ ! -f ${O1}.success ]]
then

    bcftools view -Ou -c 0 ${I1}.vcf.gz | bcftools norm -Ou -m -any | \
	bcftools annotate -Ou -x ID | \
	bcftools view -c 0 | \
	java -Xmx2g -jar ${snpeff}/SnpSift.jar dbnsfp -v -f $(echo ${dbnsfp_fields}|tr ' ' ',') -db ${dbnsfp} - | \
	grep -v "^##bcftools\|^##SnpEff\|^##SnpSift" | \
	bgzip > ${O1}.vcf.gz && tabix -f ${O1}.vcf.gz

    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${O1}.success
    else
	touch ${O1}.fail
	exit 1
    fi

fi

if [[ ! -f ${out}.BCFQUERYOUT.success ]]
then
    echo "CHROM POS REF ALT Ensembl_geneid Ensembl_transcriptid MPC_score gnomAD_exomes_AC gnomAD_exomes_AF gnomAD_exomes_AN gnomAD_exomes_AFR_AC gnomAD_exomes_AFR_AF gnomAD_exomes_AFR_AN gnomAD_exomes_AMR_AC gnomAD_exomes_AMR_AF gnomAD_exomes_AMR_AN gnomAD_exomes_ASJ_AC gnomAD_exomes_ASJ_AF gnomAD_exomes_ASJ_AN gnomAD_exomes_EAS_AC gnomAD_exomes_EAS_AF gnomAD_exomes_EAS_AN gnomAD_exomes_FIN_AC gnomAD_exomes_FIN_AF gnomAD_exomes_FIN_AN gnomAD_exomes_NFE_AC gnomAD_exomes_NFE_AF gnomAD_exomes_NFE_AN gnomAD_exomes_SAS_AC gnomAD_exomes_SAS_AF gnomAD_exomes_SAS_AN gnomAD_genomes_AC gnomAD_genomes_AF gnomAD_genomes_AN gnomAD_genomes_AFR_AC gnomAD_genomes_AFR_AF gnomAD_genomes_AFR_AN gnomAD_genomes_AMR_AC gnomAD_genomes_AMR_AF gnomAD_genomes_AMR_AN gnomAD_genomes_ASJ_AC gnomAD_genomes_ASJ_AF gnomAD_genomes_ASJ_AN gnomAD_genomes_EAS_AC gnomAD_genomes_EAS_AF gnomAD_genomes_EAS_AN gnomAD_genomes_FIN_AC gnomAD_genomes_FIN_AF gnomAD_genomes_FIN_AN gnomAD_genomes_NFE_AC gnomAD_genomes_NFE_AF gnomAD_genomes_NFE_AN" | tr ' ' '\t' > ${out}.BCFQUERYOUT
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/dbNSFP_Ensembl_geneid\t%INFO/dbNSFP_Ensembl_transcriptid\t%INFO/dbNSFP_MPC_score\t%INFO/dbNSFP_gnomAD_exomes_AC\t%INFO/dbNSFP_gnomAD_exomes_AF\t%INFO/dbNSFP_gnomAD_exomes_AN\t%INFO/dbNSFP_gnomAD_exomes_AFR_AC\t%INFO/dbNSFP_gnomAD_exomes_AFR_AF\t%INFO/dbNSFP_gnomAD_exomes_AFR_AN\t%INFO/dbNSFP_gnomAD_exomes_AMR_AC\t%INFO/dbNSFP_gnomAD_exomes_AMR_AF\t%INFO/dbNSFP_gnomAD_exomes_AMR_AN\t%INFO/dbNSFP_gnomAD_exomes_ASJ_AC\t%INFO/dbNSFP_gnomAD_exomes_ASJ_AF\t%INFO/dbNSFP_gnomAD_exomes_ASJ_AN\t%INFO/dbNSFP_gnomAD_exomes_EAS_AC\t%INFO/dbNSFP_gnomAD_exomes_EAS_AF\t%INFO/dbNSFP_gnomAD_exomes_EAS_AN\t%INFO/dbNSFP_gnomAD_exomes_FIN_AC\t%INFO/dbNSFP_gnomAD_exomes_FIN_AF\t%INFO/dbNSFP_gnomAD_exomes_FIN_AN\t%INFO/dbNSFP_gnomAD_exomes_NFE_AC\t%INFO/dbNSFP_gnomAD_exomes_NFE_AF\t%INFO/dbNSFP_gnomAD_exomes_NFE_AN\t%INFO/dbNSFP_gnomAD_exomes_SAS_AC\t%INFO/dbNSFP_gnomAD_exomes_SAS_AF\t%INFO/dbNSFP_gnomAD_exomes_SAS_AN\t%INFO/dbNSFP_gnomAD_genomes_AC\t%INFO/dbNSFP_gnomAD_genomes_AF\t%INFO/dbNSFP_gnomAD_genomes_AN\t%INFO/dbNSFP_gnomAD_genomes_AFR_AC\t%INFO/dbNSFP_gnomAD_genomes_AFR_AF\t%INFO/dbNSFP_gnomAD_genomes_AFR_AN\t%INFO/dbNSFP_gnomAD_genomes_AMR_AC\t%INFO/dbNSFP_gnomAD_genomes_AMR_AF\t%INFO/dbNSFP_gnomAD_genomes_AMR_AN\t%INFO/dbNSFP_gnomAD_genomes_ASJ_AC\t%INFO/dbNSFP_gnomAD_genomes_ASJ_AF\t%INFO/dbNSFP_gnomAD_genomes_ASJ_AN\t%INFO/dbNSFP_gnomAD_genomes_EAS_AC\t%INFO/dbNSFP_gnomAD_genomes_EAS_AF\t%INFO/dbNSFP_gnomAD_genomes_EAS_AN\t%INFO/dbNSFP_gnomAD_genomes_FIN_AC\t%INFO/dbNSFP_gnomAD_genomes_FIN_AF\t%INFO/dbNSFP_gnomAD_genomes_FIN_AN\t%INFO/dbNSFP_gnomAD_genomes_NFE_AC\t%INFO/dbNSFP_gnomAD_genomes_NFE_AF\t%INFO/dbNSFP_gnomAD_genomes_NFE_AN\n' ${O1}.vcf.gz >> ${out}.BCFQUERYOUT
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${out}.BCFQUERYOUT.success
    else
	touch ${out}.BCFQUERYOUT.fail
	exit 1
    fi
fi

## STEP2: loftee -------------------------------------------------------------------

if [[ ! -f ${out}.tab.gz.success ]]
then
    vep -i ${I2}.vcf.gz --tab -o ${out}.tab --cache --dir_cache ${cachedir} --offline --force_overwrite \
	--plugin LoF,loftee_path:${lofteedir},human_ancestor_fa:${lofteeres}/human_ancestor.fa.gz

    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${out}.tab.success
    else
	touch ${out}.tab.fail
	exit 1
    fi
    
    cat ${out}.tab | bgzip > ${out}.tab.gz && tabix ${out}.tab.gz
    
    STATUS=$?
    if [[ ${STATUS} -eq 0 ]]
    then
	touch ${out}.tab.gz.success
    else
	touch ${out}.tab.gz.fail
	exit 1
    fi
fi

##STEP4: clean up
rm ${tmp}*


