#! /bin/bash

#$ -N ancestQC
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

## A script for processing cohort data files for use in LAI with 1000 genomes reference individuals
## Written by Elizabeth Atkinson. 1/23/18
# sub-scripts for processing steps courtesy of Meng Lin and Chris Gignoux

## Required parameters:
## 1. Binary plink files for the cohort data in question. Will just have to put in the bed stem (not including file extension)
## 2. plink installed and on the path
## 3. R installed and on the path
## 4. python installed and on the path
## Result: A new set of binary plink files where all allele rsIDs are renamed to dbsnp144, sites are oriented to 1000 genomes, with non-matching sites, indels, duplicates, sex chromosomes, and triallelic sites removed. Output plink files will be the input DATA name suffixed with QC.

## Usage is CohortDataQC.sh <data-stem> <dbSNP-bedfile> <ref-panel-legend, i.e. from 1kG>


if [ $HOSTNAME == 'sofia.udd.cl' ] || [[ $HOSTNAME == compute*-1-*.local ]]
then
    baseDir="/hpcudd/ICIM/boris/projects/local_ancestry/00_newpipeline"
    export PATH="/hpcudd/home/boris/miniconda3/bin:$PATH"
elif [ $HOSTNAME == 'mendel' ]
then
    baseDir="/home/boris/storage/00_papers/local_ancestry_pipeline/00_newpipeline"
else
    echo "Unrecognized host $HOSTNAME"
    echo "can't locate genome references"
    exit 1
fi


## Mendel directories
#baseDir="/home/boris/storage/00_papers/local_ancestry_pipeline/00_newpipeline"

Legend="1000GP_Phase3/autosome_1000GP_Phase3_GRCh37_combined.legend.gz"
snp151="dbSNPv151_20180423/dbSNPv151.bed.gz"
scripts="${baseDir}/Post-QC"

##  Unpack the parameters into labelled variables
DATA=$1
CHROM=$2

## Make working directory
prefix=chr${CHROM}_${DATA}
mkdir -p $prefix

## Get ref snps and legend

tabix ${baseDir}/${snp151} ${CHROM} > ${prefix}/chr${CHROM}.db151.bed
DBSNP=chr${CHROM}.db151.bed

tabix ${baseDir}/${Legend} ${CHROM} > ${prefix}/chr${CHROM}.1000GP.legend
LEG=chr${CHROM}.1000GP.legend

## keep only the autosomes in the data file
plink --bfile $DATA --chr $CHROM --make-bed --out ${prefix}/chr${CHROM}_$DATA.auto

DATA=${prefix}
cd ${prefix}

## Fix "." in bim
mv $DATA.auto.bim $DATA.auto.bim.old
awk 'BEGIN{OFS="\t"}{if($2~/rs/) print $0; else print $1, $1":"$4,$3,$4,$5,$6}' $DATA.auto.bim.old > $DATA.auto.bim


## Find and get rid of duplicate loci in the bim file
#then keep the good SNPs in the plink file

cut -f2,4 $DATA.auto.bim| uniq -f1 > $DATA.NonDupSNPs
cut -f2,4 $DATA.auto.bim| uniq -D -f1 > $DATA.DuplicateSNPs
cat $DATA.DuplicateSNPs | uniq -f1 > $DATA.FirstDup
cat $DATA.NonDupSNPs $DATA.FirstDup > $DATA.SNPstoKeep

#set up environment to run plink again
plink --bfile $DATA.auto --extract $DATA.SNPstoKeep --make-bed --out $DATA.auto.nodup


## Update SNP IDs to dbsnp 144
python ${scripts}/update_rsID_bim_arg.py --bim $DATA.auto.nodup.bim --bed $DBSNP --format T --codechr F --out $DATA.auto.nodup.dbsnp.bim

#copy the other files over to this name
cp $DATA.auto.nodup.bed $DATA.auto.nodup.dbsnp.bed
cp $DATA.auto.nodup.fam $DATA.auto.nodup.dbsnp.fam

#orient to 1000G
python ${scripts}/match_against_1000g_v2.py --bim $DATA.auto.nodup.dbsnp.bim --legend $LEG --out $DATA.1kg
# This script has three outputs: ##modified to be suffixed to allow for full paths to be input
# 1) [outfile].Indel.txt: a bim file of indels 
# 2) [outfile].NonMatching.txt: a bim file containing loci not found in 1000 genome, or has different coding alleles than 1000 genome (tri-allelic, for example). They should be removed.
# 3) [outfile].FlipStrand.txt: a bim file containing loci to flip. 

#combine the lists of indels and triallelic/non-matching sites into one list of bad SNPs to remove
cat $DATA.1kg.Indel.txt $DATA.1kg.NonMatching.txt > $DATA.1kg.badsites.txt

## Flip strands for flipped sites and remove non-matching loci using plink. 
plink --bfile $DATA.auto.nodup.dbsnp --exclude $DATA.1kg.badsites.txt --make-bed --out $DATA.auto.nodup.dbsnp.1ksites
plink --bfile $DATA.auto.nodup.dbsnp.1ksites --flip $DATA.1kg.FlipStrand.txt --make-bed --out $DATA.auto.nodup.dbsnp.1ksites.flip

##export a warning flag if there are too many mismatched sites compared to 1000 genomes

wc -l $DATA.1kg.NonMatching.txt > nonCount
wc -l $DATA.1kg.FlipStrand.txt > flipcount
wc -l $DATA.auto.nodup.dbsnp.bim > totalsites
paste nonCount flipcount totalsites > SiteCounts

awk '{if ($1/$5 > 0.01) print "WARNING: "$1/$5*100"% of sites are problematic when compared to 1000G. This could be indicative of a different reference build or other date file incompatibility." }' SiteCounts > Warnings.out

## Find and remove A/T C/G loci
python ${scripts}/find_cg_at_snps.py $DATA.auto.nodup.dbsnp.1ksites.flip.bim > $DATA.ATCGsites

plink --bfile $DATA.auto.nodup.dbsnp.1ksites.flip --exclude $DATA.ATCGsites --make-bed --out $DATA.QCed

#cohort data is now formatted to merge properly with the 1000G reference panel

## clean up
if [ -s $DATA.QCed.bed ]
then
    mkdir -p tmp
    mv $DATA.QCed* tmp
    mv Warnings.out tmp
    rm -f chr${CHROM}*
    rm -f *ount*
    rm -f totalsites
    mv tmp/* .
    rm -fr tmp
fi


