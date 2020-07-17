#! /bin/bash

#$ -N rfmix
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe smp 8


### Elizabeth G. Atkinson 
### 2/22/18
## post-processing shapeit haps/sample files to be input into RFmix for cohort local ancestry inference

## Usage is sh Merge_Phase_RFmix.sh <data-stem> <ref-stem> <genetic-recomb-map> <ancestry-ref-map>
## the script will only consider the autosomes unless modified
## data stem is the core filename before "*chr{1-22}.haps/sample"
## shapeit2, RFMix v2, plink, and vcftools are expected to be in the path 
## GEN expected to be in the format of the hapmap recombination map
## ancestry reference map assigns reference individuals to their ancestral populations of origin, as described in the RFMix manual


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


##  Unpack the parameters into labelled variables
BASE=$1
CHROM=$2

DATA=chr${CHROM}_${BASE}
cd $DATA

REF="${baseDir}/1000GP_Phase3_3WayAdmixture_reference/1000GP_Phase3_3WayAdmixture_ref"

#Genetic map
GEN=${baseDir}/1000GP_Phase3/genetic_map_chr${CHROM}_combined_b37.txt

#Sample map
awk '{gsub("CEU", "EUR"); gsub("YRI", "AFR");print $2"\t"$1}' $REF.fam > ref.sample.map
MAP=ref.sample.map

#cohort data is now formatted to merge properly with the 1000G reference panel from step 1 - suffixed .QCed
#merge 1000G and cohort data
plink --bfile $REF --bmerge $DATA.QCed --chr $CHROM --make-bed --out $DATA.QCed.1kmerge

#filter merged dataset to only include well-genotyped sites present on both the cohort and 1000G platforms - 90% genotyping rate and MAF >= 0.5%
plink --bfile $DATA.QCed.1kmerge --allow-no-sex --make-bed --geno 0.1 --maf 0.005 --out $DATA.QCed.1kmerge.filt

#separate out the chromosomes for phasing
#for i in {1..22}; do plink --bfile $DATA.QCed.1kmerge.filt --allow-no-sex --chr ${i} --make-bed --out $DATA.QCed.1kmerge.filt.chr${i} ;done

#then phase them with SHAPEIT2
#assuming all chroms are present in the same genetic map file linked to in initial command
###NOTE - this conducts joint phasing with the reference panel. In many cases you'll want to phase the cohort data using the reference panel as a separate flag.
#That can be instead implemented in SHAPEIT2 with a flag similar to:
#  --input-ref reference.haplotypes.gz reference.legend.gz reference.sample \

# Shapeit v2
shapeit --input-bed $DATA.QCed.1kmerge.filt -M $GEN -O $DATA.QCed.1kmerge.filt.phased --thread 8

##also make a list of the individuals 
cut -d' ' -f2 $DATA.QCed.fam > $DATA.indivs.txt

#convert the shapeit output into VCF format to put into RFmixv2...
shapeit -convert --input-haps $DATA.QCed.1kmerge.filt.phased --output-vcf $DATA.QCed.1kmerge.filt.phased.vcf

#make a vcf file of just the cohort individuals
vcftools --vcf $DATA.QCed.1kmerge.filt.phased.vcf --keep $DATA.indivs.txt --recode --out $DATA.QCed.1kmerge.filt.phased.cohort

#make a vcf file of just the ref individuals, assuming they're everyone who wasn't in the cohort
vcftools --vcf $DATA.QCed.1kmerge.filt.phased.vcf --remove $DATA.indivs.txt --recode --out $DATA.QCed.1kmerge.filt.phased.ref

#bgzip these
bgzip $DATA.QCed.1kmerge.filt.phased.cohort.recode.vcf
bgzip $DATA.QCed.1kmerge.filt.phased.ref.recode.vcf

#and run RFmix. Split for each chromosome separately
#the recombination map might need to be further processed to make RFMix happy depending on the format

awk -v chrom=${CHROM} '$1!~/position/{print chrom"\t"$1"\t"$3}' ${GEN} > genmap_${CHROM}

rfmix -f $DATA.QCed.1kmerge.filt.phased.cohort.recode.vcf.gz -r $DATA.QCed.1kmerge.filt.phased.ref.recode.vcf.gz --chromosome=$CHROM -m $MAP -g genmap_${CHROM} -n 5 -e 1 --reanalyze-reference --n-threads=8 -o $DATA.rfmix

## clean up
if [ -s $DATA.rfmix.msp.tsv ]
then
    mkdir -p tmp
    mv *recode* tmp/
    mv *rfmix* tmp/
    rm -f shapeit_*
    rm -f $DATA*
    rm -f genmap_$CHROM
    rm -f ref.sample.map
    mv tmp/* .
    rm -fr tmp
fi


