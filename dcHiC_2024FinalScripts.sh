##############################################################################
################## Interactive nodes used, but option for script #############
##############################################################################

# #!/bin/bash
# ##########
# # set -euo pipefail
# export LC_ALL=C

topDir="~/"
# genome ID
genomeID="hg38+TB40E"
# scripts directory, contains scripts/
scriptsDir="~/"
# Reference genome direction, contains fasta files 
refDir="~/${genomeID}"
refSeq="${refDir}/${genomeID}.fna"

# default account, $account are the $queue nodes
account="XXX"
# default queue
queue="XXX"
# default queue time
queue_time="24:00:00"
# email for completion
emailaddress="${USER}@.XXX"

#output messages for log files/debugging
debugdir="$topDir/debug"
## Create output directory, used for reporting commands output
if [ ! -d "$debugdir" ]; then
        mkdir "$debugdir"
        chmod 777 "$debugdir"
fi

biasesDir="$topDir/biases"
## Create output directory, used for reporting commands output
if [ ! -d "$biasesDir" ]; then
        mkdir "$biasesDir"
        chmod 777 "$biasesDir"
fi

software="~/"
fithic="~/"
hicpro="~/"
dchic="~/"

hg38_blacklist="~/hg38+TB40E/hg38-blacklist.v2sorted.bed"
centromere="~/hg38+TB40E/extra/Modeled_regions_for_hg38MOD.txt"
# centromere="~/hg38+TB40E/hg38_CENTROMEREmod.bed"


## Make scripts for different resolutions
# resolution="500000" 
# res="500kb"
# resolution="250000"
# res="250kb"
resolution="100000"
res="100kb"
# resolution="50000"
# res="50kb"
# resolution="10000"
# res="10kb"

## --bdgfile = add bedgraph to viz file generation

##############################################################################
##############################################################################
##############################################################################

#### biases files were regenerated for use in fithic (if applicable)
# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/MOCK72rep1.${resolution}_iced.matrix \
# -b $hicpro/MOCK72rep1.${resolution}_abs.bed \
# -s $hicpro/MOCK72rep1.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > MOCK72rep1.biases.gz
# rm fithic*

# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/MOCK120rep1.${resolution}_iced.matrix \
# -b $hicpro/MOCK120rep1.${resolution}_abs.bed \
# -s $hicpro/MOCK120rep1.${resolution}_iced.matrix.biases \
# -r $resolution 
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > MOCK120rep1.biases.gz
# rm fithic*

# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/MOCK72rep2.${resolution}_iced.matrix \
# -b $hicpro/MOCK72rep2.${resolution}_abs.bed \
# -s $hicpro/MOCK72rep2.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > MOCK72rep2.biases.gz
# rm fithic*

# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/MOCK120rep2.${resolution}_iced.matrix \
# -b $hicpro/MOCK120rep2.${resolution}_abs.bed \
# -s $hicpro/MOCK120rep2.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > MOCK120rep2.biases.gz
# rm fithic*


# cd /projects/b1042/WalshLab/cdr5028/hic/fithic
# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/INF72rep1.${resolution}_iced.matrix \
# -b $hicpro/INF72rep1.${resolution}_abs.bed \
# -s $hicpro/INF72rep1.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > INF72rep1.biases.gz
# rm fithic*

# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/INF96rep3.${resolution}_iced.matrix \
# -b $hicpro/INF96rep3.${resolution}_abs.bed \
# -s $hicpro/INF96rep3.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > INF96rep3.biases.gz
# rm fithic*

# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/INF72rep2.${resolution}_iced.matrix \
# -b $hicpro/INF72rep2.${resolution}_abs.bed \
# -s $hicpro/INF72rep2.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > INF72rep2.biases.gz
# rm fithic*

# python ~//utils/HiCPro2FitHiC.py \
# -i $hicpro/INF96rep4.${resolution}_iced.matrix \
# -b $hicpro/INF96rep4.${resolution}_abs.bed \
# -s $hicpro/INF96rep4.${resolution}_iced.matrix.biases \
# -r $resolution
# zcat fithic.biases.gz | awk -F '\t' '!($3=="-1")' | gzip > INF96rep4.biases.gz
# rm fithic*

# topDir="~/"
# cp *biases* $topDir



### add blacklist information
bedtools intersect -wa -c -f 1.0 -a $hicpro/MOCK72rep1.100000_abs.bed -b $hg38_blacklist > $hicpro/MOCK72rep1.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/MOCK72rep2.100000_abs.bed -b $hg38_blacklist > $hicpro/MOCK72rep2.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/INF72rep1.100000_abs.bed -b $hg38_blacklist > $hicpro/INF72rep1.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/INF72rep2.100000_abs.bed -b $hg38_blacklist > $hicpro/INF72rep2.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/MOCK120rep1.100000_abs.bed -b $hg38_blacklist > $hicpro/MOCK120rep1.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/MOCK120rep2.100000_abs.bed -b $hg38_blacklist > $hicpro/MOCK120rep2.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/INF96rep3.100000_abs.bed -b $hg38_blacklist > $hicpro/INF96rep3.100000_absclean.bed
bedtools intersect -wa -c -f 1.0 -a $hicpro/INF96rep4.100000_abs.bed -b $hg38_blacklist > $hicpro/INF96rep4.100000_absclean.bed


### using unmodified bed/matrix files from nfcore/HiCpro output
## this did work, but would be better with blacklist info in bed file
conda activate ~/dchic

## removing some chromosomes lacking parms
## have to split arms for better compartment calling
## PC1 may just be arms in infection usually

perl $software/utility/Chromosome_ArmWise_PCA/run_dcHiC_chrArms_pca_step1_tmp.pl \
--cmd=${dchic}/input.COMBINED_HiC${res}.txt \
--cen=${centromere} \
--pfx=dcHiC_${res} \
--pexcl chr13,chr14,chr15,chr21,chr22,chrY,chrM,chrEBV,chrTB40E \
--qexcl chrY,chrM,chrEBV,chrTB40E \
--cpu 4 --mem=180 --hrs=2

bash ${dchic}/parm/parm_dcHiC_${res}_job.sh

bash ${dchic}/qarm/qarm_dcHiC_${res}_job.sh

perl $software/utility/Chromosome_ArmWise_PCA/run_dcHiC_chrArms_combine_step2.pl \
--cmd=${dchic}/input.COMBINED_HiC${res}.txt

Rscript $software/dchicf.r --file ${dchic}/input.COMBINED_HiC${res}.txt --pcatype analyze --diffdir MOCKvINF_${res}

Rscript $software/dchicf.r --file ${dchic}/input.COMBINED_HiC${res}.txt --pcatype subcomp --subnum 6 --diffdir MOCKvINF_${res}


Rscript $software/dchicf.r --file input.COMBINED_HiC${res}.txt --pcatype fithic --diffdir MOCKvINF_${res} \
--fithicpath $fithic/fithic.py --pythonpath $pythonv

####
## copy .txt files needed for fithic loop calling (from both arms?)
# for file in $list;
# do rm ~//${file}_100kb_pca/intra_pca/${file}_100kb_mat/chr*.txt
# done

list="MOCK_R1 MOCK_R2 MOCK_R3 MOCK_R4 INF_R1 INF_R2 INF_R3 INF_R4"
chromlistQ="chr13 chr14 chr15 chr21 chr22"
chromlist="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr16 chr17 chr18 chr19 chr20 chrX"

for file in $list;
do cp ~//parm/${file}_100kb_pca/intra_pca/${file}_100kb_mat/chr*.txt \
~//${file}_100kb_pca/intra_pca/${file}_100kb_mat
done

for file in $list;
do
for chrom in $chromlistQ;
do cp ~//qarm/${file}_100kb_pca/intra_pca/${file}_100kb_mat/${chrom}.txt \
~//${file}_100kb_pca/intra_pca/${file}_100kb_mat
done
done

## don't combine q arm with missing parm chromosomes (chr13,chr14,chr15,chr21,chr22,chrY,chrM,chrEBV,chrTB40E)
for file in $list;
do
for chrom in $chromlist;
do tail -n+2 ~//qarm/${file}_100kb_pca/intra_pca/${file}_100kb_mat/${chrom}.txt >> \
~//${file}_100kb_pca/intra_pca/${file}_100kb_mat/${chrom}.txt
done
done
####

Rscript $software/dchicf.r --file input.COMBINED_HiC${res}.txt --pcatype dloop --diffdir MOCKvINF_${res}

Rscript $software/dchicf.r --file ${dchic}/input.COMBINED_HiC${res}.txt --pcatype enrich --region both --exclA FALSE --pcscore TRUE --diffdir MOCKvINF_${res} --genome hg38

Rscript $software/dchicf.r --file ${dchic}/input.COMBINED_HiC${res}.txt --pcatype viz --diffdir MOCKvINF_${res} --genome hg38




##############################################################################
##############################################################################
##############################################################################
##############################################################################
