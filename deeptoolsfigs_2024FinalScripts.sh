######################################################################################################################################################################################################################################################################################################
##################################################################################################
##################################################################################################

srun --pty -N 1 -n 4 --mem=180G --account=b1042 --time=12:00:00 --partition=genomics bash -l


workingFiles="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT"
refDir="/projects/b1154/referenceGenomes/hg38+TB40E"
cores="2"
TT="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/transcription/cov_normbw"
genedata_v1="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/transcription/beds"

module load deeptools


##################################################################################################
##################################################################################################
### CUT&RUN
# Need to combine H3K9me3 R1 and R2 (two reps)
## files are in 2024finalfigs
# conda activate test (ucsc conda does not work for averaging)
bigwigAverage -b MOCK_H3K9me3_R1.bigWig MOCK_H3K9me3_R2.bigWig \
-o MOCK_H3K9me3_avg.bw

bigwigAverage -b INF_H3K9me3_R1.bigWig INF_H3K9me3_R2.bigWig \
-o INF_H3K9me3_avg.bw

bigwigAverage -b MOCK_H3K4me3_R1.bigWig MOCK_H3K4me3_R2.bigWig \
-o MOCK_H3K4me3_avg.bw

bigwigAverage -b INF_H3K4me3_R1.bigWig INF_H3K4me3_R2.bigWig \
-o INF_H3K4me3_avg.bw


######################################################################################################################################################################################################################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
################################## TRANSCRIPTION COVERAGE TRACKS ##################################
##################################################################################################

########### looking at transcription (coverage) between mock/inf
## module load deeptools

######## 96hr TTseq
bamCoverage -b TT09_nodups_Processed.out.bam -o TT09_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TT10_nodups_Processed.out.bam -o TT10_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TT11_nodups_Processed.out.bam -o TT11_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TT12_nodups_Processed.out.bam -o TT12_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC

######## 72hr TTseq
bamCoverage -b TT17_nodups_Processed.out.bam -o TT17_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TT18_nodups_Processed.out.bam -o TT18_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TT19_nodups_Processed.out.bam -o TT19_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TT20_nodups_Processed.out.bam -o TT20_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC

######## 96hr TRseq
bamCoverage -b TR09_nodups_Processed.out.bam -o TR09_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TR10_nodups_Processed.out.bam -o TR10_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TR11_nodups_Processed.out.bam -o TR11_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TR12_nodups_Processed.out.bam -o TR12_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC

######## 72hr TRseq
bamCoverage -b TR17_nodups_Processed.out.bam -o TR17_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TR18_nodups_Processed.out.bam -o TR18_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TR19_nodups_Processed.out.bam -o TR19_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC
bamCoverage -b TR20_nodups_Processed.out.bam -o TR20_cov_norm.bw \
--blackListFileName $refDir/hg38-blacklist.v2sorted.bed --effectiveGenomeSize 3100159023 --normalizeUsing RPGC


### Transcriptional coverage average bw
# Need to combine (two reps)
# conda activate /home/cdr5028/software/anaconda3/envs/test (ucsc conda does not work for averaging)

bigwigAverage -b  TT09_cov_norm.bw  TT10_cov_norm.bw \
-o MOCK_TT96_avg.bw
bigwigAverage -b  TT11_cov_norm.bw  TT12_cov_norm.bw \
-o INF_TT96_avg.bw
bigwigAverage -b  TT17_cov_norm.bw  TT18_cov_norm.bw \
-o MOCK_TT72_avg.bw
bigwigAverage -b  TT19_cov_norm.bw  TT20_cov_norm.bw \
-o INF_TT72_avg.bw

bigwigAverage -b  TR09_cov_norm.bw  TR10_cov_norm.bw \
-o MOCK_TR96_avg.bw
bigwigAverage -b  TR11_cov_norm.bw  TR12_cov_norm.bw \
-o INF_TR96_avg.bw
bigwigAverage -b  TR17_cov_norm.bw  TR18_cov_norm.bw \
-o MOCK_TR72_avg.bw
bigwigAverage -b  TR19_cov_norm.bw  TR20_cov_norm.bw \
-o INF_TR72_avg.bw


######################################################################################################################################################################################################################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
############################ Deeptools Figures ###################################################
##################################################################################################

workingFiles="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT"
refDir="/projects/b1154/referenceGenomes/hg38+TB40E"
cores="2"
TT="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/transcription/cov_normbw"
genedata_v1="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/transcription/beds"

####################################################
############ A / B compartments w/ virus contacts ############
####################################################
# # ## A compartment
computeMatrix scale-regions \
-S $workingFiles/hic/viraldomains/INF-TB40E_combined.bw \
-R $workingFiles/hic/compartments/INFposPCm.bed \
--regionBodyLength 100000 \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--binSize 10 \
--blackListFileName $refDir/hg38-blacklist.v2.bed \
--numberOfProcessors $cores \
--startLabel start --endLabel end \
-o matrix.100kb-TB40Ecombined-INFcompA.gz
# # 
# # 
# # ## B compartment
# # 
computeMatrix scale-regions \
-S $workingFiles/hic/viraldomains/INF-TB40E_combined.bw \
-R $workingFiles/hic/compartments/INFnegPCm.bed \
--regionBodyLength 100000 \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--binSize 10 \
--blackListFileName $refDir/hg38-blacklist.v2.bed \
--numberOfProcessors $cores \
--startLabel start --endLabel end \
-o matrix.100kb-TB40Ecombined-INFcompB.gz
# 
# # ###### Plotting
# # 
plotProfile -m matrix.100kb-TB40Ecombined-INFcompA.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best --colors purple \
-o TB40Ecombined-100kb-INFcompA_meanline.pdf

plotProfile -m matrix.100kb-TB40Ecombined-INFcompB.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best --colors purple \
-o TB40Ecombined-100kb-INFcompB_meanline.pdf

###### ###### ###### ###### ###### 
###### Subcompartments
###### ###### ###### ###### ###### 

computeMatrix scale-regions \
-S $workingFiles/hic/viraldomains/viral_COMBINED/INFc-TB40E-capture-sorted-100bpmerg.bw \
-R $workingFiles/hic/compartments/Inf_subcomp_A3m.bed $workingFiles/hic/compartments/Inf_subcomp_A2m.bed $workingFiles/hic/compartments/Inf_subcomp_A1m.bed $workingFiles/hic/compartments/Inf_subcomp_A0m.bed \
$workingFiles/hic/compartments/Inf_subcomp_B3m.bed $workingFiles/hic/compartments/Inf_subcomp_B2m.bed $workingFiles/hic/compartments/Inf_subcomp_B1m.bed $workingFiles/hic/compartments/Inf_subcomp_B0m.bed \
--regionBodyLength 100000 \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--binSize 10 \
--blackListFileName $refDir/hg38-blacklist.v2.bed \
--startLabel start --endLabel end \
-o matrix.100kb-TB40Ecombined-INFSUBcomp.gz 

plotProfile -m matrix.100kb-TB40Ecombined-INFSUBcomp.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best \
-o TB40Ecombined-100kb-INFSUBcomp_meanline.pdf \
--colors "#e31a1c" "#fb9a99" "#ff7f00" "#fdbf6f" "#1f78b4" "#a6cee3" "#33a02c" "#b2df8a"

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###
### subcompartment colors from Flourish
## red = #e31a1c, pink = #fb9a99, orange = #ff7f00, lightorange = #fdbf6f, 
## lightgreen = #b2df8a, green = #33a02c, lightblue = #a6cee3, blue = 1f78b4
###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


####################################################
############ A / B compartments w/ CUT&RUN ############
####################################################


# ####################################
# ####################################
# ###### INF
# 
# ## A compartment
# 
computeMatrix scale-regions \
-S $workingFiles/cutrun/INF_H3K9me3_avg.bw \
-R $workingFiles/hic/compartments/INFposPCm.bed \
--regionBodyLength 100000 \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--binSize 10 \
--blackListFileName $refDir/hg38-blacklist.v2.bed \
--numberOfProcessors $cores \
--startLabel start --endLabel end \
-o matrix.100kb-INF_H3K9me3-INFcompA.gz
# # 
# # ## B compartment
# # 
computeMatrix scale-regions \
-S $workingFiles/cutrun/INF_H3K9me3_avg.bw \
-R $workingFiles/hic/compartments/INFnegPCm.bed \
--regionBodyLength 100000 \
--beforeRegionStartLength 5000 \
--afterRegionStartLength 5000 \
--binSize 10 \
--blackListFileName $refDir/hg38-blacklist.v2.bed \
--numberOfProcessors $cores \
--startLabel start --endLabel end \
-o matrix.100kb-INF_H3K9me3-INFcompB.gz
# # 
# # ###### Plotting
# # 
plotProfile -m matrix.100kb-INF_H3K9me3-INFcompA.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best --colors darkgrey \
-o INF_H3K9me3-INFcompA_meanline.pdf


plotProfile -m matrix.100kb-INF_H3K9me3-INFcompB.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best --colors darkgrey \
-o INF_H3K9me3-INFcompB_meanline.pdf
# 
# 


####################################################
############ Gene sets w/ CUT&RUN ############
####################################################

#### For H3K4me3 (96hpi) over gene sets (72hpi and 96hpi)
# module load deeptools # deeptools/3.5.1
# 
# 
# 
# ############################
# ####### Look at upregulated genes (top / read-in gene lists)
# 
# 
computeMatrix reference-point --referencePoint TSS -S  $workingFiles/cutrun/MOCK_H3K4me3_avg.bw  $workingFiles/cutrun/INF_H3K4me3_avg.bw \
                              -R  $genedata_v1/All.genes_96.bed $genedata_v1/commonUPgenes.bed $genedata_v1/commonDOWNgenes.bed $genedata_v1/ARTD_RIgenesALL.bed \
                              --beforeRegionStartLength 2000 \
                              --afterRegionStartLength 2000 \
                              --binSize 10 \
                              --blackListFileName $refDir/hg38-blacklist.v2.bed \
                              -o matrix.H3K4me3-genes_refpointTSS.gz -p $cores 
#                               --skipZeros 

plotProfile -m matrix.H3K4me3-genes_refpointTSS.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best \
-o H3K4me3-genes_meanline.pdf \
--colors "#000000" "#33a02c" "#ff7f00" "#b2df8a"
     



####################################################
############ Gene sets w/ Transcriptional Coverage ############
####################################################
#### 96
## MOCK
# $genedata_v1/ARTD_RIgenes.bed \
computeMatrix scale-regions -S $TT/MOCK_TT96_avg.bw $TT/INF_TT96_avg.bw \
                              -R $genedata_v1/All.genes_96.bed $genedata_v1/commonUPgenes.bed $genedata_v1/commonDOWNgenes.bed \
                              --regionBodyLength 3000 \
                              --beforeRegionStartLength 1000 \
                              --afterRegionStartLength 1000 \
                              --binSize 10 \
                              --blackListFileName $refDir/hg38-blacklist.v2.bed \
                              -o matrix.TT96_avg-allgenesets96.gz
#                               --skipZeros \
#                               --unscaled5prime 100 --unscaled3prime 500

plotProfile -m matrix.TT96_avg-allgenesets96.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best \
-o TT96_avg-allgenesets96_meanline.pdf \
--colors "#000000" "#33a02c" "#ff7f00" 


computeMatrix scale-regions -S $TT/MOCK_TR96_avg.bw $TT/INF_TR96_avg.bw \
                              -R $genedata_v1/All.genes_96.bed $genedata_v1/commonUPgenes.bed $genedata_v1/commonDOWNgenes.bed \
                              --regionBodyLength 3000 \
                              --beforeRegionStartLength 1000 \
                              --afterRegionStartLength 1000 \
                              --binSize 10 \
                              --blackListFileName $refDir/hg38-blacklist.v2.bed \
                              -o matrix.TR96_avg-allgenesets96.gz 
#                               --skipZeros \
#                               --unscaled5prime 500 --unscaled3prime 500

plotProfile -m matrix.TR96_avg-allgenesets96.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best \
-o TR96_avg-allgenesets96_meanline.pdf \
--colors "#000000" "#33a02c" "#ff7f00" 



# $TT/TT17_cov_norm.bw  $TT/TT18_cov_norm.bw $TT/TT19_cov_norm.bw  $TT/TT20_cov_norm.bw
computeMatrix scale-regions -S $TT/MOCK_TT72_avg.bw $TT/INF_TT72_avg.bw \
                              -R $genedata_v1/All.genes_96.bed $genedata_v1/commonUPgenes.bed $genedata_v1/commonDOWNgenes.bed \
                              --regionBodyLength 3000 \
                              --beforeRegionStartLength 1000 \
                              --afterRegionStartLength 1000 \
                              --binSize 10 \
                              --blackListFileName $refDir/hg38-blacklist.v2.bed \
                              -o matrix.TT72_avg-allgenesets72.gz 
#                               --skipZeros \
#                               --unscaled5prime 500 --unscaled3prime 500

plotProfile -m matrix.TT72_avg-allgenesets72.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best \
-o TT72_avg-allgenesets72_meanline.pdf \
--colors "#000000" "#33a02c" "#ff7f00" 


computeMatrix scale-regions -S $TT/MOCK_TR72_avg.bw $TT/INF_TR72_avg.bw \
                              -R $genedata_v1/All.genes_96.bed $genedata_v1/commonUPgenes.bed $genedata_v1/commonDOWNgenes.bed \
                              --regionBodyLength 3000 \
                              --beforeRegionStartLength 1000 \
                              --afterRegionStartLength 1000 \
                              --binSize 10 \
                              --blackListFileName $refDir/hg38-blacklist.v2.bed \
                              -o matrix.TR72_avg-allgenesets72.gz 
#                               --skipZeros \
#                               --unscaled5prime 500 --unscaled3prime 500

plotProfile -m matrix.TR72_avg-allgenesets72.gz \
--averageType mean \
--plotType lines \
--startLabel start --endLabel end \
--legendLocation best \
-o TR72_avg-allgenesets72_meanline.pdf \
--colors "#000000" "#33a02c" "#ff7f00" 

