######################################################################################################################################################################################################################################################################################################
##################################################################################################
##################################################################################################

#!/bin/bash
##########
set -euo pipefail
export LC_ALL=C

###
# This is to use the nextflow hic pipeline
# Run with all fastqs for hic in subfolder

topDir="~/"
# genomeID="hg38+TB40E+dm6"
genomeID="hg38+TB40E"
# scripts directory, contains scripts/
scriptsDir="~/"
# Reference genome direction, contains fasta files 
refDir="~/${genomeID}"
refSeq="${refDir}/${genomeID}.fna"
account="XXX"
# default queue
queue="XXX"
# default queue time
queue_time="48:00:00"
# processors per node for alignment
# ppnAlign=24
# ppnAlign=32
# email for completion
emailaddress="${USER}@XXX.edu"

#output messages for log files/debugging
outDir="$topDir/debug"
## Create output directory, used for reporting commands output
if [ ! -d "$outDir" ]; then
        mkdir "$outDir"
        chmod 777 "$outDir"
fi

## Where to store all the files
outputdir=${topDir}"/output"
## Create output directory
if [ ! -d "$outputdir" ]; then
        mkdir "$outputdir"
        chmod 777 "$outputdir"
fi

## https://github.com/nf-core/configs/blob/master/docs/nu_genomics.md
# Used the following for nextflow
## cd bin
## curl -s https://get.nextflow.io | bash
## export PATH=~/bin:$PATH

#     name            = 'nf-core/hic'
#     author          = """Nicolas Servant"""
#     homePage        = 'https://github.com/nf-core/hic'
#     description     = """Analysis of Chromosome Conformation Capture data (Hi-C)"""
#     mainScript      = 'main.nf'
#     nextflowVersion = '!>=22.10.1'
#     version         = '2.1.0'
#     doi             = ''

## when run below: N E X T F L O W  ~  version 23.10.1

# Example
# nextflow run nf-core/cutandrun --input samplesheet.csv  --outdir <OUTDIR> --genome GRCh37 -profile docker
## https://nf-co.re/docs/usage/configuration#running-nextflow-on-your-system


jid=`sbatch <<- COMBINE-NF-HIC | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/COMBINE-NF-HIC-%j.out
#SBATCH -e $outDir/COMBINE-NF-HIC-%j.err
#SBATCH -t $queue_time
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH -J "COMBINE-NF-HIC"
sleep 1; chmod -R 777 $topDir/debug
date

module purge
module load singularity/latest
# module load nextflow/22.10.1
## already loaded nextflow version 23.10.1.5891
module load graphviz/2.40.1
module load java/jdk11.0.10

nextflow run nf-core/hic \
-r 2.1.0 \
--input $samplesheet_1 \
--outdir $outputdir \
--split_fastq TRUE --fastq_chunks_size 20000000 \
--fasta $refSeq \
--bwt2_index ${refDir}/bowtie/bowtie2 \
--digestion arima \
--chromosome_size $refDir/${genomeID}.chrom.sizes_trimmed \
--restriction_fragments $refDir/${genomeID}_Arima_restriction_fragments.bed \
--hicpro_maps TRUE \
--bin_size 50000,100000,250000,500000 \
-profile nu_genomics \
--email $emailaddress \
--skip_tads FALSE --res_tads 40000,20000 --tads_caller hicexplorer,insulation \
--skip_compartments FALSE

date
COMBINE-NF-HIC`


jid=`sbatch <<- COMBINE_2-NF-HIC | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/COMBINE_2-NF-HIC-%j.out
#SBATCH -e $outDir/COMBINE_2-NF-HIC-%j.err
#SBATCH -t $queue_time
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH -J "COMBINE_2-NF-HIC"
sleep 1; chmod -R 777 $topDir/debug
date

module purge
module load singularity/latest
# module load nextflow/22.10.1
## already loaded nextflow version 23.10.1.5891
module load graphviz/2.40.1
module load java/jdk11.0.10

nextflow run nf-core/hic \
-r 2.1.0 \
--input $samplesheet_2 \
--outdir $outputdir \
--split_fastq TRUE --fastq_chunks_size 20000000 \
--fasta $refSeq \
--bwt2_index ${refDir}/bowtie/bowtie2 \
--digestion arima \
--chromosome_size $refDir/${genomeID}.chrom.sizes_trimmed \
--restriction_fragments $refDir/${genomeID}_Arima_restriction_fragments.bed \
--hicpro_maps TRUE \
--bin_size 50000,100000,250000,500000 \
-profile nu_genomics \
--email $emailaddress \
--skip_tads FALSE --res_tads 40000,20000 --tads_caller hicexplorer,insulation \
--skip_compartments FALSE

date
COMBINE_2-NF-HIC`

jid=`sbatch <<- COMBINE_3-NF-HIC | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/COMBINE_3-NF-HIC-%j.out
#SBATCH -e $outDir/COMBINE_3-NF-HIC-%j.err
#SBATCH -t $queue_time
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH -J "COMBINE_3-NF-HIC"
sleep 1; chmod -R 777 $topDir/debug
date

module purge
module load singularity/latest
# module load nextflow/22.10.1
## already loaded nextflow version 23.10.1.5891
module load graphviz/2.40.1
module load java/jdk11.0.10

nextflow run nf-core/hic \
-r 2.1.0 \
--input $samplesheet_3 \
--outdir $outputdir \
--split_fastq TRUE --fastq_chunks_size 20000000 \
--fasta $refSeq \
--bwt2_index ${refDir}/bowtie/bowtie2 \
--digestion arima \
--chromosome_size $refDir/${genomeID}.chrom.sizes_trimmed \
--restriction_fragments $refDir/${genomeID}_Arima_restriction_fragments.bed \
--hicpro_maps TRUE \
--bin_size 50000,100000,250000,500000 \
-profile nu_genomics \
--email $emailaddress \
--skip_tads FALSE --res_tads 40000,20000 --tads_caller hicexplorer,insulation \
--skip_compartments FALSE

date
COMBINE_3-NF-HIC`

jid=`sbatch <<- KD-COMBINE-NF-HIC | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/KD-COMBINE-NF-HIC-%j.out
#SBATCH -e $outDir/KD-COMBINE-NF-HIC-%j.err
#SBATCH -t $queue_time
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH -J "KD-COMBINE-NF-HIC"
sleep 1; chmod -R 777 $topDir/debug
date

module purge
module load singularity/latest
# module load nextflow/22.10.1
module load graphviz/2.40.1
module load java/jdk11.0.10

nextflow run nf-core/hic \
-r 2.1.0 \
--input $samplesheet_KD \
--outdir $outputdir \
--split_fastq TRUE --fastq_chunks_size 20000000 \
--fasta $refSeq \
--bwt2_index ${refDir}/bowtie/bowtie2 \
--digestion arima \
--chromosome_size $refDir/${genomeID}.chrom.sizes_trimmed \
--restriction_fragments $refDir/${genomeID}_Arima_restriction_fragments.bed \
--hicpro_maps TRUE \
--bin_size 50000,100000,250000,500000 \
-profile nu_genomics \
--email $emailaddress \
--skip_tads FALSE --res_tads 40000,20000 --tads_caller hicexplorer,insulation --skip_compartments FALSE 

date
KD-COMBINE-NF-HIC`

jid=`sbatch <<- KD2-COMBINE-NF-HIC | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/KD2-COMBINE-NF-HIC-%j.out
#SBATCH -e $outDir/KD2-COMBINE-NF-HIC-%j.err
#SBATCH -t $queue_time
#SBATCH -N 1
#SBATCH --ntasks=8
#SBATCH --mem=120G
#SBATCH -J "KD2-COMBINE-NF-HIC"
sleep 1; chmod -R 777 $topDir/debug
date

module purge
module load singularity/latest
# module load nextflow/22.10.1
module load graphviz/2.40.1
module load java/jdk11.0.10

nextflow run nf-core/hic \
-r 2.1.0 \
--input $samplesheet_KD2 \
--outdir $outputdir \
--split_fastq TRUE --fastq_chunks_size 20000000 \
--fasta $refSeq \
--bwt2_index ${refDir}/bowtie/bowtie2 \
--digestion arima \
--chromosome_size $refDir/${genomeID}.chrom.sizes_trimmed \
--restriction_fragments $refDir/${genomeID}_Arima_restriction_fragments.bed \
--hicpro_maps TRUE \
--bin_size 50000,100000,250000,500000 \
-profile nu_genomics \
--email $emailaddress \
--skip_tads FALSE --res_tads 40000,20000 --tads_caller hicexplorer,insulation --skip_compartments FALSE 

date
KD2-COMBINE-NF-HIC`
