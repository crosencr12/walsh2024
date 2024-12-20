######################################################################################################################################################################################################################################################################################################
##################################################################################################
##################################################################################################



#!/bin/bash
##########
#
# This is to use the nextflow cut&run pipeline
# Run with all fastqs for cut&run are in subfolder: /projects/b1042/WalshLab/cdr5028/cutrun/fastqs
## cd /projects/b1042/WalshLab/cdr5028/cutrun/fastqs
# Then run this code use the following: bash /projects/b1154/cdr5028/scripts/cutrun_nf-pipeline_v2.sh

topDir="/projects/b1042/WalshLab/cdr5028/cutrun"
# genome ID
genomeID="hg38+TB40E"
# scripts directory, contains scripts/
scriptsDir="/projects/b1154/cdr5028/scripts"
# Reference genome direction, contains fasta files 
refDir="/projects/b1154/referenceGenomes/${genomeID}"
refSeq="${refDir}/${genomeID}.fna"
# Sample sheet located: /projects/b1042/WalshLab/cdr5028/cutrun
samplesheet="${topDir}/fastqs/nfcore_cutrun-samplesheet_v2.csv"

# default account, $account are the $queue nodes
account="b1042"
# default queue
queue="genomics"
# default queue time
queue_time="48:00:00"
# processors per node for alignment
ppnAlign=24
# email for completion
emailaddress="${USER}@e.northwestern.edu"

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

# Example
# nextflow run nf-core/cutandrun --input samplesheet.csv  --outdir <OUTDIR> --genome GRCh37 -profile docker

# If your data includes IgG controls, these will additionally be de-duplicated. 
# It is not the normal protocol to de-duplicate the target reads, however, if this is required, use the --dedup_target_reads true switch.

jid=`sbatch <<- NF-CUTRUN | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/NF-CUTRUN-%j.out
#SBATCH -e $outDir/NF-CUTRUN-%j.err
#SBATCH -t $queue_time
#SBATCH -n 1
#SBATCH -c ${ppnAlign}
#SBATCH --mem-per-cpu=5g
#SBATCH -J "NF-CUTRUN"
sleep 1; chmod -R 777 $topDir/debug
date

nextflow run nf-core/cutandrun --input $samplesheet --outdir $outputdir --fasta $refSeq \
--bowtie2 ${refDir}/bowtie/ --gtf ${refDir}/$genomeID.gtf --spikein_genome K12-MG1655 -profile nu_genomics \
--save_merged_fastq true \
--normalisation_mode Spikein --use_control true --igg_scale_factor 0.5 \
--peakcaller SEACR,MACS2 --email $emailaddress  --seacr_stringent 'relaxed' \
--macs2_narrow_peak true

date
NF-CUTRUN`
#dependSTAR="afterok:$jid"
