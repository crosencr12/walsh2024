################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

#!/bin/bash
##########
#
# [topDir]/fastq  - Should contain the fastq files. This code assumes that
#                   there is an "R" in the appropriate files, i.e. *R*.fastq
#
# To run this code use the following: bash /projects/b1154/cdr5028/scripts/ttseq_pipeline_v4.sh
shopt -s extglob
version_TTSeq="CR_4.0" 
export LANG=C #In order for proper perl locale
emailaddress="${USER}@e.northwestern.edu"

# top level directory, can also be set in options
topDir=$(pwd)
# unique name for jobs in this run
samplename=$(basename $topDir)
# output data to this directory
outputdir=${topDir}"/output"
# genome ID
genomeID="hg38+TB40E"
# scripts directory, contains scripts/
scriptsDir="/projects/b1154/cdr5028/scripts"
# trimmomatic directory
trimmomaticPath="/home/cdr5028/software/anaconda3/bin/trimmomatic/share/trimmomatic-0.39-2/trimmomatic"
# fastq_screen directory (but runs on (base) environment so no need to load or direct scripts here)
#fastq_screenDir="/home/cdr5028/software/fastq_screen_v0.15.2"
# Reference genome direction, contains fasta files 
refDir="/projects/b1154/referenceGenomes"
refSeq="${refDir}/${genomeID}/${genomeID}.fna"
# default account, $account are the $queue nodes
account="b1042"
# default queue
queue="genomics"
# default queue time
queue_time="12:00:00"
# processors per node for alignment
ppnAlign=24
#output messages for log files/debugging
outDir="$topDir/debug"
# fastq files should look like filename_R1.fastq and filename_R2.fastq 
# if your fastq files look different, change this value
read1str="_R1" 
read2str="_R2" 

## Read arguments                                                     
usageHelp="Usage: ${0##*/} [list] [-g genomeID] [-h]"
genomeHelp="* [genomeID] must be defined in the script"
helpHelp="* -h: print this help and exit"

printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$genomeHelp"
    echo -e "$helpHelp"
    exit "$1"
}

while getopts "g:h" opt; do
    case $opt in
	g) genomeID=$OPTARG ;;
	h) printHelpAndExit 0;;
	[?]) printHelpAndExit 1;;
    esac
done

## Check that index for refSeq exists
if [ ! -e "${refDir}/${genomeID}/STAR/SA" ]; then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run STAR genomeGenerate on this file before running.";
    exit 100;
fi

if [ -z "$genomePath" ]
then
        #If no path to genome is give, set it here.
        genomePath=$refDir/$genomeID/$genomeID.chrom.sizes
fi

## Directories to be created and regex strings for listing files
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/output"
tmpdir=${topDir}"/tmp"

## Create output directory, used for reporting commands output
if [ ! -d "$outDir" ]; then
        mkdir "$outDir"
        chmod 777 "$outDir"
fi

## Create temporary directory, used for sort later
if [ ! -d "$tmpdir" ]; then
	mkdir "$tmpdir"
	chmod 777 "$tmpdir"
fi

## Check that fastq directory exists and has proper fastq files
if [ ! -d "$topDir/fastq" ]; then
	echo "Directory \"$topDir/fastq\" does not exist."
	echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
	exit 100
else 
    if stat -t ${fastqdir} >/dev/null 2>&1
    then
	echo "(-: Looking for fastq files...fastq files exist"
    else
	echo "***! Failed to find any files matching ${fastqdir}"
	exit 1		
    fi
fi

# Check if fastq files are trimmed and if so remove the trimmed files so they don't get re-trimmed later:
if ls $topDir/fastq/*trimmed* 1> /dev/null 2>&1; then
    echo "*** Trimmed fastq files exist, removing them..."
    trimmedFiles=$(ls $topDir/fastq/*trimmed*)
    for file in $trimmedFiles; do
    	rm $file
    done
fi

## Create output directory
if [[ -d "$outputdir" ]] 
then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    exit 100			
else
    mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 100; } 
fi

## Arguments have been checked and directories created
## NOW THE FUN BEGINS

# Get the fastq files and check to make sure paired:
fastqR1files=${topDir}"/fastq/*_R1*.fastq*"

for read1file in ${fastqR1files[@]}; do
	if [[ ! -f ${read1file/_R1/_R2} ]]; then
		echo "***! Cannot find read 2 for fastq $read1file.  Please make sure read 2 exists before proceeding."
		exit 1
	else
		echo "Read 2 file for read 1 file $read1file is ${read1file/_R1/_R2}"
		# Get the name without the read1str for future use.
		name=$(basename ${read1file%$read1str*})
		echo "   Basename for these fastqs will be $name."
	fi
done

# Add header containing command executed and timestamp:
jid=`sbatch <<- HEADER | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/${samplename}_cmd-%j.out
#SBATCH -e $outDir/${samplename}_cmd-%j.err
#SBATCH -t 00:02:00
#SBATCH -c 1
#SBATCH -J "${samplename}_cmd"
sleep 1; chmod -R 777 ${topDir}/debug/
date
echo "$0 $@"
echo "TT-seq version:$version_TTSeq" 
date
HEADER`

# First do QC.  Must have FastQ screen installed for this (FastQC comes with biobuilds).

if [ ! -d "$topDir/qc_output" ]; then
	mkdir $topDir/qc_output
fi

echo -e "(-: Running FastQC and FastQ Screen on $fastqdir in queue $queue"

for file in $fastqdir; do
sleep 0.5
fname=$(basename $file)

# Perform QC on the fastqs
jid=`sbatch <<- QC | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=4gb
#SBATCH -o $outDir/${samplename}_qc_${fname}-%j.out
#SBATCH -e $outDir/${samplename}_qc_${fname}-%j.err
#SBATCH -t 24:00:00
#SBATCH -J "${samplename}_qc_${fname}"
sleep 1; chmod -R 777 $topDir/debug

module load fastqc/0.11.5

fastqc -t 4 $file --outdir=$topDir/qc_output
fastqc -t 4 ${file/_R1/_R2} --outdir=$topDir/qc_output

module load bwa

fastq_screen --aligner bwa --force $file --outdir=$topDir/qc_output
fastq_screen --aligner bwa --force ${file/_R1/_R2} --outdir=$topDir/qc_output

QC`

done

# Trim adaptor sequences from fastq

trimR1files=()
trimR2files=()
dependtrim="afterok"
for read1file in ${fastqR1files[@]}; do
sleep 0.5
name=$(basename ${read1file%$read1str*})
trimR1files+=($topDir/fastq/${name}_trimmed_1P.fastq.gz)
trimR2files+=($topDir/fastq/${name}_trimmed_2P.fastq.gz)

trimmomaticVersion=$(trimmomatic -version)
trimParameters="ILLUMINACLIP:${trimmomaticPath}/adapters/TruSeq3-PE.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

jid=`sbatch <<- TRIM | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/${samplename}_trim_${name}-%j.out
#SBATCH -e $outDir/${samplename}_trim_${name}-%j.err
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH -c ${ppnAlign}
#SBATCH --mem-per-cpu=4g
#SBATCH -J "${samplename}_trim_${name}"
sleep 1; chmod -R 777 $topDir/debug
date


echo "Trimming fastq files using trimmomatic from $trimmomaticVersion and trimming parameters $trimParameters"

trimmomatic PE -threads ${ppnAlign} -basein $read1file -baseout $topDir/fastq/${name}_trimmed.fastq.gz $trimParameters

date

TRIM`
dependtrim="$dependtrim:$jid"
done

trimR1list=$( IFS=$','; echo "${trimR1files[*]}" )
trimR2list=$( IFS=$','; echo "${trimR2files[*]}" )

# Now align the files
# Use STAR for alignment
# Use the ENCODE options in the STAR manual 2.7.5c
# See also https://www.encodeproject.org/rna-seq/long-rnas/

jid=`sbatch <<- STAR | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/${samplename}_STAR-%j.out
#SBATCH -e $outDir/${samplename}_STAR-%j.err
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH -c ${ppnAlign}
#SBATCH -d ${dependtrim}
#SBATCH --mem-per-cpu=5g
#SBATCH -J "${samplename}_STAR"
sleep 1; chmod -R 777 $topDir/debug
date

conda activate tt-seq

echo 'Aligning paired-end reads with STAR.'

if ! STAR \
--runThreadN $ppnAlign \
--runMode alignReads \
--genomeDir $refDir/${genomeID}/STAR \
--readFilesIn $trimR1list $trimR2list \
--readFilesCommand zcat \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--quantMode TranscriptomeSAM GeneCounts \
--outReadsUnmapped Fastq \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM $(expr $ppnAlign \* 4000000000)\
--outSAMmultNmax 1 \
--outFileNamePrefix $outputdir/${samplename}_
then
	exit 100
else
	samtools index $outputdir/${samplename}_Aligned.sortedByCoord.out.bam
	# Retain only properly paired reads
	samtools view -@ $ppnAlign -b -f 2 $outputdir/${samplename}_Aligned.sortedByCoord.out.bam > $outputdir/${samplename}_Aligned_properlyPaired.sortedByCoord.out.bam
	samtools index $outputdir/${samplename}_Aligned_properlyPaired.sortedByCoord.out.bam
	echo "(-: STAR align of ${samplename}_Aligned.sortedByCoord.out.bam created successfully"
fi		

echo 'Marking duplicates and multimappers with STAR.'
# STAR de-duping only works with PE data
if ! STAR \
--runThreadN $ppnAlign \
--runMode inputAlignmentsFromBAM \
--inputBAMfile $outputdir/${samplename}_Aligned_properlyPaired.sortedByCoord.out.bam \
--bamRemoveDuplicatesType UniqueIdentical \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 96000000000 \
--outFileNamePrefix $outputdir/${samplename}_markdups_
then
	exit 100
else
	# STAR only marks duplicates, so need to use samtools to remove the dups from the bam file.
	samtools view -@ $ppnAlign -b -F 1024 $outputdir/${samplename}_markdups_Processed.out.bam > $outputdir/${samplename}_nodups_Processed.out.bam
	samtools index $outputdir/${samplename}_nodups_Processed.out.bam
	echo "(-: STAR deduping of ${samplename}_markdups_Processed.out.bam created successfully"
fi

echo "Generating signal tracks with STAR."

if ! STAR \
--runThreadN $ppnAlign \
--runMode inputAlignmentsFromBAM \
--inputBAMfile $outputdir/${samplename}_nodups_Processed.out.bam \
--outWigType bedGraph \
--outWigStrand Stranded \
--outWigNorm RPM \
--outFileNamePrefix $outputdir/${samplename}_nodups_tracks_
then
	exit 100
else
	echo "(-: STAR signal tracks of $outputdir/${samplename}_nodups_tracks_Signal.Unique.str1.out.bg	$outputdir/${samplename}_nodups_tracks_Signal.UniqueMultiple.str1.out.bg	$outputdir/${samplename}_nodups_tracks_Signal.Unique.str2.out.bg	$outputdir/${samplename}_nodups_tracks_Signal.UniqueMultiple.str2.out.bg created successfully"
fi

date
STAR`
dependSTAR="afterok:$jid"

# Make tracks

jid=`sbatch <<- TRACKS | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/${samplename}_make_tracks-%j.out
#SBATCH -e $outDir/${samplename}_make_tracks-%j.err
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -d ${dependSTAR}
#SBATCH --mem-per-cpu=8g
#SBATCH -J "${samplename}_make_tracks"
sleep 1; chmod -R 777 $topDir/debug
date

conda activate ucsc

module load coreutils/9.0
LC_COLLATE=C

echo 'Converting tracks from bedGraph to bigWig.'
echo 'Using negative values for the minus strand.'

awk -F $'\t' 'BEGIN {OFS = FS} {\\\$4 = \\\$4 * -1; print \\\$0}' $outputdir/${samplename}_nodups_tracks_Signal.Unique.str1.out.bg > $outputdir/tmp.bg
sort -S 7G -k 1,1 -k 2,2n $outputdir/tmp.bg -o $outputdir/${samplename}_nodups_tracks_Signal.Unique.str1.signFlipped.out.bg
sort -S 7G -k 1,1 -k 2,2n $outputdir/${samplename}_nodups_tracks_Signal.Unique.str2.out.bg -o $outputdir/${samplename}_nodups_tracks_Signal.Unique.str2.out.bg

bedGraphToBigWig $outputdir/${samplename}_nodups_tracks_Signal.Unique.str1.signFlipped.out.bg $refDir/${genomeID}/${genomeID}.chrom.sizes $outputdir/${samplename}_minus.bw
bedGraphToBigWig $outputdir/${samplename}_nodups_tracks_Signal.Unique.str2.out.bg $refDir/${genomeID}/${genomeID}.chrom.sizes $outputdir/${samplename}_plus.bw

if [[ ! -f $outputdir/${samplename}_minus.bw && ! -f $outputdir/${samplename}_plus.bw ]]
then
	exit 100
else
	echo "(-: STAR signal tracks of $outputdir/${samplename}_minus.bw $outputdir/${samplename}_plus.bw created successfully"
fi

date
TRACKS`
dependTracks="afterok:$jid"


# Run htseqcounts to get gene counts file
jid=`sbatch <<- HTSEQ | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH --mem-per-cpu=4g
#SBATCH -o $topDir/debug/${samplename}_htseq-count_allGenes-%j.out
#SBATCH -e $topDir/debug/${samplename}_htseq-count_allGenes-%j.err
#SBATCH -t 24:00:00
#SBATCH -d ${dependTracks}
#SBATCH -c 1
#SBATCH -J "${samplename}_htseq-count_allGenes"
sleep 1; chmod -R 777 $topDir/debug
date

conda activate tt-seq

# My gtf file is based on exons, and I want to extract the gene_name for easy output!
htseq-count -f bam -r pos -s reverse -t exon -i gene_id --additional-attr=gene_name $outputdir/${samplename}_nodups_Processed.out.bam $refDir/${genomeID}/${genomeID}.gtf > $outputdir/${samplename}_gene_counts.txt

# Remove the special counters and save them to a separate file:
sed -n -e '/^__/p' $outputdir/${samplename}_gene_counts.txt > $outputdir/${samplename}_gene_counts_HTSeqStats.txt
sed -i '/^__/d' $outputdir/${samplename}_gene_counts.txt

date
HTSEQ`
dependHTSEQ="afterok:$jid"


jid=`sbatch <<- SENDEMAIL | egrep -o -e "\b[0-9]+$"
#!/bin/bash
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $topDir/debug/${samplename}_TTseq-pipeline_complete-%j.out
#SBATCH -e $topDir/debug/${samplename}_TTseq-pipeline_complete-%j.err
#SBATCH -t 24:00:00
#SBATCH -c 1
#SBATCH -J "${samplename}_TTseq-pipeline_complete"
#SBATCH -d ${dependHTSEQ}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$emailaddress
sleep 1; chmod -R 777 $topDir/debug
date

SENDEMAIL`


echo "(-: Finished adding all jobs... Now is a good time to get that cup of coffee.."


##########################################################################################################################################
