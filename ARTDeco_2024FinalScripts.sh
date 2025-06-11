##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################
######################### ARTDeco Analysis (Transcriptional readthrough) #########################
##################################################################################################

# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03551-0
# https://github.com/sjroth/ARTDeco/

# conda create --name ARTD3
# conda install -n ARTD3 python=3.6.6 bedops=2.4.40 homer=4.10 samtools=1.9 # R=4.1.2
# conda install -n ARTD3 numpy=1.16.6 pandas=0.24.2
# 
# conda install -n ARTD3 bx-python=0.8.9
# conda install -n ARTD3 rseqc=3.0.1
# 
# conda install -n ARTD3 rpy2=2.9.4
# conda install -n ARTD3 networkx=2.5.1
# conda install -n ARTD3 bioconductor-deseq2=1.20.0
# 
# conda list -n ARTD3

# This creates an environment called ARTDeco.
# Installation
# Once the prerequisites are installed, go into the same directory as setup.py and run the following code:
# wget https://github.com/sjroth/ARTDeco/blob/master/setup.py
# cd ~/software
# git clone https://github.com/sjroth/ARTDeco.git
# python setup.py install

conda activate ARTD3

### running with no MODE option
# ARTDeco -gtf-file ~/final.GRCh38.gtf \
# -chrom-sizes-file ~/hg38.chrom.sizes_mod \
# -bam-files-dir ~/TTseq_bams \
# -home-dir ~/home \
# -layout PE -stranded TRUE -orientation Reverse -cpu 8 -skip-bam-summary \
# -meta-file ~/ARTDmeta_final.txt \
# -comparisons-file ~/ARTDcomparisons_final.txt
# 
# ARTDeco -gtf-file ~/final.GRCh38.gtf \
# -chrom-sizes-file ~/hg38.chrom.sizes_mod \
# -bam-files-dir ~/TRseq_bams \
# -home-dir ~/home \
# -layout PE -stranded TRUE -cpu 8 -skip-bam-summary \
# -meta-file ~/ARTDmetaTR_final.txt \
# -comparisons-file ~/ARTDcomparisons_final.txt

################
ARTDeco -mode preprocess -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TTseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -orientation Reverse -cpu 8 -skip-bam-summary \
-meta-file ~/ARTDmeta_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt

ARTDeco -mode preprocess -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TRseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -cpu 8 -skip-bam-summary \
-meta-file ~/ARTDmetaTR_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt

###
ARTDeco -mode readthrough -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TTseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -orientation Reverse -cpu 8 -skip-bam-summary \
-meta-file ~/ARTDmeta_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt \
-summary-genes 3000

ARTDeco -mode readthrough -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TRseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -cpu 8 -skip-bam-summary \
-meta-file ~/ARTDmetaTR_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt \
-summary-genes 3000

###
ARTDeco -mode get_dogs -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TTseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -orientation Reverse -cpu 8 -skip-bam-summary \
-meta-file ~/ARTDmeta_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt

ARTDeco -mode get_dogs -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TRseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -cpu 8 -skip-bam-summary \
-meta-file ~/ARTDmetaTR_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt

###
ARTDeco -mode diff_exp_read_in -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TTseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -orientation Reverse -cpu 8 -skip-bam-summary \
-overwrite \
-meta-file ~/ARTDmeta_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt

ARTDeco -mode diff_exp_read_in -gtf-file ~/final.GRCh38.gtf \
-chrom-sizes-file ~/hg38.chrom.sizes_mod \
-bam-files-dir ~/TRseq_bams \
-home-dir ~/home \
-layout PE -stranded TRUE -cpu 8 -skip-bam-summary \
-overwrite \
-meta-file ~/ARTDmetaTR_final.txt \
-comparisons-file ~/ARTDcomparisons_final.txt


