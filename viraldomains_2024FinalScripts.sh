##################################################################################################
##################################################################################################
############################ Viral Hi-C Contacts #################################################
##################################################################################################


# ##### Need hg38+TB40E digested genome by HiC-Pro/digest_genome.py (Arima: restriction_site='^GATC,G^ANTC, ligation_site='GATCGATC,GATCANTC,GANTGATC,GANTANTC)
# # Generate frag site for genome (both hg38 and TB40E)
# python /home/cdr5028/software/HiC-Pro/bin/utils/digest_genome.py -r G^ANTC ^GATC  -o hg38+TB40E_arima.bed $refDir/hg38+TB40E.fna
# # crop the chrUn and alt chromosomes out (saves time)
# grep -v chrUN hg38+TB40E_arima.bed > hg38+TB40E_arima-tmp.bed
# grep -v random hg38+TB40E_arima-tmp.bed > hg38+TB40E_arima-clean.bed
# ##################################################################################################
# #### Generate a viewpoints from capture site
# # Make capture site regions for TB40E as a bed file
# # TB40E = 236482bp
# # This generates TB40E capture sites from fragmented genome (by Arima enzymes)
# grep chrTB40E hg38+TB40E_arima.bed > TB40E.bed
# ### To combine all TB40E regions into one bed file (lrg = large):
# # head -1 TB40E.bed > TB40E_tmp.bed
# # awk -v FS='\t' -v OFS='\t' '{print $1, $2, 236482, $4, $5, $6}'  TB40E_tmp.bed > TB40E_lrg.bed
### Make larger windows w/ bedtools makewindows 
# bedtools makewindows -g TB40E.chrom.sizes -w 1000 -i winnum -s 500 > TB40E_1kbins.500bw.bed
##################################################################################################
################### Infected 
dataDir="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/hicpro_allValidPairs"
cores="2"
software="/home/cdr5028/software"
refDir="/projects/b1154/referenceGenomes/hg38+TB40E"
workingFiles="/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT"

### 8/28/24
## making new combined (all infected samples) TB40E capture sites
# pwd /projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/viraldomains
python $software/HiC-Pro/bin/utils/make_viewpoints.py -i $dataDir/INFc2.allValidPairs \
-f $refDir/hg38+TB40E_arima.bed  -t $refDir/TB40E.bed \
-v -o INFc2-TB40E-capture.bedgraph
cat INFc2-TB40E-capture* > INFc2-TB40E_combine.bedgraph
sort -k1,1 -k2,2n INFc2-TB40E_combine.bedgraph > INFc2-TB40E-capture-sorted.bedgraph

bedtools merge -i INFc2-TB40E-capture-sorted.bedgraph -c 4 -o sum,count > INFc2-TB40E-capture-sorted-nomerg.bedgraph
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4}' INFc2-TB40E-capture-sorted-nomerg.bedgraph > INFc2-TB40E-capture-sorted-nomerg1.bedgraph
sort -k1,1 -k2,2n INFc2-TB40E-capture-sorted-nomerg1.bedgraph > INFc2-TB40E-capture-sorted-nomerg_sorted.bedgraph
grep -v random INFc2-TB40E-capture-sorted-nomerg_sorted.bedgraph > INFc2-TB40E-capture-sorted-nomerg_sorted-clean1.bedgraph
grep -v chrUn INFc2-TB40E-capture-sorted-nomerg_sorted-clean1.bedgraph > INFc2-TB40E-capture-sorted-nomerg_sorted-cleanfinal.bedgraph
# Finally, use the UCSC bedGraphToBigWig tool:
bedGraphToBigWig INFc2-TB40E-capture-sorted-nomerg_sorted-cleanfinal.bedgraph $refDir/hg38+TB40E.chrom.sizes_trimmed INFc2-TB40E-capture-sorted-nomerg.bw

bedtools merge -i INFc2-TB40E-capture-sorted.bedgraph -c 4 -o sum,count -d 100 > INFc2-TB40E-capture-sorted-100bpmerg.bedgraph
# ## need bw file
awk '{printf "%s\t%d\t%d\t%2.3f\n" , $1,$2,$3,$4}' INFc2-TB40E-capture-sorted-100bpmerg.bedgraph  > INFc2-TB40E-capture-sorted-100bpmerg1.bedgraph
sort -k1,1 -k2,2n INFc2-TB40E-capture-sorted-100bpmerg1.bedgraph > INFc2-TB40E-capture-sorted-100bpmerg_sorted.bedgraph
grep -v random INFc2-TB40E-capture-sorted-100bpmerg_sorted.bedgraph > INFc2-TB40E-capture-sorted-100bpmerg_sorted-clean1.bedgraph
grep -v chrUn INFc2-TB40E-capture-sorted-100bpmerg_sorted-clean1.bedgraph > INFc2-TB40E-capture-sorted-100bpmerg_sorted-cleanfinal.bedgraph
# Finally, use the UCSC bedGraphToBigWig tool:
bedGraphToBigWig INFc2-TB40E-capture-sorted-100bpmerg_sorted-cleanfinal.bedgraph $refDir/hg38+TB40E.chrom.sizes_trimmed INFc2-TB40E-capture-sorted-100bpmerg.bw


######################################################################################################################################################################################################################################################################################################
##################################################################################################
##################################################################################################
######################### Viral Domains ####################
## Viral Enriched Regions (VERs) were identified by taking the file "INFc2-TB40E-capture-sorted-nomerg_sorted-cleanfinal.bedgraph"
## and mapping reads to 100kb bins of the human genome, then taking the top 10% of bins with the most contacts
