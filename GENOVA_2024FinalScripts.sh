################################################################################################
################################################################################################

#### For some analysis I needed combined cool files
# ## going to use cooler tools to just merge cool files
# ###### 
# ### Make sure you run "cooler balance" before downstream analysis

# conda activate /home/cdr5028/envs/pentads-env

# cooler merge INF72merg.100000_balanced.cool INF_72.100000_balanced.cool INF_72-2.100000_balanced.cool
# cooler merge MOCK72merg.100000_balanced.cool MOCK_72.100000_balanced.cool MOCK_72-2.100000_balanced.cool
# cooler merge INF120merg.100000_balanced.cool INF_120.100000_balanced.cool INF_120-2.100000_balanced.cool
# cooler merge MOCK120merg.100000_balanced.cool MOCK_120.100000_balanced.cool MOCK_120-2.100000_balanced.cool
# cooler merge INF96merg.100000_balanced.cool INF96rep3.100000_balanced.cool INF96rep4.100000_balanced.cool
# cooler balance INF72merg.100000_balanced.cool
# cooler balance MOCK72merg.100000_balanced.cool
# cooler balance INF120merg.100000_balanced.cool
# cooler balance MOCK120merg.100000_balanced.cool
# cooler balance INF96merg.100000_balanced.cool


# cooler merge INFmerg.100000_balanced.cool INF72merg.100000_balanced.cool INF120merg.100000_balanced.cool
# cooler merge MOCKmerg.100000_balanced.cool MOCK72merg.100000_balanced.cool MOCK120merg.100000_balanced.cool
# cooler balance INFmerg.100000_balanced.cool
# cooler balance MOCKmerg.100000_balanced.cool

# cooler merge INFmerg2.100000_balanced.cool INFmerg.100000_balanced.cool INF96rep3.100000_balanced.cool
# cooler merge INFmergALL.100000_balanced.cool INFmerg2.100000_balanced.cool INF96rep4.100000_balanced.cool
# cooler balance INFmergALL.100000_balanced.cool

# cooler merge siCINF72merg1.100000_balanced.cool siCINF72.100000_balanced.cool siCINF72rep3.100000_balanced.cool
# cooler merge siCINF72merg2.100000_balanced.cool siCINF72merg1.100000_balanced.cool siCINF72_2.100000_balanced.cool
# cooler balance siCINF72merg1.100000_balanced.cool
# cooler balance siCINF72merg2.100000_balanced.cool
# 
# cooler merge siSINF72merg1.100000_balanced.cool siSINF72.100000_balanced.cool siSINF72rep3.100000_balanced.cool
# cooler merge siSINF72merg2.100000_balanced.cool siSINF72merg1.100000_balanced.cool siSINF72_2.100000_balanced.cool
# cooler balance siSINF72merg1.100000_balanced.cool
# cooler balance siSINF72merg2.100000_balanced.cool
# 
# cooler merge siAINF72merg1.100000_balanced.cool siAINF72.100000_balanced.cool siAINF72rep3.100000_balanced.cool
# cooler merge siAINF72merg2.100000_balanced.cool siAINF72merg1.100000_balanced.cool siAINF72_2.100000_balanced.cool
# cooler balance siAINF72merg1.100000_balanced.cool
# cooler balance siAINF72merg2.100000_balanced.cool


################################################################################################
## GENOVA
# written in R
# can load cooler files and HiC-pro files AND juicer HiC files
# can do saddle (on custom?)

# install.packages("remotes")
# remotes::install_github("robinweide/GENOVA")
library(GENOVA)

# install.packages("devtools")
# devtools::install_github("robinweide/GENOVA", ref = 'dev')

# install.packages("BiocManager")
# BiocManager::install("rtracklayer")
# library(rtracklayer)
# 
# remotes::install_github("aidenlab/straw/R")
# library(strawr)
# BiocManager::install("rhdf5")


## Modeled_regions_for_hg38.txt

centromeres = read.delim('/projects/b1154/referenceGenomes/hg38+TB40E/hg38.centromere.bed',
                         sep = '\t',
                         h = F,
                         stringsAsFactors = F)
head(centromeres)

######
## Z norm
MOCK_100kb_coolZ <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/MOCKmerg.100000_balanced.cool',
  sample_name = "MOCKc",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "black",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = TRUE)

INF_100kb_coolZ <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/INFmerg.100000_balanced.cool',
  sample_name = "INFc",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "red",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = TRUE)

siC_100kb_coolZ <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/siCINF72merg2.100000_balanced.cool',
  sample_name = "siCmerg",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "darkred",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = TRUE)

siS_100kb_coolZ <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/siSINF72merg2.100000_balanced.cool',
  sample_name = "siSmerg",
  centromeres = centromeres,
  resolution = 100000,
  balancing =  FALSE, # TRUE is the default
  colour = "darkorange",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = TRUE)

siA_100kb_coolZ <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/siAINF72merg2.100000_balanced.cool',
  sample_name = "siAmerg",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "darkorange",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = TRUE)
#######
#######
MOCK_100kb_cool <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/MOCKmerg.100000_balanced.cool',
  sample_name = "MOCKc",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "black",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = FALSE)

INF_100kb_cool <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/INFmerg.100000_balanced.cool',
  sample_name = "INFc",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "red",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = FALSE)

siC_100kb_cool <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/siCINF72merg2.100000_balanced.cool',
  sample_name = "siCmerg",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "darkred",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = FALSE)

siS_100kb_cool <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/siSINF72merg2.100000_balanced.cool',
  sample_name = "siSmerg",
  centromeres = centromeres,
  resolution = 100000,
  balancing =  FALSE, # TRUE is the default
  colour = "darkorange",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = FALSE)

siA_100kb_cool <- load_contacts(
  signal_path = '/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/hic/hicmatrices/cool/siAINF72merg2.100000_balanced.cool',
  sample_name = "siAmerg",
  centromeres = centromeres,
  resolution = 100000,
  balancing = FALSE, # TRUE is the default
  colour = "darkorange",
  # scale_bp = NULL,
  scale_bp = 1e9,
  scale_cis = TRUE,
  z_norm = FALSE)


######################################################
######################################################
# trans interactions
# trans_matrixplot(
#   exp1 = siC_100kb_cool,
#   chrom_up = "chr2", start_up = 0, end_up = 240e6,
#   chrom_down = "chrTB40E", start_down = 0, end_down = 280e6,
#   # colour_lim = c(0, 1000),
#   colour_lim = NULL,
#   rasterise = FALSE,
#   colour_bar = TRUE
# )
# trans_matrixplot(
#   exp1 = siC_100kb_cool,
#   chrom_up = "chr2", start_up = 0, end_up = 240e6,
#   chrom_down = "chrTB40E", start_down = 0, end_down = 280e6,
#   # colour_lim = c(0, 1000),
#   colour_lim = NULL,
#   rasterise = FALSE,
#   colour_bar = TRUE
# )

########################
### make sure znorm is ON
### 10(c) - 5(a) = 5 (pos)
### 5(c) - 10(a) = -5 (neg)

trans_matrixplot(
  exp1 = siC_100kb_coolZ, exp2 = siA_100kb_coolZ,
  chrom_down = "chr5", start_down = 0, end_down = 280e6,
  chrom_up = "chrTB40E", start_up = 0, end_up = 200e6,
  colour_lim = c(-5, 5),
  # colour_lim = NULL,
  rasterise = FALSE,
  colour_bar = TRUE
)

########################
#### 
### make sure znorm is OFF
RCP_out <- RCP(explist = list(MOCK_100kb_cool, INF_100kb_cool))
visualise(RCP_out, contrast = 1, metric = 'lfc')

RCP_outSvM <- RCP(explist = list(MOCK_100kb_cool, siS_100kb_cool))
visualise(RCP_outSvM, contrast = 1, metric = 'lfc')

RCP_outAvM <- RCP(explist = list(MOCK_100kb_cool, siA_100kb_cool))
visualise(RCP_outAvM, contrast = 1, metric = 'lfc')

RCP_outS <- RCP(explist = list(siC_100kb_cool, siS_100kb_cool))
visualise(RCP_outS, contrast = 1, metric = 'lfc')

RCP_outA <- RCP(explist = list(siC_100kb_cool, siA_100kb_cool))
visualise(RCP_outA, contrast = 1, metric = 'lfc')



######## 
############# genome wide analysis

# cm = chromosome_matrix(
#   list(MOCK_100kb_cool, INF_100kb_cool),
#   expected = "trans",
#   include_chr = "all",
#   exclude_chr = c("chrM", "chrY", "chrEBV"),
#   sort_chr = FALSE
#   )
# visualise(cm)
# expected = c("bins", "sums", "trans", "cis", "regress"),
# RCP_out <- RCP(explist = list(MOCK_100kb_cool, INF_100kb_cool),
#                chromsToUse = 'chr1')
# visualise(RCP_out)
# 
# # Plot RCP: combined


####################################################
####################################################
##########
# cisChrom_out <- cis_trans( list(MOCK_100kb_cool, INF_100kb_cool) )
# barplot(cisChrom_out$cis, names.arg = cisChrom_out$sample, ylim = c(0, 100) )
# abline(h = cisChrom_out$cis[1], col = 'red', lty = 3)
# abline(h = cisChrom_out$cis[2], col = 'red', lty = 3)
# 
# p_arms <- data.frame('chromosome' = centromeres[,1], 'start' = 0, 'end' = centromeres[,2])
# cisChrom_outP <- cis_trans( list(MOCK_100kb_cool, INF_100kb_cool) , bed = p_arms)
# barplot(cisChrom_outP$cis, names.arg = cisChrom_outP$sample, ylim = c(0, 100))
# abline(h = 90, col = 'red', lty = 3)
# abline(h = 93, col = 'red', lty = 3)
# 
# q_arms <- data.frame('chromosome' = centromeres[,1], 'start' = 0, 'end' = centromeres[,3])
# cisChrom_outQ <- cis_trans( list(MOCK_100kb_cool, INF_100kb_cool) , bed = q_arms)
# barplot(cisChrom_outQ$cis, names.arg = cisChrom_outQ$sample, ylim = c(0, 100))
# abline(h = 90, col = 'red', lty = 3)
# abline(h = 93, col = 'red', lty = 3)

## compartment score generated for input into saddle plot
MOCK_H3K4me3_seacr = read.delim('/projects/b1042/WalshLab/cdr5028/2024_figurefiles_CURRENT/cutrun/MOCK_H3K4me3_R1.seacr.peaks.relaxed.bed',
                                header = FALSE)

### find GENOVA Compartment score
CS_out = compartment_score(list(MOCK_100kb_cool, INF_100kb_cool),
                           bed = MOCK_H3K4me3_seacr)
CS_out2 = compartment_score(list(MOCK_100kb_cool, siC_100kb_cool),
                           bed = MOCK_H3K4me3_seacr)
CS_out_SKD = compartment_score(list(siC_100kb_cool, siS_100kb_cool),
                           bed = MOCK_H3K4me3_seacr)
CS_out_AKD = compartment_score(list(siC_100kb_cool, siA_100kb_cool),
                           bed = MOCK_H3K4me3_seacr)

# visualise(CS_out, chr = "chr5")

# saddle plot
saddle_out = saddle(list(MOCK_100kb_cool, INF_100kb_cool),
                    CS_discovery = CS_out,
                    bins = 10)
visualise(saddle_out)

saddle_out2 = saddle(list(MOCK_100kb_cool, siC_100kb_cool),
                    CS_discovery = CS_out2,
                    bins = 10)
visualise(saddle_out2)

saddle_out_SKD = saddle(list(siC_100kb_cool, siS_100kb_cool),
                    CS_discovery = CS_out_SKD,
                    bins = 10)
visualise(saddle_out_SKD)

saddle_out_AKD = saddle(list(siC_100kb_cool, siA_100kb_cool),
                        CS_discovery = CS_out_AKD,
                        bins = 10)
visualise(saddle_out_AKD)

