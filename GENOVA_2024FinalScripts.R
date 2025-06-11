################################################################################################
################################################################################################
######## Previously combined and balanced cooler files were used as input for GENOVA
################################################################################################
## GENOVA
# install.packages("remotes")
# remotes::install_github("robinweide/GENOVA")
library(GENOVA)

# install.packages("BiocManager")
# BiocManager::install("rtracklayer")
# library(rtracklayer)
# 
# remotes::install_github("aidenlab/straw/R")
# library(strawr)
# BiocManager::install("rhdf5")


centromeres = read.delim('~/hg38.centromere.bed',
                         sep = '\t',
                         h = F,
                         stringsAsFactors = F)
head(centromeres)

######
## Z norm
MOCK_100kb_coolZ <- load_contacts(
  signal_path = '~/MOCKmerg.100000_balanced.cool',
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
  signal_path = '~/INFmerg.100000_balanced.cool',
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
  signal_path = '~/siCINF72merg2.100000_balanced.cool',
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
  signal_path = '~/siSINF72merg2.100000_balanced.cool',
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
  signal_path = '~/siAINF72merg2.100000_balanced.cool',
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
  signal_path = '~/MOCKmerg.100000_balanced.cool',
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
  signal_path = '~/INFmerg.100000_balanced.cool',
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
  signal_path = '~/siCINF72merg2.100000_balanced.cool',
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
  signal_path = '~/siSINF72merg2.100000_balanced.cool',
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
  signal_path = '~/siAINF72merg2.100000_balanced.cool',
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


## compartment score generated for input into saddle plot
MOCK_H3K4me3_seacr = read.delim('~/MOCK_H3K4me3.seacr.peaks.relaxed.bed',
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

