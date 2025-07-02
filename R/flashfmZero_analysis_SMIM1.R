library(readr)
library(flashfmZero)

# read in GWAS summary statistics files for SMIM1
gwasfiles <- list.files(pattern = "_SMIM1\\.tsv\\.gz$")
traits <- unlist(strsplit(gwasfiles,"_SMIM1.tsv.gz"))
obsgwas <- vector("list",length(gwasfiles))
for(i in 1:length(gwasfiles)) obsgwas[[i]] <- as.data.frame(read_tsv(gwasfiles[i], col_names=TRUE))
names(obsgwas) <- traits


# read in trait covariance matrix
covY <- read_tsv("trait_covariance.tsv.gz", col_names = TRUE)
covY  <- as.data.frame(covY)
rownames(covY) <- covY[,1]
covY  <- as.matrix(covY[,-1])


# harmonise GWAS summary statistics
# gwas files are in two formats (different column names) so need to do in stages:
# identify gwas files with column names containing "INT" - these are the 36 classical blood cell traits
colnames <- lapply(obsgwas,names)
tmp <- lapply(colnames,function(x) grep("INT",x))
ind <- which(sapply(tmp, function(x) length(x)>0))
obsgwas36 <- obsgwas[ind]

# the #VARIANT column is chosen for snpID because there are missing values in the ID_dbSNP49 column
obsgwas36 <- harmoniseGWAS(obsgwas36, minMAF=0.005, minINFO=0.4, beta_colname="EFFECT_INT", 
se_colname="SE_INT", snpID_colname="#VARIANT", EA_colname="ALT", NEA_colname="REF", 
EAfreq_colname="ALT_FREQ_INT", BP_colname="BP", pvalue_colname="MLOG10P_INT", INFO_colname="INFO_INT")
# convert -log10p to p
for(i in 1:36) obsgwas36[[i]]$p_value <- 10^(-obsgwas36[[i]]$p_value)

# 63 non-classical blood cell traits
obsgwas63 <- obsgwas[(1:99)[-ind]]

obsgwas63 <- harmoniseGWAS(obsgwas63, minMAF=0.005, minINFO=0.4, beta_colname="EFFECT", 
se_colname="SE", snpID_colname="#VARIANT", EA_colname="ALT", NEA_colname="REF", 
EAfreq_colname="ALT_FREQ", BP_colname="BP", pvalue_colname="P", INFO_colname="INFO")

# combine all gwas and harmonise
obsgwas <- c(obsgwas36,obsgwas63)

obsgwas <- harmoniseGWAS(obsgwas, minMAF=0.005, minINFO=0.4, beta_colname="beta", 
se_colname="SE", snpID_colname="rsID", EA_colname="EA", NEA_colname="NEA", 
EAfreq_colname="EAF", BP_colname="BP", pvalue_colname="p_value", INFO_colname="INFO")


# factor analysis of observed traits
facheck <- psych::fa.parallel(covY, fm="ml", n.obs=43000,fa="fa", plot=FALSE)
faY <- psych::fa(r=covY, nfactors=facheck$nfact, fm="ml", rotate="varimax")

# calculate latent factor contributions to the traits
latent_terms <- flashfmZero::factor_contributions(faY$loadings) 
Contributions <- latent_terms[[1]] 
FAloadings <- latent_terms[[2]]

# Calculate latent factor GWAS summary statistics:
covY <- covY[rownames(FAloadings),rownames(FAloadings)]
obsGWAS <- obsgwas[rownames(FAloadings)]
latent_gwas <- flashfmZero::latentGWAS(obsGWAS = obsGWAS,
covY = covY, L = FAloadings, beta_colname="beta", se_colname="SE", snpID_colname="rsID", EAF = obsGWAS[[1]]$EAF)


# Identify the latent factors with a genome-wide significant (p<5✕10-8) signal in the region, 
# and if there are no such latent factors, then stop the analysis
check <- sapply(latent_gwas,function(x) any(x$p_value<5E-8))
fm_traits_latent <- names(check[check==TRUE])
if(length(fm_traits_latent)>0) sig_gwas_latent <- latent_gwas[fm_traits_latent]
try(if(length(fm_traits_latent)==0) stop("No latent factor has a signal"))


# read in LD matrix
ld_matrix <- read_tsv("ld_matrix.tsv.gz", col_names = TRUE)
ld_matrix  <- as.data.frame(ld_matrix)
rownames(ld_matrix) <- ld_matrix[,1]
corX  <- as.matrix(ld_matrix[,-1])


# Re-name the variant IDs in the latent factor GWAS summary statistics and the LD matrix 
# so that they start with a character and do not contain a “:”.
corXnames <- paste0("chr",rownames(corX))
rownames(corX) <- colnames(corX) <- corXnames
sig_gwas_latent_names <- paste0("chr",sig_gwas_latent[[1]]$rsID)
for(i in 1:length(sig_gwas_latent)) sig_gwas_latent[[i]]$rsID <- sig_gwas_latent_names
names(sig_gwas_latent) <- fm_traits_latent


if(length(fm_traits_latent)==1){
   FM_latent <-JAMdwithGroups(sig_gwas_latent[[1]], N=N, corX,   save.path="tmpDIR",r2.minmerge=0.8, cred=0.95)
   FM_latent_CS <- FM_latent$CS
 } else {
     FM_latent <-FLASHFMZEROwithJAMd(sig_gwas_latent, corX, 
  N = 43000, save.path="tmpDIR", NCORES=length(sig_gwas_latent), r2.minmerge=0.8)                                                 
  FM_latent_CS <-allcredsetsPP(FM_latent$mpp.pp, cred=.95) 
  names(FM_latent_CS$fm) <- fm_traits_latent
  names(FM_latent_CS$flashfm) <- fm_traits_latent
 }


