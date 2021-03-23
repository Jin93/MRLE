# MRLE

R package for Mendelian Randomization analysis for latent exposures (MRLE)<sup>1</sup> . The method uses GWAS summary-level association statistics of the outcome and K>=3 observable traits/biomarkers on a set of SNPs (instrumental variables) that are associated with at least two of the traits at specified thresholds of significance. The method uses an underlying structural equation model to describe causal paths between the SNPs, the latent exposure, the traits co-regulated by the exposure, and the outcome. A series of estimating functions are then constructed by equating the second-order sample moments of the summary-level association statistics with the corresponding theoretical moments, based on which inference for the model parameters can be done via Generalized Method of Moments.


## Contents

## Installation

``` r
# install.packages("devtools")
devtools::install_github("Jin93/MRLE")
library(readr)
library(dplyr)
library(genio)
library(data.table)
library(MendelianRandomization) # for conducting test based on the IVW estimator
```

## An example - testing the causal effect of chronic inflammation on rheumatoid arthritis (RA)

### Step 1: data preparation

First, generate the following matrices and save them in the data/ folder:
    1. GWAS summary-level association statistics for the biomarkers (with A1 and A2 matched across biomarkers)
    ``` r
    SNPID	Chr	Position(hg19)	A1	A2	OR(A1)	OR_95%CIlow	OR_95%CIup	pvalue_raw	pvalue	z	beta_scale	sigMsq
    rs3094315	1	752566	A	G	1.14	1.03	1.26	0.0093	0.0108228574562146	2.54836690974419	0.0244961759708397	9.24001490901232e-05
    rs3131972	1	752721	A	G	0.88	0.79	0.97	0.009	0.0146348973490792	-2.44128999686958	-0.0234668991857106	9.24001490901232e-05
    rs3131969	1	754182	A	G	0.85	0.75	0.96	0.0088	0.00985974473267805	-2.58070972418631	-0.0248070712626192	9.24001490901232e-05
    ```
    3. 
##### data sources:
    1. GWAS summary data for RA   
    Okada, Y., Wu, D., Trynka, G., Raj, T., Terao, C., Ikari, K., Kochi, Y., Ohmura, K., Suzuki, A., Yoshida, S. and Graham, R.R., 2014. Genetics of rheumatoid arthritis contributes to biology and drug discovery. Nature, 506(7488), pp.376-381.
    2. GWAS summary data for CRP   
    GWAS on 320041 unrelated UK Biobank individuals
    3. GWAS summary data for IL1, IL6, IL8, TNP-alpha and MCP-1  
    Ahola-Olli, A.V., Würtz, P., Havulinna, A.S., Aalto, K., Pitkänen, N., Lehtimäki, T., Kähönen, M., Lyytikäinen, L.P., Raitoharju, E., Seppälä, I. and Sarin, A.P., 2017. Genome-wide association study identifies 27 loci influencing concentrations of circulating cytokines and growth factors. The American Journal of Human Genetics, 100(1), pp.40-50.

All GWAS samples are of European ancestry.

``` r
library(MRLE)
traitvec = c('IL6','IL8','TNF','MCP1','CRP') # vector of K biomarkers
K = length(traitvec)
exposure = 'inflammation'
outcome = 'ra'
N.biomarker = c(8189,3526,3454,8293,320041) # GWAS sample size of the traits/biomarkers
names(N.biomarker) = traitvec
thetak.sign = rep(1,K) # directions of the causal effect of chronic inflammation on the inflammatory biomarkers
alpha0 = 5e-2
z0 = qnorm(p=alpha0/2,lower.tail=FALSE)
maf.thr = 0.01

####################################################################################################
##################################### Step 1: Data Preprocessing ###################################
################## Merge: match GWAS summary statistics of the biomarkers ##########################
######## SNP filtering - remove SNPs that have one or more of the following conditions
############### 1. MAF < 0.01
############### 2. effective sample size < 0.67 * 0.9 percentile
############### 3. within the major histocompatibility complex (MHC) region ( 26Mb \~ 34Mb on chromosome 6)
############### 4. alleles do not match those in the 1000 Genomes Project.
######## match A1 and A2 across biomarkers and save them in columns A1, A2
# setwd('~/MRLE/')
# GWAS data for cytokines
sumdata1 = bigreadr::fread2('data/sumdata_cytokines.txt')
# GWAS data for CRP
sumdata2 = bigreadr::fread2('data/sumdata_crp.txt')
sumbiomarkers = merge(sumdata1,sumdata2,by='rsid')
#### match A1 and A2 of CRP GWAS data (columns A1.CRP, A2.CRP) to A1 and A2 of the rest of the biomarkers (columns A1, A2)
#### re
## scenario 1: A1 A2 flipped between CRP and the other biomarkers
inds.flipped = ((sumbiomarkers[["A1"]] == sumbiomarkers[[paste0('A2.CRP')]])&(sumbiomarkers[["A2"]] == sumbiomarkers[[paste0('A1.CRP')]]))
## match to biomarkers A1 A2
if (sum(inds.flipped) > 0){
  sumbiomarkers[[paste0("beta.CRP")]][inds.flipped] = -sumbiomarkers[[paste0("beta.CRP")]][inds.flipped]
  sumbiomarkers[[paste0("A1.CRP")]][inds.flipped] = sumbiomarkers[['A1']][inds.flipped]
  sumbiomarkers[[paste0("A2.CRP")]][inds.flipped] = sumbiomarkers[['A2']][inds.flipped]
}
## scenario 2: A1 A2 flipped; modifiable strand-ambiguous:
combine.alleles = function(x) paste(x,collapse='')
alleles = cbind(sumbiomarkers[["A1"]], sumbiomarkers[["A2"]], sumbiomarkers[[paste0('A1.CRP')]], sumbiomarkers[[paste0('A2.CRP')]])
alleles = apply(alleles,1,combine.alleles)
inds.ambiguous= which(alleles %in% c('ACGT','AGCT','TCGA','TGCA','CATG','CTAG','GATC','GTAC'))
if (length(inds.ambiguous) > 0){
  sumbiomarkers[[paste0("beta.CRP")]][inds.ambiguous] = -sumbiomarkers[[paste0("beta.CRP")]][inds.ambiguous]
  sumbiomarkers[[paste0("A1.CRP")]][inds.ambiguous] = sumbiomarkers[['A1']][inds.ambiguous]
  sumbiomarkers[[paste0("A2.CRP")]][inds.ambiguous] = sumbiomarkers[['A2']][inds.ambiguous]
}
###### now decide which ones to keep
inds.ambiguous.keep = alleles %in% c('ACTG','AGTC','TCAG','TGAC','CAGT','CTGA','GACT','GTCA')
if (sum(inds.ambiguous.keep) > 0){
  sumbiomarkers[[paste0("A1.CRP")]][inds.ambiguous.keep] = sumbiomarkers[['A1']][inds.ambiguous.keep]
  sumbiomarkers[[paste0("A2.CRP")]][inds.ambiguous.keep] = sumbiomarkers[['A2']][inds.ambiguous.keep]
}
inds.match.keep = ((sumbiomarkers[["A1"]] == sumbiomarkers[[paste0('A1.CRP')]])&(sumbiomarkers[["A2"]] == sumbiomarkers[[paste0('A2.CRP')]]))
if (sum(inds.ambiguous.keep)>0){
  inds.keep = c(which(inds.ambiguous.keep),which(inds.match.keep))
}
if (sum(inds.ambiguous.keep)==0){
  inds.keep = which(inds.match.keep)
}
sumbiomarkers = sumbiomarkers[inds.keep,]
sumbiomarkers = sumbiomarkers[,-which(colnames(sumbiomarkers) %in% c('A1.CRP','A2.CRP'))]
write_delim(sumbiomarkers, file='data/sumbiomarkers.txt', delim='\t')

# ------------- Load GWAS summary data for the outcome
sumoutcome = bigreadr::fread2('data/raw_sumdata_ra.txt')
sumoutcome$SE = log(sumoutcome$`OR(A1)`)/sumoutcome$z
sumoutcome = sumoutcome[,c('Chr','Position(hg19)','SNPID','A1','A2','OR(A1)','SE','pvalue')]
sumoutcome = sumoutcome[complete.cases(sumoutcome),]
colnames(sumoutcome) = c('CHR','POS','SNP','A1','A2','beta','se','p')
sumoutcome = sumoutcome[,c('SNP','A1','A2','beta','se','p')]
colnames(sumoutcome) = c("MarkerName","A1","A2","beta","se","P")
# Remove SNPs with sample size < 0.67 * (90 percentile)
if ('N' %in% colnames(sumoutcome)){
  sumoutcome = sumoutcome %>% filter(N>0.67*quantile(N,0.9))
}
sumoutcome$A1 = toupper(sumoutcome$A1)
sumoutcome$A2 = toupper(sumoutcome$A2)
colnames(sumoutcome)[which(colnames(sumoutcome) == 'MarkerName')] = 'rsid'
sumoutcome = sumoutcome[,c('rsid', 'A1', 'A2', 'beta', 'se', 'P')]
names(sumoutcome)[2:ncol(sumoutcome)] = paste0(names(sumoutcome)[2:ncol(sumoutcome)],'.',outcome)
sumoutcome$rsid = as.character(sumoutcome$rsid)
# ------------- combine GWAS summary data for biomarkers with GWAS summary data for the outcome
sumdata = sumbiomarkers %>% inner_join(sumoutcome, by = 'rsid')
sumdata = sumdata[complete.cases(sumdata),]

#### match A1 and A2 of the outcome (columns A1.outcome, A2.outcome) to A1 and A2 of the biomarkers (columns A1, A2)
#### re
## scenario 1: A1 A2 flipped between the biomarker and the outcome
inds.flipped = ((sumdata[["A1"]] == sumdata[[paste0('A2.',outcome)]])&(sumdata[["A2"]] == sumdata[[paste0('A1.',outcome)]]))
## match to biomarkers A1 A2
if (sum(inds.flipped) > 0){
  sumdata[[paste0("beta.",outcome)]][inds.flipped] = -sumdata[[paste0("beta.",outcome)]][inds.flipped]
  sumdata[[paste0("A1.",outcome)]][inds.flipped] = sumdata[['A1']][inds.flipped]
  sumdata[[paste0("A2.",outcome)]][inds.flipped] = sumdata[['A2']][inds.flipped]
}
## scenario 2: A1 A2 flipped between the biomarker and the outcome; modifiable strand-ambiguous:
alleles = cbind(sumdata[["A1"]], sumdata[["A2"]], sumdata[[paste0('A1.',outcome)]], sumdata[[paste0('A2.',outcome)]])
alleles = apply(alleles,1,combine.alleles)
inds.ambiguous= which(alleles %in% c('ACGT','AGCT','TCGA','TGCA','CATG','CTAG','GATC','GTAC'))
if (length(inds.ambiguous) > 0){
  sumdata[[paste0("beta.",outcome)]][inds.ambiguous] = -sumdata[[paste0("beta.",outcome)]][inds.ambiguous]
  sumdata[[paste0("A1.",outcome)]][inds.ambiguous] = sumdata[['A1']][inds.ambiguous]
  sumdata[[paste0("A2.",outcome)]][inds.ambiguous] = sumdata[['A2']][inds.ambiguous]
}
###### now decide which ones to keep
inds.ambiguous.keep = alleles %in% c('ACTG','AGTC','TCAG','TGAC','CAGT','CTGA','GACT','GTCA')
if (sum(inds.ambiguous.keep) > 0){
  sumdata[[paste0("A1.",outcome)]][inds.ambiguous.keep] = sumdata[['A1']][inds.ambiguous.keep]
  sumdata[[paste0("A2.",outcome)]][inds.ambiguous.keep] = sumdata[['A2']][inds.ambiguous.keep]
}
inds.match.keep = ((sumdata[["A1"]] == sumdata[[paste0('A1.',outcome)]])&(sumdata[["A2"]] == sumdata[[paste0('A2.',outcome)]]))
if (sum(inds.ambiguous.keep)>0){
  inds.keep = c(which(inds.ambiguous.keep),which(inds.match.keep))
}
if (sum(inds.ambiguous.keep)==0){
  inds.keep = which(inds.match.keep)
}
sumdata = sumdata[inds.keep,]

# Extract trait-specific column names
trait.spec = NULL
for (trait in c(traitvec,outcome)){
  trait.spec = c(trait.spec, paste0(c("beta.", "se.", "P."), trait))
}
sumdata = sumdata[,c("loc","rsid", "A1", "A2", trait.spec)]
sumdata <- sumdata %>%
  mutate(position = strsplit(loc, split = ':')) %>%
  mutate(Chr = as.integer(sapply(position, function(x) x[1])), Pos = as.integer(sapply(position, function(x) x[2])))
sumdata = sumdata[,c("loc","Chr","Pos","rsid", "A1", "A2", trait.spec)]
```


### Step 2: select instrumental variables (IVs)
``` r
alpha1 = 1e-3 # SNP significance level of the association with the rest of the biomarkers
alpha2 = 5e-6 # SNP significance level of the association with CRP
z = qnorm(p=alpha1/2,lower.tail=FALSE)
z2 = qnorm(p=alpha2/2,lower.tail=FALSE)
ind1 = ifelse(sumdata[[paste0('P.',traitvec[1])]] <= alpha1, 1, 0)
ind2 = ifelse(sumdata[[paste0('P.',traitvec[2])]] <= alpha1, 1, 0)
ind3 = ifelse(sumdata[[paste0('P.',traitvec[3])]] <= alpha1, 1, 0)
ind4 = ifelse(sumdata[[paste0('P.',traitvec[4])]] <= alpha1, 1, 0)
ind5 = ifelse(sumdata[[paste0('P.',traitvec[5])]] <= alpha1, 1, 0)
ind6 = ifelse(sumdata[[paste0('P.',traitvec[6])]] <= alpha2, 1, 0)
ind0 = ind1 + ind2 + ind3 + ind4 + ind5 + ind6
sumdata = sumdata[ind0 >= 2,]

conf.rsid = readRDS(paste0('data/rsid_confounding.rds'))
confounders = c('bmi','sbp','diabetes','smoking','alcohol','hdl','ldl')
conf.rsid = unlist(conf.rsid[c(1:7)])
### do not remove SNPs associated with the confounders: o.w. there is no significant SNP left for CRP
sumdata = sumdata[!(sumdata[['rsid']] %in% conf.rsid),]
write_delim(sumdata, file = 'data/sumdata.txt', delim='\t')
```



### Step 3: LD clumping
``` r
# ----------------------- Create SNP list -----------------------
chr.list = numeric()
for (chr in 1:22){
  sum.chrspecific = sumdata[sumdata[['Chr']] == chr,]
  if (nrow(sum.chrspecific) > 0){
    chr.list = c(chr.list,chr)
    snpinfo <- make_bim(n = nrow(sum.chrspecific))
    # add the "chr" prefix to the chromosome values so that we recognize them when we see them later.
    snpinfo$chr <- chr
    # Make SNP IDs look like "rs" IDs
    snpinfo$id <- sum.chrspecific[['rsid']]
    snpinfo$posg <- 0
    snpinfo$pos <- 0
    snpinfo$ref <- 0
    snpinfo$alt <- 0
    ### delete the one duplicated SNP
    if (chr == 5){
      snpinfo = snpinfo[-which(snpinfo$id == 'rs12186596'),]
    }
    snpinfo <- snpinfo[,c('chr','id','posg','pos','alt','ref')]
    write_bim(paste0("data/LD-clumping/snpinfo.",exposure,".",outcome,".",chr,".bim"),snpinfo)
  }
}

# ----------------------- Select SNPs that are present in 1000 Genome reference panel -----------------------
# mydir = 'data/LD-clumping/'
# name<-"select_intersect"
# for (chr in chr.list){
#   filen<-paste0(mydir, "sh/", name, '_',exposure,'.', outcome,".", chr, ".sh")
#   file.create(filen)
#   zz <- file(filen, "w")
#   cat("#$ -cwd", "", file = zz, sep = "\n")
#   cat(paste0('#$ -o ',mydir,'logfile'), file = zz, sep = "\n")
#   cat(paste0('#$ -e ',mydir,'logfile'), file = zz, sep = "\n")
#   cat(paste0('cd ',mydir), "\n", file=zz)
#   cat("\n", file=zz)
#   cmd1 <- paste(paste0('plink2 --bfile 1000G_bial_nochild'),
#                 paste0('--extract data/LD-clumping/snpinfo.',exposure,'.', outcome,'.',chr,".bim"),
#                 paste0('--keep 1000G_Europeans.txt'),
#                 paste0('--make-bed'),
#                 paste0('--out data/LD-clumping/1000G.',exposure,'.',outcome,".",chr))
#   system(cmd1)
# }


# ----------------------- prepare the summary data as the input for LD Clumping -----------------------
second.min = function(x) sort(x,decreasing=F)[2]
for (chr in 1:22){
  sum.chrspecific = sumdata[sumdata[['Chr']] == chr,]
  temfile = paste0('data/LD-clumping/1000G.',exposure,".",outcome,'.',chr,'.bim')
  if(file.exists(temfile)){
    sump <- read_delim(temfile, delim='\t',col_names = F)
    sump <- sum.chrspecific[(sum.chrspecific[['rsid']] %in% sump[[2]]),c("rsid",paste0('P.',traitvec))]
    sum.p <- apply(sump[,paste0('P.',traitvec)], 1, second.min)
    sump <- cbind(sump$rsid,sum.p)
    colnames(sump) <- c('SNP','P')
    fwrite(sump,file=paste0('data/LD-clumping/snplist.',exposure,'.',outcome,'.',chr,'.rsid'),row.names = F, quote = F, sep=' ')
    print(chr)
  }
}


# ----------------------- LD Clumping -----------------------
# parameters for LD clumping
r2=0.05; kb=1024
name<-"ldclump"
for (chr in c(1:22)){
  filen<-paste0(mydir, "sh/", name, '_', exposure,'.', outcome,'.', chr, ".sh")
  file.create(filen)
  zz <- file(filen, "w")
  cat("#$ -cwd", "", file = zz, sep = "\n")
  #cat(paste0('#$ -o data/LD-clumping/logfile'), file = zz, sep = "\n")
  #cat(paste0('#$ -e data/LD-clumping/logfile'), file = zz, sep = "\n")
  cat(paste0('cd data/LD-clumping/'), "\n", file=zz)
  cat("\n", file=zz)
  cmd1 <- paste0('/dcl01/chatterj/data/jin/software/plink --bfile data/LD-clumping/1000G.',exposure,'.',outcome,'.',chr,
                 ' --clump data/LD-clumping/snplist.',exposure,'.',outcome,'.',chr,'.rsid',
                 ' --clump-p1 ',alpha1,
                 ' --clump-p2 ',alpha1,
                 ' --clump-r2 ',r2,
                 ' --clump-kb ',kb,
                 ' --out data/LD-clumping/',exposure,'.',outcome,'.',chr)
  system(cmd1)
}
# ------------------ create clumped SNP list ------------------
select.type = setting1[set1,'select.type']
#load(paste0('/dcl01/chatterj/data/jin/mr/sumdat/sumdata.inflammation-rhema-adj3-',select.type,'-p1=',alpha2,'-p2=',alpha1,'.RData'))
chr = 1
sumclumped = NULL
temfile = paste0('data/LD-clumping/',exposure,'.', outcome,'.', chr, '.clumped')
if(file.exists(temfile)){
  snp.clumped = read_delim(temfile,delim=' ',progress=F) # data.frame
  colnames(snp.clumped) = trimws(colnames(snp.clumped))
  snp.clumped = trimws(snp.clumped$SNP)
  sumclumped = sumdata[((sumdata$Chr == chr)&(sumdata$rsid %in% snp.clumped)),]
}
for (chr in 2:22){
  temfile = paste0('data/LD-clumping/',exposure,'.', outcome,'.', chr, '.clumped')
  if(file.exists(temfile)){
    temp = read_delim(temfile,delim=' ',progress=F) # data.frame
    colnames(temp) = trimws(colnames(temp))
    temp = trimws(temp$SNP)
    sumclumpedtemp = sumdata[((sumdata$Chr == chr)&(sumdata$rsid %in% temp)),]
    sumclumped = rbind(sumclumped,sumclumpedtemp)
  }
}
snp.clumped = sumclumped[['rsid']]
save(snp.clumped,file=paste0('data/LD-clumping/snp.clumped.',exposure, outcome,'.RData'))

```


### Step 4: estimate between-trait covariance matrix of the GWAS summary statistics

Please install the [LDSC](https://github.com/bulik/ldsc) software that performs LD score regression<sup>2,3</sup> to estimate heritability and genetic correlation. Python and all the Python packages are required for running LDSC. LDSC recommends installing Anaconda 2 which comes with all the required Python packages. LD score files and HapMap 3 SNP list provided by LD Hub can be downloaded and unzipped by running the following command and put them in the ldsc/ folder.

``` r
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 eur_w_ld_chr.tar.bz2
tar -xvf eur_w_ld_chr.tar
bunzip2 w_hm3.snplist.bz2
```

The following R code can then be used to estimate the between-trait correlation matrix, which will be used later to estimate between-trait covariance matrices for the GWAS summary association statistics.

``` r
ldsc.out.path="data/ldsc.out/"
ldsc.path = "data/ldsc/"
python.path = "~/anaconda2/bin/python" # change to the path to your python software
ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/") # default, does not need to be specified
# create input list for covinfo_inflammation function:
sumstats = list()
for (k in 1:K){
  sumstats[[k]] = sumbiomarkers[,c('rsid','Chr','Pos','A2','A1',paste0(c('beta','se','P'),'.',traitvec[k]))]
  colnames(sumstats[[k]]) = c('rsid','Chr','Pos','A2','A1','beta','se','P')
  sumstats[[k]]$N = N.biomarker[k]
}
covresults = covinfo_inflammation(sumstats, out.path=ldsc.out.path, ldsc.path, python.path, ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), maf.thr, mergeallele = TRUE, K, screening = F)
#coherit.mat: Genetic covariance matrix
#ldscint.mat: LD score regression intercept matrix.
cor.mat = covresults$ldscint.mat
colnames(cor.mat) = rownames(cor.mat) = traitvec
save(cor.mat, file = 'cor.mat.RData')
```


### Step 5: test the causal effect of chronic inflammation on RA

``` r
############ read clumped snp info
sumdata = sumdata[sumdata[['rsid']] %in% snp.clumped,]
p.exposure = apply(sumdata[,paste0('P.',traitvec)], 1, second.min)
p.outcome = sumdata[[paste0('P.',outcome)]]


######### Filtering: select the SNPs that reach genome-wide significance for Bks:
Bkind = list()
sign.level = ifelse(traitvec == 'CRP', alpha2, alpha1)
for (i in 1:K){
  Bkind[[i]] = ifelse(sumdata[[paste0('P.',traitvec[i])]]<= sign.level[i],1,0)
}
Bkinds = colSums(matrix(unlist(Bkind),K,nrow(sumdata),byrow=T))
ivind = which(Bkinds > 1) # index of the SNPs that are associated with at least one of the Bks
ivindk = lapply(1:K,FUN=function(x){which(Bkind[[x]]>0)}) # index of the SNPs that are associated with each Bk
sumtable = sumdata[ivind,]

which.significant = which(sapply(1:K,function(x){length(ivindk[[x]])})>0)
traitvec = traitvec[which.significant]
K=length(traitvec)
thetak.sign = thetak.sign[which.significant]

##### may exclude some biomarkers:
trait.spec = NULL
for (trait in c(traitvec,outcome)){
  trait.spec = c(trait.spec, paste0(c("beta.", "se.", "P."), trait))
}
sumtable = sumtable[,c("Chr","Pos","rsid", "A1", "A2", trait.spec)]
##### load between-biomarker correlation matrix: cor.mat
load(paste0('data/cor.mat.RData'))
cov.mat = matrix(0,K+1,K+1)
cov.mat[2:(K+1),2:(K+1)] = cor.mat[traitvec,traitvec]
diag(cov.mat) = 1
cor.mat = cov.mat # adding the correlation between biomarkers and the outcome
Cov.mat = list()
for (ni in 1:nrow(sumtable)){
  beta.sd = diag(sapply(c(outcome,traitvec), function(x) {(sumtable[[paste0('se.',x)]][ni])}))
  Cov.mat[[ni]] = beta.sd %*% cor.mat %*% beta.sd
}
beta.sd = diag(sapply(c(outcome,traitvec), function(x) median(sumtable[[paste0('se.',x)]])))
c.mat = beta.sd %*% cor.mat %*% beta.sd

########## Data Analysis
### data for ivw analysis
Bkind = list()
ivw.ind = list()
sign.level = ifelse(traitvec == 'CRP', alpha2, alpha1)
for (i in 1:K){
  Bkind[[i]] = ifelse(sumtable[[paste0('P.',traitvec[i])]]<= sign.level[i],1,0)
}
Bkinds = colSums(matrix(unlist(Bkind),K,nrow(sumtable),byrow=T))
ivind = which(Bkinds > 1) # index of the SNPs that are associated with at least one of the Bks
ivindk = lapply(1:K,FUN=function(x){which(Bkind[[x]]>0)}) # index of the SNPs that are associated with each Bk
names(ivindk) = traitvec
sumtable0=sumtable

set.seed(2020)
output = mrle(sumtable, thetak.sign, Cov.mat, alpha0)
```


####  Reference:
      1. Jin, J., Qi, G., Yu, Z. and Chatterjee, N., 2021. Mendelian Randomization Analysis Using Multiple Biomarkers of an Underlying Common Exposure. bioRxiv.  [https://doi.org/10.1101/2021.02.05.429979](https://doi.org/10.1101/2021.02.05.429979)

