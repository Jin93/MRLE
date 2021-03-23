# MRLE

R package for Mendelian Randomization analysis for latent exposures (MRLE)<sup>1</sup> . The method uses GWAS summary-level association statistics of the outcome and K>=3 observable traits/biomarkers on a set of SNPs (instrumental variables) that are associated with at least two of the traits at specified thresholds of significance. The method uses an underlying structural equation model to describe causal paths between the SNPs, the latent exposure, the traits co-regulated by the exposure, and the outcome. A series of estimating functions are then constructed by equating the second-order sample moments of the summary-level association statistics with the corresponding theoretical moments, based on which inference for the model parameters can be done via Generalized Method of Moments.


## Installation

``` r
# install.packages("devtools")
devtools::install_github("Jin93/MRLE")
library(readr)
library(dplyr)
library(genio)
library(data.table)
library(MendelianRandomization) # for conducting test based on the IVW estimator
# setwd('~/MRLE/') # set path to the Github directory
```

## Example data analysis
### Testing the causal effect of chronic inflammation on rheumatoid arthritis (RA)

### Step 1: data preparation

#### Data sources (All GWAS samples are of European ancestry):
       1. GWAS summary data for RA   
    Okada, Y., Wu, D., Trynka, G., Raj, T., Terao, C., Ikari, K., Kochi, Y., Ohmura, K., Suzuki, A., Yoshida, S. and Graham, R.R., 2014. Genetics of rheumatoid arthritis contributes to biology and drug discovery. Nature, 506(7488), pp.376-381.
       2. GWAS summary data for CRP   
    GWAS on 320041 unrelated UK Biobank individuals
       3. GWAS summary data for IL1, IL6, IL8, TNP-alpha and MCP-1  
    Ahola-Olli, A.V., Würtz, P., Havulinna, A.S., Aalto, K., Pitkänen, N., Lehtimäki, T., Kähönen, M., Lyytikäinen, L.P., Raitoharju, E., Seppälä, I. and Sarin, A.P., 2017. Genome-wide association study identifies 27 loci influencing concentrations of circulating cytokines and growth factors. The American Journal of Human Genetics, 100(1), pp.40-50.



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
```

#### Before running the following R code, please unzip the following example data files and store them in data/ folder.
       1. data/sumdata_IL6.txt.zip  
       2. data/sumdata_IL8.txt.zip  
       3. data/sumdata_TNF.txt.zip  
       4. data/sumdata_MCP1.txt.zip  
       5. data/sumdata_CRP.txt.zip 

#### Please note that these input GWAS data should be pre-processed by a SNP filtering step, where the SNPs that have one or more of the following conditions are removed:
       1. MAF < 0.01  
       2. Effective sample size < 0.67 * 0.9 percentile
       3. Within the major histocompatibility complex (MHC) region ( 26Mb \~ 34Mb on chromosome 6)
       4. Alleles do not match those in the 1000 Genomes Project.

#### We then use the following code to merge GWAS summary statistics of the biomarkers:
``` r
sumdata1 = bigreadr::fread2(paste0('data/sumdata_',traitvec[1],'.txt'))
for (k in 2:K){
  tem = bigreadr::fread2(paste0('data/sumdata_',traitvec[k],'.txt'))
  tem = tem[,c('rsid', 'A1', 'A2', paste0(c('beta', 'se', 'P'), '.', traitvec[k]))]
  sumdata1 = sumdata_merge(sumdata1, tem, 'A1', 'A2', 'A1', 'A2', paste0('beta.',traitvec[k]))
}
sumbiomarkers = sumdata1; rm(sumdata1)
# write_delim(sumbiomarkers, file='data/sumbiomarkers.txt', delim='\t')
```

#### Load GWAS summary data for the outcome
``` r
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
```

#### Combine it with GWAS summary data for the biomarkers
```r
sumdata = sumdata_merge(sumbiomarkers, sumoutcome, 'A1', 'A2', paste0('A1.',outcome), paste0('A2.',outcome), paste0('beta.',outcome))

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
ind5 = ifelse(sumdata[[paste0('P.',traitvec[5])]] <= alpha2, 1, 0)
ind0 = ind1 + ind2 + ind3 + ind4 + ind5
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
    if ((chr == 5)&('rs12186596' %in% snpinfo$id)){
      snpinfo = snpinfo[-which(snpinfo$id == 'rs12186596'),]
    }
    snpinfo <- snpinfo[,c('chr','id','posg','pos','alt','ref')]
    write_bim(paste0("data/LD-clumping/snpinfo.",exposure,".",outcome,".",chr,".bim"),snpinfo)
  }
}
```

#### Download [1000 Genomes genotype data](https://www.internationalgenome.org/data/) and save it in data/ folder. Construct plink format genotype data for 1000G individudals of European ancestry.
```r
# mydir = 'data/LD-clumping/'
# name<-"select_intersect"
# for (chr in chr.list){
#   cmd1 <- paste(paste0('/dcl01/chatterj/data/jin/software/plink2 --bfile data/1000G_bial_nochild_europeans'),
#                 paste0('--extract data/LD-clumping/snpinfo.',exposure,'.', outcome,'.',chr,".bim"),
#                 #paste0('--keep data/1000G_Europeans.txt'),
#                 paste0('--make-bed'),
#                 paste0('--out data/LD-clumping/1000G.',exposure,'.',outcome,".",chr))
#   system(cmd1)
# }
```

#### Generate the summary data files which will be used as the input for LD Clumping
```r
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
```

#### LD clumping
```r
# parameters for LD clumping
r2=0.05; kb=1024
for (chr in c(1:22)){
  cmd1 <- paste0('/dcl01/chatterj/data/jin/software/plink --bfile data/LD-clumping/1000G.',exposure,'.',outcome,'.',chr,
                 ' --clump data/LD-clumping/snplist.',exposure,'.',outcome,'.',chr,'.rsid',
                 ' --clump-p1 ',alpha1,
                 ' --clump-p2 ',alpha1,
                 ' --clump-r2 ',r2,
                 ' --clump-kb ',kb,
                 ' --out data/LD-clumping/',exposure,'.',outcome,'.',chr)
  system(cmd1)
}
```

#### Create clumped SNP list
```r
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

#### The following R code can then be used to estimate the between-trait correlation matrix, which will be used later to estimate between-trait covariance matrices for the GWAS summary association statistics.

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
covresults = covinfo(sumstats, out.path=ldsc.out.path, ldsc.path, python.path, ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), maf.thr, mergeallele = TRUE, K)
#coherit.mat: Genetic correlation matrix
#ldscint.mat: LD score regression intercept matrix.
cor.mat = covresults$ldscint.mat
colnames(cor.mat) = rownames(cor.mat) = traitvec
save(cor.mat, file = 'cor.mat.RData')
```


### Step 5: test the causal effect of chronic inflammation on RA

#### Read clumped snp information.
``` r
############ read clumped snp info
sumdata = sumdata[sumdata[['rsid']] %in% snp.clumped,]
p.exposure = apply(sumdata[,paste0('P.',traitvec)], 1, second.min)
p.outcome = sumdata[[paste0('P.',outcome)]]


Bkind = list()
sign.level = ifelse(traitvec == 'CRP', alpha2, alpha1)
for (i in 1:K){
  Bkind[[i]] = ifelse(sumdata[[paste0('P.',traitvec[i])]]<= sign.level[i],1,0)
}
Bkinds = colSums(matrix(unlist(Bkind),K,nrow(sumdata),byrow=T))
ivind = which(Bkinds > 1) # index of the SNPs that are associated with at least one of the Bks
ivindk = lapply(1:K,FUN=function(x){which(Bkind[[x]]>0)}) # index of the SNPs that are associated with each biomarker
sumtable = sumdata[ivind,]

which.significant = which(sapply(1:K,function(x){length(ivindk[[x]])})>0)
traitvec = traitvec[which.significant]
K=length(traitvec)
thetak.sign = thetak.sign[which.significant]

trait.spec = NULL
for (trait in c(traitvec,outcome)){
  trait.spec = c(trait.spec, paste0(c("beta.", "se.", "P."), trait))
}
sumtable = sumtable[,c("Chr","Pos","rsid", "A1", "A2", trait.spec)]
```

#### Estimate between-biomarker covariance matrices.
```r
##### load between-biomarker genetic correlation matrix: cor.mat
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
```

#### Conduct hypothesis testing.
```r
### Data preparation for obtaining the IVW estimators.
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

output = mrle(sumtable, thetak.sign, Cov.mat, alpha0, ivindk)
```
This gives the output
```r
          M       P-value Reject Direction
IVW.IL6  25  0.000000e+00      1        -1
IVW.IL8  33  0.000000e+00      1         1
IVW.TNF  31  0.000000e+00      1         1
IVW.MCP1 17  0.000000e+00      1        -1
IVW.CRP   7 3.790054e-271      1         1
MRLE     55  9.313389e-23      1         1
```

###  References:
      1. Bulik-Sullivan, B., Finucane, H.K., Anttila, V., Gusev, A., Day, F.R., Loh, P.R., Duncan, L., Perry, J.R., Patterson, N., Robinson, E.B. and Daly, M.J., 2015. An atlas of genetic correlations across human diseases and traits. Nature genetics, 47(11), p.1236. [https://www.nature.com/articles/ng.3406.pdf?origin=ppub](https://www.nature.com/articles/ng.3406.pdf?origin=ppub)
      2. Bulik-Sullivan, B.K., Loh, P.R., Finucane, H.K., Ripke, S., Yang, J., Patterson, N., Daly, M.J., Price, A.L. and Neale, B.M., 2015. LD Score regression distinguishes confounding from polygenicity in genome-wide association studies. Nature genetics, 47(3), pp.291-295. [https://www.nature.com/articles/ng.3211](https://www.nature.com/articles/ng.3211)
      3. Jin, J., Qi, G., Yu, Z. and Chatterjee, N., 2021. Mendelian Randomization Analysis Using Multiple Biomarkers of an Underlying Common Exposure. bioRxiv.  [https://doi.org/10.1101/2021.02.05.429979](https://doi.org/10.1101/2021.02.05.429979)

