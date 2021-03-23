#' Combine alleles in two files into one string
combine.alleles = function(x) paste(x,collapse='')

#' Merge GWAS summary-level association statistics of two traits
#'
#' This function merge GWAS summary-level association statistics of two traits according to rsid.
#' It then matches alleles of the second trait to those of the first trait with the following steps:
#' 1. Align flipped alleles
#' 2. Align alleles with modifiable strand-ambiguous
#' @param sumdata1 a matrix storing GWAS summary-level association statistics of the first trait (name: traitname1)
#'  that contain columns 'loc', 'Chr', 'Pos', 'rsid', A1.col.1, A2.col.1, and paste0(c('beta', 'se', 'P'), '.', traitname1).
#'  Here, A1.col.1, A2.col.1, beta, se and P represent allele1, allele2, the estimated effect size
#'  (corresponding to allele1), standard error and p-value obtained from GWAS.
#' @param sumdata2 a matrix storing GWAS summary-level association statistics of the first trait (name: traitname2)
#'  with columns A1.col.2, A2.col.2, 'rsid', and paste0(c('beta', 'se', 'P'), '.', traitname2).
#'  Here, A1.col.1, A2.col.1, beta, se and P represent allele1, allele2, the estimated effect size
#'  (corresponding to allele1), standard error and p-value obtained from GWAS.
#' @param A1.col.1 column name in sumdata1 that corresponds to allele1 of the summary data for the first trait.
#' @param A2.col.1 column name in sumdata1 that corresponds to allele2 of the summary data for the first trait.
#' @param A1.col.2 column name in sumdata1 that corresponds to allele1 of the summary data for the second trait.
#' @param A2.col.2 column name in sumdata1 that corresponds to allele2 of the summary data for the second trait.
#' @param beta.col.2 column name in sumdata2 that corresponds to beta of the summary data for the second trait.
#' @return Merged summary data matrix
#' @export
sumdata_merge = function(sumdata1, sumdata2, A1.col.1, A2.col.1, A1.col.2, A2.col.2, beta.col.2){
  colnames(sumdata1)[which(colnames(sumdata1) == A1.col.1)] = 'A1'
  colnames(sumdata1)[which(colnames(sumdata1) == A2.col.1)] = 'A2'
  colnames(sumdata2)[which(colnames(sumdata2) == A1.col.2)] = 'A1'
  colnames(sumdata2)[which(colnames(sumdata2) == A2.col.2)] = 'A2'
  sumdata = merge(sumdata1,sumdata2,by='rsid')
  ## scenario 1: A1 A2 flipped
  inds.flipped = ((sumdata[["A1.x"]] == sumdata[['A2.y']])&(sumdata[["A2.x"]] == sumdata[[paste0('A1.y')]]))
  if (sum(inds.flipped) > 0){
    sumdata[[beta.col.2]][inds.flipped] = -sumdata[[beta.col.2]][inds.flipped]
    sumdata[["A1.y"]][inds.flipped] = sumdata[['A1.x']][inds.flipped]
    sumdata[["A2.y"]][inds.flipped] = sumdata[['A2.x']][inds.flipped]
  }
  ## scenario 2: Modifiable strand-ambiguous:
  alleles = cbind(sumdata[["A1.x"]], sumdata[["A2.x"]], sumdata[[paste0('A1.y')]], sumdata[[paste0('A2.y')]])
  alleles = apply(alleles,1,combine.alleles)
  inds.ambiguous= which(alleles %in% c('ACGT','AGCT','TCGA','TGCA','CATG','CTAG','GATC','GTAC'))
  if (length(inds.ambiguous) > 0){
    sumdata[[beta.col.2]][inds.ambiguous] = -sumdata[[beta.col.2]][inds.ambiguous]
    sumdata[["A1.y"]][inds.ambiguous] = sumdata[['A1.x']][inds.ambiguous]
    sumdata[["A2.y"]][inds.ambiguous] = sumdata[['A2.x']][inds.ambiguous]
  }
  ###### Decide which ones to keep
  inds.ambiguous.keep = alleles %in% c('ACTG','AGTC','TCAG','TGAC','CAGT','CTGA','GACT','GTCA')
  if (sum(inds.ambiguous.keep) > 0){
    sumdata[[paste0("A1.y")]][inds.ambiguous.keep] = sumdata[['A1.x']][inds.ambiguous.keep]
    sumdata[[paste0("A2.y")]][inds.ambiguous.keep] = sumdata[['A2.x']][inds.ambiguous.keep]
  }
  inds.match.keep = ((sumdata[["A1.x"]] == sumdata[[paste0('A1.y')]])&(sumdata[["A2.x"]] == sumdata[[paste0('A2.y')]]))
  if (sum(inds.ambiguous.keep)>0){
    inds.keep = c(which(inds.ambiguous.keep),which(inds.match.keep))
  }
  if (sum(inds.ambiguous.keep)==0){
    inds.keep = which(inds.match.keep)
  }
  sumdata = sumdata[inds.keep,]
  sumdata = sumdata[,-which(colnames(sumdata) %in% c('A1.y','A2.y'))]
  colnames(sumdata)[which(colnames(sumdata) == 'A1.x')] = 'A1'
  colnames(sumdata)[which(colnames(sumdata) == 'A2.x')] = 'A2'
  sumdata
}


#' Data preprocessing
#'
#' @param sumstats a list of matrices storing the GWAS summary-level association statistics of the biomarkers.
#' Each matrix contains the following columns: rsid, Chr, Pos, A2, A1, beta, se, P, N.
#' @param maf.thr MAF threshold used for filtering out SNPs that have MAF<maf.thr.
#' @param K number of biomarkers used in the analysis (K>=3).
#' @return a list containing (1) sumstats.all: a matrix that contains the merged summary data;
#' (2) traitvec: a vector storing the names of the biomarkers
#' @export
preprocess = function(sumstats,maf.thr,K){
  if (is.null(names(sumstats))){
    names(sumstats) = traitvec
  }
  ###### For the estimation of between-trait correlation:
  # construct sumstats[[1:K]] as the data matrices for biomarkers
  trait = names(sumstats)[1]
  sumstats[[trait]] = sumstats[[trait]] %>%
    mutate(z = beta/se, A1 = toupper(A1), A2 = toupper(A2)) %>%
    select(Chr, Pos, rsid, A1, A2, N, z, P)
  ### delete duplicated rows (with the same rsid and beta, sd, A1, A2 etc)
  for (trait in names(sumstats)[1:K]){
    duplicated.rsid = duplicated(sumstats[[trait]][,'rsid'])
    sumstats[[trait]] = sumstats[[trait]][-which(duplicated.rsid),]
  }
  for (trait in names(sumstats)[2:K]){
    sumstats[[trait]] = sumstats[[trait]] %>%
      mutate(z = beta/se, A1 = toupper(A1), A2 = toupper(A2)) %>%
      select(rsid, A1, A2, N, z, P)
  }

  # NA: Remove SNPs with sample size < 0.67 * (90 percentile)
  # Remove SNPs within the major histocompatibility complex (MHC) region (26Mb~34Mb on chromosome 6)
  # Remove SNPs with MAF (or 1-MAF) < maf.thr
  for (trait in names(sumstats)){
    if ('N' %in% names(sumstats[[trait]])){
      sumstats[[trait]] = sumstats[[trait]] %>% filter(N>0.67*quantile(N,0.9))
    }
    if (("Chr"%in%colnames(sumstats[[trait]]) & ("Pos"%in%colnames(sumstats[[trait]])))){
      sumstats[[trait]] = sumstats[[trait]] %>% filter(!(Chr==6 & Pos>26e6 & Pos<34e6))
    }
    if ("Freq1"%in%colnames(sumstats[[trait]])){
      sumstats[[trait]] = sumstats[[trait]] %>% filter((Freq1>maf.thr) & (Freq1<1-maf.thr))
    }
  }
  ### Trait-specific names
  trait = names(sumstats)[1]
  names(sumstats[[trait]])[4:ncol(sumstats[[trait]])] = paste0(names(sumstats[[trait]])[4:ncol(sumstats[[trait]])],'.',trait)
  for (trait in names(sumstats)[2:K]){
    names(sumstats[[trait]])[2:ncol(sumstats[[trait]])] = paste0(names(sumstats[[trait]])[2:ncol(sumstats[[trait]])],'.',trait)
  }
  # Remove missing data, make A1 represent minor allele, remove rare variants
  # Merge data set: only keep SNPs that are available for all traits.
  sumstats.all = sumstats[[1]]
  for (i in 2:length(sumstats)){
    sumstats.all = sumstats.all %>% inner_join(sumstats[[i]], by = 'rsid')
  }

  # Remove missing data
  sumstats.all = sumstats.all[complete.cases(sumstats.all),]

  # If Freq1 is provided, align the alleles of the first trait such that A1 is the minor allele
  trait = names(sumstats)[1]
  if (paste0("Freq1.",trait) %in% colnames(sumstats.all)){
    tempA1 = sumstats.all[[paste0('A1.',trait)]]
    tempA2 = sumstats.all[[paste0('A2.',trait)]]
    tempz = sumstats.all[[paste0("z.",trait)]]
    # Retrieve indices that needs to be changed
    inds = sumstats.all[[paste0("Freq1.",trait)]] > 0.5
    tempA1[inds] = sumstats.all[[paste0('A2.',trait)]][inds]
    tempA2[inds] = sumstats.all[[paste0('A1.',trait)]][inds]
    tempz[inds] = -tempz[inds]
    sumstats.all = sumstats.all %>% mutate(A1 = tempA1, A2 = tempA2)
    sumstats.all[[paste0("z.",trait)]] = tempz
  } else{
    sumstats.all = sumstats.all %>% rename(A1 = paste0('A1.',trait), A2 = paste0('A2.',trait))
  }

  ###### Remove strand-ambiguous SNPs + Align all alleles to the first trait (A1 is minor allele)
  for (i in 2:length(sumstats)){
    trait = names(sumstats)[i]
    if ('Freq1' %in% names(sumstats[[trait]])){
      inds.keep = (((sumstats.all[['A1']] == sumstats.all[[paste0('A1.',trait)]])&(sumstats.all[['A2']] == sumstats.all[[paste0('A2.',trait)]])&(sumstats.all[[paste0('Freq1.',trait)]] < 0.5)) | ((sumstats.all[['A1']] == sumstats.all[[paste0('A2.',trait)]])&(sumstats.all[['A2']] == sumstats.all[[paste0('A1.',trait)]])&(sumstats.all[[paste0('Freq1.',trait)]]>0.5)))
    }
    if (!'Freq1' %in% names(sumstats[[trait]])){
      inds.keep = (((sumstats.all[['A1']] == sumstats.all[[paste0('A1.',trait)]])&(sumstats.all[['A2']] == sumstats.all[[paste0('A2.',trait)]])) | ((sumstats.all[['A1']] == sumstats.all[[paste0('A2.',trait)]])&(sumstats.all[['A2']] == sumstats.all[[paste0('A1.',trait)]])))
    }
    sumstats.all = sumstats.all[inds.keep,]
    inds = (sumstats.all[[paste0("A1.",trait)]] != sumstats.all[["A1"]])
    sumstats.all[[paste0("z.",trait)]][inds] = -sumstats.all[[paste0("z.",trait)]][inds]
  }
  trait.spec = NULL
  for (trait in names(sumstats)){
    trait.spec = c(trait.spec, paste0(c("N.", "z.", "P."), trait))
  }
  sumstats.all = sumstats.all[,c("rsid", "A1", "A2", trait.spec)]
  return(list(sumstats.all = sumstats.all, traitvec = names(sumstats)))
}
