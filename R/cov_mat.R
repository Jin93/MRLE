#' Estimate between-biomarker genetic correlation matrix.
#'
#' @param sumstats.all a matrix storing the GWAS summary-level association statistics of
#' the biomarkers.
#' @param traitvec a vector storing the names of the biomarkers.
#' @param out.path path to the folder where the output from LD score regression (LDSC) will be saved.
#' @param ldsc.path path to the folder where the required files for LDSC are stored.
#' @param python.path path to the python software.
#' @param ldscore.path path to the ldscore files. Default: file.path(ldsc.path,"eur_w_ld_chr/").
#' @param mergeallele a logical constant for whether or not the alleles should be merged.
#' Default: TRUE.
#' @return a list containing (1) coherit.mat: the estimated genetic correlation matrix;
#' (2) ldscint.mat: the calculated LD score regression intercept matrix.
#' @export
covmat = function(sumstats.all, traitvec, out.path, ldsc.path, python.path = NULL,
                  ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), mergeallele = TRUE){
  # Add anaconda python to path: for ldsc (no need to do this if the path is already set by the system).
  if (!is.null(python.path)){
    Sys.setenv(PATH = paste(python.path, Sys.getenv('PATH'), sep = ':'))
  }
  # Check if out.path exists
  if (!dir.exists(out.path)){
    dir.create(out.path)
  }

  # Further preprocess the data: randomly simulate signs for z.trait1 (LDSC requires median(z) is small).
  inds = sample(nrow(sumstats.all), round(nrow(sumstats.all))/2)
  tempA1 = sumstats.all$A1
  sumstats.all$A1[inds] = sumstats.all$A2[inds]
  sumstats.all$A2[inds] = tempA1[inds]
  for (trait in traitvec){
    sumstats.all[[paste0("z.",trait)]][inds] = -sumstats.all[[paste0("z.",trait)]][inds]
  }

  # Write files for LDSC
  for (trait in traitvec){
    sumstats.all %>%
      rename(N = paste0("N.",trait), z = paste0("z.",trait), P = paste0("P.",trait)) %>%
      select(rsid, A1, A2, N, z, P) %>%
      write_delim(file.path(out.path, paste0(trait,'_for_ldsc.txt')), delim = '\t')
  }

  # Call ldsc to estimate heritability, coheritability
  # Change summary level data to .sumstats format
  for (trait in traitvec){
    munge.sumstats.code = paste(python.path,file.path(ldsc.path,"munge_sumstats.py"),
                                "--sumstats",
                                file.path(out.path, paste0(trait,'_for_ldsc.txt')),
                                "--out", file.path(out.path, paste0(trait,"_ldsc_format")))
    if (mergeallele==TRUE){
      munge.sumstats.code = paste(munge.sumstats.code, "--merge-alleles", file.path(ldsc.path,"w_hm3.snplist"))
    }
    system(munge.sumstats.code)
    print(paste("ldsc_sumstasts", trait))
  }

  # Fit LD score regression
  for (i in 1:(length(traitvec)-1)){
    traits = traitvec[i:(length(traitvec))]
    sumstats.files = paste(file.path(out.path, paste0(traits,"_ldsc_format.sumstats.gz")), collapse = ",")
    ldsc.code = paste(python.path, file.path(ldsc.path,"ldsc.py"),
                      "--rg", sumstats.files,
                      "--ref-ld-chr", ldscore.path,
                      "--w-ld-chr", ldscore.path,
                      "--out", file.path(out.path, paste0(paste(traits, collapse = "_"),"_ldsc_results")))
    system(ldsc.code)
    print(paste("ldsc", traits[1]))
  }

  # Retrieve heritability-coheritability matrix from the log files
  coherit.mat = matrix(NA, nrow = length(traitvec), ncol = length(traitvec))
  ldscint.mat = matrix(NA, nrow = length(traitvec), ncol = length(traitvec)) # Stores LDSC intercepts

  for (i in 1:(length(traitvec)-1)){
    traits = traitvec[i:length(traitvec)]
    logfile = readLines(file.path(out.path, paste0(paste(traits, collapse = "_"),"_ldsc_results.log")))
    # Get heritability information
    if (i==1){
      for (j in 1:length(traitvec)){
        if (j == 1){  # When j==1 the format is slightly different
          ind = which(logfile == "Heritability of phenotype 1") # Two lines below the title is the heritability
        } else{
          ind = which(logfile == paste0("Heritability of phenotype ", j, "/", length(traitvec)))
        }
        coherit.mat[j,j] = as.numeric(strsplit(logfile[ind+2], split = " ")[[1]][5])
        ldscint.mat[j,j] = as.numeric(strsplit(logfile[ind+5], split = " ")[[1]][2])
      }
    }

    # Get coheritability information
    gencov.inds = which(logfile == "Genetic Covariance")
    gencovs = sapply(strsplit(logfile[gencov.inds+2], split = " "), function(x) x[5])
    coherit.mat[i,(i+1):length(traitvec)] = as.numeric(gencovs)
    coherit.mat[(i+1):length(traitvec),i] = as.numeric(gencovs)

    ldscints = sapply(strsplit(logfile[gencov.inds+4], split = " "), function(x) x[2])
    ldscint.mat[i,(i+1):length(traitvec)] = as.numeric(ldscints)
    ldscint.mat[(i+1):length(traitvec),i] = as.numeric(ldscints)
  }

  return(list(coherit.mat = coherit.mat, ldscint.mat = ldscint.mat))
}


#' Provide information on between-biomarker correlation matrix.
#'
#' @param sumstats a list of matrices storing the GWAS summary-level association statistics of
#' the biomarkers.
#' @param out.path path to the folder where the output from LD score regression (LDSC) will be saved.
#' @param ldsc.path path to the folder where the required files for LDSC are stored.
#' @param python.path path to the python software.
#' @param ldscore.path path to the ldscore files. Default: file.path(ldsc.path,"eur_w_ld_chr/").
#' @param mergeallele a logical constant for whether or not the alleles should be merged.
#' Default: TRUE.
#' @return a list containing (1) coherit.mat: the estimated genetic correlation matrix;
#' (2) ldscint.mat: the calculated LD score regression intercept matrix.
#' @export
covinfo = function(sumstats, out.path, ldsc.path, python.path = NULL,
                   ldscore.path = file.path(ldsc.path,"eur_w_ld_chr/"), maf.thr, mergeallele = TRUE, K){
  temp = preprocess(sumstats, maf.thr, K = K)
  sumstats.all = temp[[1]]
  traitvec = temp[[2]]
  rm(temp)

  temp = covmat(sumstats.all, traitvec, out.path, ldsc.path = ldsc.path, python.path = python.path, ldscore.path = ldscore.path, mergeallele = mergeallele)
  coherit.mat = temp[[1]]
  ldscint.mat = temp[[2]]

  return(list(coherit.mat = coherit.mat, ldscint.mat = ldscint.mat))
}

