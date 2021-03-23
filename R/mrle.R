#' Mendelian Randomization analysis using multiple Biomarkers of an underlying common exposure
#'
#' This function conducts Mendelian Randomization (MR) analysis for a latent exposure (MRLE).
#' The method uses summary-level association statistics from genome-wide association study
#' of the outcome and K>=3 observable traits/biomarkers that are co-regulated by the latent exposure
#' on a set of SNPs (instrumental variables) that are associated with
#' at least two of the traits/biomarkers at specified thresholds of significance.
#'
#' @param sumtable the M by (K+1)*3+6 input matrix containing
#' @param thetak_sign a vector of length K indicating the directions of the effect of the latent exposure on
#' the traits/biomarkers
#' @param Cov.mat the (K+1) by (K+1) estimated covariance matrix of the GWAS summary statistics.
#' The first variable corresponds to the outcome; the 2nd to the (K+1)-th variables correspond to
#' the k-th observable traits.
#' @param alpha0 significance level of the hypothesis test.
#' @return A matrix of the infile
#' @export

mrle = function(sumtable, thetak_sign, Cov.mat, alpha0 = 5e-2){
  n.rep=200; C = 1
  K = (ncol(sumtable)-6)/3-1
  traitvec = sapply(1:K,function(x){strsplit(colnames(sumtable)[5+3*x],'\\.')[[1]][2]})
  outcome = strsplit(colnames(sumtable)[5+3*(K+1)],'\\.')[[1]][2]
  ##### GWAS summary statistics for Y (outcome):
  sbetay = sumtable[[paste0("beta.",outcome)]]
  ssigmay = sumtable[[paste0("se.",outcome)]]
  ##### GWAS summary statistics for biomarkers:
  sbetaBk = ssigmaBk = list()
  for (i in 1:K){
    sbetaBk[[i]] = sumtable[[paste0("beta.",traitvec[i])]]
    ssigmaBk[[i]] = sumtable[[paste0("se.",traitvec[i])]]
  }
  ###### Select the SNPs that reach genome-wide significance in the study associated with Bks:
  M = nrow(sumtable)
  X = list(K=K,sbetay=sbetay,ssigmay=ssigmay,sbetaBk=sbetaBk,ssigmaBk=ssigmaBk)
  x = cbind(sbetay,ssigmay,sapply(1:K,FUN=function(x){sbetaBk[[x]]}),sapply(1:K,FUN=function(x){ssigmaBk[[x]]}))
  ####### 1. IVW Estimator:
  #### assume that we know the sign of thetak:
  ivw.results = ivwout(sbetaBk,ssigmaBk,sbetay,ssigmay,multiplebiom=T)

  set.seed(2021)
  #### assume that we know the sign of thetak:
  z0 = qnorm(p=alpha0/2,lower.tail=FALSE)
  d = K*(K+3)/2
  p = 2*K+1
  # number of constraits met
  p.constraint = (K<5) * (p-2) + (K==5) * (p-4) + (K==6)*(p-4) + (K==8) * (p-7)
  eta.est = list()
  T.F = numeric(); T.F.max = 1e5
  diff.thr1 = 1e-30; diff.thr2 = 1e-30
  iter.max = 200
  n.nonzero = 0
  n.try = 1
  est.GMM = rej.GMM = 0
  se.GMM = NA
  rep = 1
  while ((n.nonzero <2) | (is.na(se.GMM)) | (is.na(est.GMM))){
    thetainit = runif(1,-sqrt(1),sqrt(1))
    thetakinit = thetak_sign*runif(K,0,sqrt(0.7))
    h2xinit = abs((mean(sbetay^2)-mean(ssigmay^2))/((runif(1,-sqrt(1),sqrt(1)))^2))
    h2kinit = sapply(1:K,FUN=function(x){mean(sbetaBk[[x]]^2)-mean(ssigmaBk[[x]]^2) - h2xinit * ((runif(1,0,sqrt(0.7)))^2)})
    init = c(thetainit*sqrt(h2xinit),thetakinit*sqrt(h2xinit),h2kinit)
    names(init) = c("theta",paste0("theta",1:K),paste0("hsq",1:K))
    etanew = init
    eta = rep(0,p)
    iter = 1
    Jn = diag(1,p)
    VGinv = diag(1,d)
    ###### GMM step 1: obtain eta1 by setting weighting matrix C = I:
    while ((sqrt(sum((eta-etanew)^2)) > diff.thr1)&(iter<=iter.max)){
      eta = etanew
      Gn = gn(eta,X,Cov.mat)
      dGn = Dg(eta)
      Tn = t(dGn) %*% Gn
      d2Gn = D2Gn(eta)
      Jn1 = t(dGn) %*% dGn
      Jn2 = sapply(1:p,FUN=function(x){t(Gn) %*% d2Gn[[x]]})
      Jn = Jn1 + Jn2
      lambda = C * sqrt(sum(Tn^2))
      epsilon = diag(lambda,p)
      Jn = Jn + epsilon
      etanew = eta - solve(Jn,tol=1e-80) %*% Tn
      iter = iter + 1
      t.f = OF(etanew,X,Cov.mat,VGinv)
    }
    eta1 = etanew
    constraints = 1
    for (k in 1:K){
      if (!is.na(sum(thetak_sign))){
        constraints = constraints + ((eta1[k+1]*thetak_sign[k]) > 1e-12) * 1
      }
    }
    for (k in 1:K){
      constraints = constraints + ((eta1[k+K+1]<1/M) & (eta1[k+K+1] > 0))*1
    }
    eta.est[[rep]] = eta1
    T.F[rep] = t.f
    VGinv = cal.Omega(eta1, sigmasqBk = sapply(1:K,FUN=function(x){mean((ssigmaBk[[x]])^2)}),
                      sigmasqY = mean(ssigmay^2), cov.mat=c.mat)
    VGinv = solve(VGinv,tol=1e-80)
    if (constraints < p.constraint){
      eta.est[[rep]] = rep(0,p)
      T.F[rep] = T.F.max
    }
    ###### GMM step 2: obtain the final estimate by setting weighting matrix to
    # the one estimated based on the estimate from step 1.
    if (constraints >=(p.constraint-2)){
      etanew = eta1
      eta = rep(0,p)
      iter = 1
      t.f=numeric()
      Jn = diag(1,p)
      while((sqrt(sum((eta-etanew)^2)) > diff.thr2)&(iter<=iter.max))
      {
        eta = etanew
        VGinv = cal.Omega(eta, sigmasqBk = sapply(1:K,FUN=function(x){mean((ssigmaBk[[x]])^2)}),
                          sigmasqY = mean(ssigmay^2), cov.mat=c.mat)
        VGinv = solve(VGinv,tol=1e-80)
        Gn = gn(eta,X,Cov.mat)
        dGn = Dg(eta)
        Tn = t(dGn) %*% VGinv %*% Gn
        d2Gn = D2Gn(eta) # a list
        Jn1 = t(dGn) %*% VGinv %*% dGn
        Jn2 = sapply(1:p,FUN=function(x){t(Gn) %*% VGinv %*% d2Gn[[x]]})
        Jn = Jn1 + Jn2
        lambda = C * sqrt(sum(Tn^2))
        epsilon = diag(lambda,p)
        Jn = Jn + epsilon
        etanew = eta - solve(Jn,tol=1e-80) %*% Tn
        t.f = OF(etanew,X,Cov.mat,VGinv)
        iter = iter + 1
      }
      eta.est[[rep]] = etanew
      T.F[rep] = t.f
      constraints = 1
      for (k in 1:K){
        if (!is.na(sum(thetak_sign))){
          constraints = constraints + ((eta.est[[rep]][k+1]*thetak_sign[k]) > 1e-12) * 1
        }
      }
      for (k in 1:K){
        constraints = constraints + ((eta.est[[rep]][k+K+1] < 3/M) & (eta.est[[rep]][k+K+1]) > 0)*1
      }
      if ((constraints <p.constraint)|T.F[rep]<0){#|(sum(eta.est[[rep]][(K+2):(K*2+1)]<0)>0)
        eta.est[[rep]] = rep(0,p)
        T.F[rep] = T.F.max
      }
    }
    if (sum(eta.est[[rep]]) != 0) n.nonzero = n.nonzero + 1
    if (rep %% 30 == 0){
      etahat0 = eta.est[[which.min(abs(T.F))]]
      n.try = n.try + 1
      if (sum(etahat0) == 0) etahat0 = rep(NA,p)
      etahat = etahat0
      names(etahat) = c("theta",paste0("theta",1:K),paste0("hsq",1:K))
      est.GMM = se.GMM = p.GMM = rej.GMM = 0
      est.GMM = etahat[1]
      if (!is.na(sum(etahat))){
        asym.cov = asymcov.gmm(etahat, sigmasqBk = sapply(1:K,FUN=function(x){mean((ssigmaBk[[x]])^2)}),
                                       sigmasqY = mean(ssigmay^2), c.mat, stage2.conduct = T)
        se.GMM = sqrt(asym.cov[1,1]/M)
        p.GMM = 2*pnorm(abs(est.GMM/se.GMM),0,1,lower.tail = F)
        rej.GMM = ifelse(abs(est.GMM/se.GMM)>z0,1,0)
      }
    }
    rep = rep + 1
  }
  etahat = etahat0
  names(etahat) = c("theta",paste0("theta",1:K),paste0("hsq",1:K))
  est.GMM = etahat[1]
  if (!is.na(sum(etahat))){
    asym.cov = asymcov.gmm(etahat, sigmasqBk = sapply(1:K,FUN=function(x){median((ssigmaBk[[x]])^2)}),
                           sigmasqY = mean(ssigmay^2), c.mat, stage2.conduct = T)
    se.GMM = sqrt(asym.cov[1,1]/M)
    p.GMM = 2*pnorm(abs(est.GMM/se.GMM),0,1,lower.tail = F)
    rej.GMM = ifelse(abs(est.GMM/se.GMM)>z0,1,0)
  }
  # results = list(P.MRLE = as.numeric(signif(p.GMM,3)),
  #                Rejection.MRLE = ifelse(as.numeric(rej.GMM)==0, 'No', 'Yes'),
  #                Direction.mu = sign(as.numeric(etahat[1])))

  results = rbind(matrix(c(unlist(ivw.results[[3]]), unlist(ivw.results[[4]]), sign(unlist(ivw.results[[1]]))),nrow=K,byrow = F),
                       matrix(c(p.GMM,rej.GMM,sign(est.GMM)),1,3))
  rownames(results) = c(paste0('IVW.',traitvec),'MRLE')
  colnames(results) = c('P-value','Reject','Direction')
  results = as.data.frame(results)
  results
}

