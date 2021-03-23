################################################################################
#' Objective function
OF = function(tet,X,Cov.mat,W){
  # x = sbetay,ssigmay,sbetaBk,ssigmaBk
  p = length(tet)
  K = X[[1]]
  sbetay = X[[2]]
  ssigmay = X[[3]]
  sbetaBk = X[[4]]
  ssigmaBk = X[[5]]
  M=length(sbetay)
  # tet: theta (1), thetak (2:(K+1)), h2gammak ((K+2):(2*K+1)), h2x (2*(K+1))
  m_y = (tet[1])^2 - mean((sbetay)^2 - (ssigmay)^2)
  m_yb = m_b = m_bjk = numeric()
  l=1
  for (k in 1:K)
  {
    m_yb[k] = tet[1] * tet[k+1] - (mean(sbetay * sbetaBk[[k]]) - mean(sapply(1:M,function(x){Cov.mat[[x]][1,k+1]})) )#cov.mat[1,k+1])
    m_b[k] = tet[K+1+k] + (tet[k+1])^2 - mean((sbetaBk[[k]])^2 - (ssigmaBk[[k]])^2)
    if (k<=(K-1))
    {
      for (k2 in (k+1):K)
      {
        m_bjk[l] = tet[k+1] * tet[k2+1] - (mean(sbetaBk[[k]] * sbetaBk[[k2]]) - mean(sapply(1:M,function(x){Cov.mat[[x]][k+1,k2+1]}) ))#cov.mat[k+1,k2+1])
        l=l+1
      }
    }
  }
  f <- c(m_y,m_yb,m_b,m_bjk)
  f = f[-1]
  t(f)%*% W %*% f
}
################################################################################
#' estimating equations
g <- function(tet,x,cov.mat){
  # x = sbetay,ssigmay,sbetaBk,ssigmaBk
  p = length(tet)
  sbetay = x[1]
  ssigmay = x[2]
  sbetaBk = ssigmaBk = list()
  for (i in 1:K){
    sbetaBk[[i]] = x[2+i]
    ssigmaBk[[i]] = x[2+K+i]
  }
  # tet: theta (1), thetak (2:(K+1)), h2gammak ((K+2):(2*K+1)), h2x (2*(K+1))
  m_y = (tet[1])^2 - ((sbetay)^2 - (ssigmay)^2)
  m_yb = m_b = m_bjk = numeric()
  l=1
  for (k in 1:K){
    m_yb[k] = tet[1] * tet[k+1] - (sbetay * sbetaBk[[k]] - cov.mat[1,k+1])
    m_b[k] = tet[K+1+k] + (tet[k+1])^2 - ((sbetaBk[[k]])^2 - (ssigmaBk[[k]])^2)
    if (k<=(K-1)){
      for (k2 in (k+1):K){
        m_bjk[l] = tet[k+1] * tet[k2+1] - (sbetaBk[[k]] * sbetaBk[[k2]] - cov.mat[k+1,k2+1])
        l=l+1
      }
    }
  }
  f <- c(m_y,m_yb,m_b,m_bjk)
  f = f[-1]
  return(f)
}
################################################################################
#' derivative matrix
Dg <- function(tet){
  # x = sbetay,ssigmay,sbetaBk,ssigmaBk
  p = length(tet)
  #### tet: theta (1), thetak (2:(K+1)), h2gammak ((K+2):(2*K+1)), h2x (2*(K+1))
  dg <- matrix(0,(K+1)*(K+2)/2,p)
  dg[1,1] = 2 * tet[1]
  for (j in 2:(K+1)){
    k = j-1 # index for Bk in Cov(Y,Bk)
    dg[j,1] = tet[k+1]
    dg[j,k+1] = tet[1]
  }
  for (j in (K+2):(2*K+1)){
    k = j - (K+1) # index for Bk in Cov(Y,Bk)
    dg[j,k+1] = 2 * tet[k+1]
    dg[j,K+1+k] = 1
  }
  l=1
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      j = 2*K + 1 + l # index of the row in dg
      ###### (k+1)-th column: thetak2*h2x
      ###### (k2+1)-th column: thetak*h2x
      ###### 2*(K+1)-th column: thetak*thetak2
      dg[j,k+1] = tet[k2+1]
      dg[j,k2+1] = tet[k+1]
      l=l+1
    }
  }
  dg = dg[-1,]
  return(dg)
}
################################################################################
#' derivative matrix
gn <- function(tet,x,Cov.mat){
  # x = K,sbetay,ssigmay,sbetaBk,ssigmaBk
  p = length(tet)
  K = x[[1]]
  sbetay = x[[2]]
  ssigmay = x[[3]]
  sbetaBk = x[[4]]
  ssigmaBk = x[[5]]
  M=length(sbetay)
  #### tet: theta (1), thetak (2:(K+1)), h2gammak ((K+2):(2*K+1)), h2x (2*(K+1))
  m_y = (tet[1])^2 - mean((sbetay)^2 - (ssigmay)^2)
  m_yb = m_b = m_bjk = numeric()
  l=1
  for (k in 1:K){
    mean(sapply(1:M,function(x){Cov.mat[[x]][1,k+1]}))
    m_yb[k] = tet[1] * tet[k+1] - (mean(sbetay * sbetaBk[[k]]) - mean(sapply(1:M,function(x){Cov.mat[[x]][1,k+1]})))#cor.mat[1,k+1]*mean(ssigmay*ssigmaBk[[k]]))
    m_b[k] = tet[K+1+k] + (tet[k+1])^2 - mean((sbetaBk[[k]])^2 - (ssigmaBk[[k]])^2)
    if (k<=(K-1)){
      for (k2 in (k+1):K){
        m_bjk[l] = tet[k+1] * tet[k2+1] - (mean(sbetaBk[[k]] * sbetaBk[[k2]]) - mean(sapply(1:M,function(x){Cov.mat[[x]][k+1,k2+1]})))#cor.mat[k+1,k2+1]*mean(ssigmaBk[[k]]*ssigmaBk[[k2]]))
        l=l+1
      }
    }
  }
  f <- c(m_y,m_yb,m_b,m_bjk)
  f = f[-1]
  f
}
################################################################################
#' second derivatives
D2Gn = function(eta){
  d2Gn = list()
  p= 2*K+1
  d = (K+1)*(K+2)/2
  d2Gn[[1]] = matrix(0,d,p)
  d2Gn[[1]][1,1] = 2
  for (k in 1:K){
    d2Gn[[1]][k+1,k+1] = 1
  }
  #### (k+1)'s:
  for (k in 1:K){
    d2Gn[[k+1]] = matrix(0,d,p)
    d2Gn[[k+1]][k+1,1] = 1
    d2Gn[[k+1]][K+1+k,k+1] = 2
    if (k > 1){
      for (j in 1){
        d2Gn[[k+1]][2*K+1+(k-j),1+j] = 1
      }
      if (2<=(k-1)){
        for (j in 2:(k-1)){
          d2Gn[[k+1]][2*K+1+sum(sapply(1:(j-1),FUN=function(x){K-x}))+(k-j),1+j] = 1
        }
      }
    }
    if (k < K){
      if (k == 1){
        for (j in (k+1):K){
          d2Gn[[k+1]][2*K+1+(j-k),1+j] = 1
        }
      }
      if (k > 1){
        for (j in (k+1):K){
          d2Gn[[k+1]][2*K+1+sum(sapply(1:(k-1),FUN=function(x){K-x}))+(j-k),1+j] = 1
        }
      }
    }
  }
  #### (K+k+1)'s:
  for (k in 1:K){
    d2Gn[[K+k+1]] = matrix(0,d,p)
  }
  for (ii in 1:length(d2Gn)){
    d2Gn[[ii]] = d2Gn[[ii]][-1,]
  }
  d2Gn
}
################################################################################
#' covariance matrix of the estimating equations
cal.Omega = function(eta, sigmasqBk, sigmasqY, cov.mat){
  # M0: the # of IVs selected
  d = (K+1)*(K+2)/2
  p = 2*(K+1)
  theta = eta[1]
  thetak = eta[2:(K+1)]
  h2gammak = eta[(K+2):(2*K+1)]
  h2x = 1
  sigmasqx = sigmasqk = sigmasqy = 1

  Omega = matrix(NA,d,d)
  ######### column, row 1
  colnames(Omega) = 1:d
  Omega[1,1] = 2*(theta^4*h2x^2 + sigmasqY^2 + 2*sigmasqY*theta^2*h2x^2)
  for (i in 2:(K+1)){
    #i = k + 1: E(g_{k+1}g_1)
    k = i-1
    Omega[i,1] = Omega[1,i] = 2*(theta^3)*thetak[k]*(h2x)^2 + 2*theta*thetak[k]*h2x*sigmasqY + 2*(theta^2)*h2x*cov.mat[1,k+1]
    #2*theta*thetak[k]*h2x*(theta^2*h2x + sigmasqY + cov.mat[1,k+1])
  }
  ### (K+k+1,1)
  for (k in 1:K){
    Omega[K+k+1,1] = Omega[1,K+k+1] = 2*theta^2*thetak[k]^2*h2x^2 + 4*theta*thetak[k]*h2x*cov.mat[1,k+1]
  }
  ### ((k,k2),1)
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      i = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (k2-k)
      # E(g_{(k,k2)}g_1)
      Omega[i,1] = Omega[1,i] = 2*theta^2*thetak[k]*thetak[k2]*h2x^2 + 2*theta*thetak[k2]*h2x*cov.mat[1,k+1] + 2*theta*thetak[k]*h2x*cov.mat[1,k2+1]
    }
  }
  ######### column, row (k+1)
  ###（k+1,k+1)
  for (k in 1:K){
    Omega[k+1,k+1] = theta^2*h2x*(2*(thetak[k]^2)*h2x + h2gammak[k]) + sigmasqBk[k]*(theta^2*h2x) + sigmasqY*(thetak[k]^2*h2x + h2gammak[k]) + 2*(theta*thetak[k]*h2x)*cov.mat[1,k+1] + sigmasqBk[k] * sigmasqY
  }
  ### (k+1,l+1)
  for (k in 1:(K-1)){
    for (l in (k+1):K){
      Omega[k+1,l+1] = Omega[l+1,k+1] = 2*thetak[k]*thetak[l]*theta^2*(h2x)^2 + thetak[k]*thetak[l]*h2x * sigmasqY + theta*thetak[k]*h2x *cov.mat[1,l+1] + theta*thetak[l]*h2x *cov.mat[1,k+1] + theta^2*h2x *cov.mat[k+1,l+1]
    }
  }
  ### (K+k+1,k+1)
  for (k in 1:K){
    Omega[K+k+1,k+1] = Omega[k+1,K+k+1] = 2*theta*thetak[k]*h2x * (thetak[k]^2*h2x + h2gammak[k]) + 2*theta*thetak[k]*h2x*sigmasqBk[k]
  }
  ### (K+k+1,l+1)
  for (k in 1:(K-1)){
    for (l in (k+1):K){
      Omega[K+k+1,l+1] = Omega[l+1,K+k+1] = 2*theta*thetak[k]^2*thetak[l]*(h2x^2) + 2*thetak[k]*thetak[l]*h2x*cov.mat[1,k+1] + 2*theta*thetak[k]*h2x*cov.mat[k+1,l+1]
      Omega[K+l+1,k+1] = Omega[k+1,K+l+1] = 2*theta*thetak[l]^2*thetak[k]*(h2x^2) + 2*thetak[l]*thetak[k]*h2x*cov.mat[1,l+1] + 2*theta*thetak[l]*h2x*cov.mat[l+1,k+1]
    }
  }
  ### ((k,k2),k+1):
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      i = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (k2-k)
      Omega[i,k+1] = Omega[k+1,i] = theta*thetak[k2]*h2x* (2*thetak[k]^2*h2x + h2gammak[k]) + thetak[k2]*theta*h2x *sigmasqBk[k] + thetak[k]*thetak[k2]*h2x*cov.mat[1,k+1] + ((thetak[k]^2)*h2x + h2gammak[k])*cov.mat[1,k2+1] + (theta*thetak[k]*h2x)*cov.mat[k+1,k2+1]
      Omega[i,k2+1] = Omega[k2+1,i] = theta*thetak[k]*h2x* (2*thetak[k2]^2*h2x + h2gammak[k2]) + (thetak[k]*theta*h2x)*sigmasqBk[k2] + (thetak[k2]*thetak[k]*h2x)*cov.mat[1,k2+1] + ((thetak[k2]^2)*h2x + h2gammak[k2])*cov.mat[1,k+1] + (theta*thetak[k2]*h2x)*cov.mat[k2+1,k+1]
    }
  }
  ### ((k,k2),l+1):
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      for (l in c(1:K)[-c(k,k2)]){
        i = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (k2-k)
        Omega[i,l+1] = Omega[l+1,i] = 2*theta*thetak[k]*thetak[k2]*thetak[l]*h2x^2 + h2x * ( (thetak[k2]*thetak[l])*cov.mat[1,k+1] + (thetak[k2]*theta)*cov.mat[k+1,l+1] + (thetak[k]*thetak[l])*cov.mat[1,k2+1] + (theta*thetak[k])*cov.mat[k2+1,l+1])
      }
    }
  }
  ########## column, row (K+k+1):
  ###（K+k+1,K+k+1)
  for (k in 1:K){
    Omega[K+k+1,K+k+1] = 2*((thetak[k]^2*h2x + h2gammak[k])^2 + 2*sigmasqBk[k]*(thetak[k]^2*h2x + h2gammak[k] )) + 2*(sigmasqBk[k]^2)  ####
  }
  ### (K+k+1,K+l+1)
  for (k in 1:K){
    for (l in c(1:K)[-k]){
      Omega[K+k+1,K+l+1] = Omega[K+l+1,K+k+1] = 2*thetak[k]^2*thetak[l]^2*(h2x^2) + 4*cov.mat[k+1,l+1]*(thetak[k]*thetak[l]*h2x )
    }
  }
  ### ((k,k2),K+k+1):
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      i = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (k2-k)
      Omega[i,K+k+1] = Omega[K+k+1,i] = 2*thetak[k]*thetak[k2]*h2x * (thetak[k]^2*h2x + h2gammak[k]) + 2*(thetak[k]*thetak[k2]*h2x )*sigmasqBk[k] + 2*((thetak[k])^2*h2x+h2gammak[k])*cov.mat[k+1,k2+1]
      Omega[i,K+k2+1] = Omega[K+k2+1,i] = 2*thetak[k2]*thetak[k]*h2x * (thetak[k2]^2*h2x + h2gammak[k2]) + 2*(thetak[k2]*thetak[k]*h2x )*sigmasqBk[k2] + 2*((thetak[k2])^2*h2x+h2gammak[k2])*cov.mat[k2+1,k+1]
    }
  }
  ### ((k,k2),K+l+1):
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      for (l in c(1:K)[-c(k,k2)]){
        i = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (k2-k)
        Omega[i,K+l+1] = Omega[K+l+1,i] = 2*thetak[k]*thetak[k2]*thetak[l]^2*h2x^2 + 2*(thetak[k2]*thetak[l]*h2x )*cov.mat[k+1,l+1] + 2*(thetak[k]*thetak[l]*h2x)*cov.mat[k2+1,l+1]
      }
    }
  }
  ########## column, row ((k,k2)):
  ###((k,k2),(k,k2))
  for (k in 1:(K-1)){
    for (k2 in (k+1):K){
      i = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (k2-k)
      Omega[i,i] = 2*thetak[k]^2*thetak[k2]^2*(h2x^2) + ((thetak[k2]^2)*h2x +h2gammak[k2])*sigmasqBk[k] + ((thetak[k]^2)*h2x +h2gammak[k])*sigmasqBk[k2] + 2*(thetak[k]*thetak[k2]*h2x)*cov.mat[k+1,k2+1]  + ((thetak[k]^2)*h2gammak[k2] + (thetak[k2]^2)*h2gammak[k])*h2x + h2gammak[k]*h2gammak[k2] + sigmasqBk[k]*sigmasqBk[k2]####
    }
  }
  ### ((k,l1),(k,l2))
  if (K >= 3){
    for (k in 1:(K-2)){
      for (l1 in (k+1):K){
        for (l2 in (k+1):K){
          if (l1!=l2){
            i1 = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (l1-k)
            i2 = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (l2-k)
            Omega[i1,i2] = Omega[i2,i1] = thetak[l1] * thetak[l2] * h2x * (2*h2x*(thetak[k]^2) + h2gammak[k]) + (thetak[l1]*thetak[l2]*h2x )*sigmasqBk[k] + (thetak[l1]*thetak[k]*h2x )*cov.mat[k+1,l2+1] + (thetak[k]*thetak[l2]*h2x )*cov.mat[l1+1,k+1] + ((thetak[k]^2)*h2x + h2gammak[k] )*cov.mat[l1+1,l2+1]
          }
        }
      }
    }
  }

  ### ((k,l1),(l1,l2))
  if (K >= 3){
    for (k in 1:(K-2)){
      for (l1 in (k+1):(K-1)){
        for (l2 in (l1+1):K){
          if (l1!=l2){
            i1 = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (l1-k)
            i2 = 2*K+1 + ifelse(l1>1,sum(sapply(1:(l1-1),FUN=function(x){K-x})),0) + (l2-l1)
            Omega[i1,i2] = Omega[i2,i1] = thetak[k]*thetak[l2]*h2x*(2*h2x*(thetak[l1]^2) + h2gammak[l1]) + (thetak[k]*thetak[l2]*h2x )*sigmasqBk[l1] + (thetak[l1]*thetak[l2]*h2x )*cov.mat[k+1,l1+1] + ((thetak[l1]^2)*h2x + h2gammak[l1] )*cov.mat[k+1,l2+1] + (thetak[k]*thetak[l1]*h2x )*cov.mat[l1+1,l2+1]
          }
        }
      }
    }
  }
  ### ((k,l1),(l2,k)): same as ((k,l1),(l1,l2))
  ### ((k,l),(k2,l)):
  if (K >= 3){
    for (l in 3:K){
      for (k in 1:(l-1)){
        for (k2 in 1:(l-1)){
          if (k!=k2){
            i1 = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (l-k)
            i2 = 2*K+1 + ifelse(k2>1,sum(sapply(1:(k2-1),FUN=function(x){K-x})),0) + (l-k2)
            Omega[i1,i2] = Omega[i2,i1] = thetak[k]*thetak[k2]*h2x*(2*h2x*(thetak[l]^2) + h2gammak[l]) + (thetak[k]*thetak[k2]*h2x )*sigmasqBk[l] + (thetak[l]*thetak[k2]*h2x )*cov.mat[k+1,l+1] + ((thetak[l]^2)*h2x + h2gammak[l] )*cov.mat[k+1,k2+1] + (thetak[k]*thetak[l]*h2x )*cov.mat[l+1,k2+1]
          }
        }
      }
    }
  }

  ### ((k,l1),(k2,l2))
  for (k in 1:(K-1)){
    for (k2 in 1:(K-1)){
      for (l1 in (1+k):K){
        for (l2 in (1+k2):K){
          if (length(unique(c(k,k2,l1,l2)))==4){
            i1 = 2*K+1 + ifelse(k>1,sum(sapply(1:(k-1),FUN=function(x){K-x})),0) + (l1-k)
            i2 = 2*K+1 + ifelse(k2>1,sum(sapply(1:(k2-1),FUN=function(x){K-x})),0) + (l2-k2)
            Omega[i1,i2] = Omega[i2,i1] = 2*thetak[k] * thetak[k2] *thetak[l1] * thetak[l2] * h2x^2 + (thetak[l1]*thetak[l2]*h2x)*cov.mat[k+1,k2+1] + (thetak[l1]*thetak[k2]*h2x)*cov.mat[k+1,l2+1] + (thetak[k]*thetak[l2]*h2x)*cov.mat[l1+1,k2+1] + (thetak[k]*thetak[k2]*h2x)*cov.mat[l1+1,l2+1]
          }
        }
      }
    }
  }
  Omega = Omega[-1,-1]
  Omega
}
################################################################################
#' estimate asymptotic covariance of the GMM estimator
asymcov.gmm = function(eta, sigmasqBk, sigmasqY, cov.mat, stage2.conduct=T){
  d = (K+1)*(K+2)/2 - 1
  p = 2*K+1
  theta = eta[1]
  thetak = eta[2:(K+1)]
  h2gammak = eta[(K+2):(2*K+1)]
  sigmasqx = sigmasqk = sigmasqy = 1
  Omega = cal.Omega(eta, sigmasqBk, sigmasqY, cov.mat)
  dGn = Dg(eta)
  if (stage2.conduct == T){
    V = t(dGn) %*% solve(Omega,tol=1e-80) %*% dGn
    V = solve(V, tol=1e-80)
  }
  if (stage2.conduct == F){
    gwginv = solve(t(dGn) %*% (dGn),tol=1e-80)
    V = gwginv %*% t(dGn) %*% Omega %*% dGn %*% gwginv
  }
  V
}

################################################################################
#' IVW estimator
ivw = function(sbetak,ssigmak,sbetayk,ssigmayk,i){
  mr.obj = mr_input(bx = sbetak, bxse = ssigmak,
                    by = sbetayk, byse = ssigmayk)
  res.ivw = mr_ivw(mr.obj,model='fixed')
  est.ivw = attr(res.ivw,"Estimate")
  se.ivw = attr(res.ivw,"StdError")
  p.ivw = attr(res.ivw,"Pvalue")#pnorm(abs(est.ivw/se.ivw),0,1,lower.tail=F)
  rej.ivw = ifelse(attr(res.ivw,"Pvalue")<alpha0,1,0)
  ###### consider the sign of thetak:
  est.ivw = est.ivw*thetak.sign[i]
  return(list(est=est.ivw,se=se.ivw,p=p.ivw,rej=rej.ivw))
}

################################################################################
#' Test based on IVW estimator
ivwout = function(sbetaBk,ssigmaBk,sbetay,ssigmay,multiplebiom){
  if (multiplebiom == T){
    ests = ses = ps = rejs = list()
    for (i in 1:K){
      if (length(ivindk[[i]])>0){
        sbetak = sbetaBk[[i]][ivindk[[i]]]
        ssigmak = ssigmaBk[[i]][ivindk[[i]]]
        sbetayk = sbetay[ivindk[[i]]]
        ssigmayk = ssigmay[ivindk[[i]]]
        ivwres = ivw(sbetak,ssigmak,sbetayk,ssigmayk,i)
        ests[[i]] = ivwres$est
        ses[[i]] = ivwres$se
        ps[[i]] = ivwres$p
        rejs[[i]] = ivwres$rej
      }
    }
  }
  if (multiplebiom == F){
    i=1
    sbetak = sbetaBk[[i]][ivindk[[i]]]
    ssigmak = ssigmaBk[[i]][ivindk[[i]]]
    sbetayk = sbetay[ivindk[[i]]]
    ssigmayk = ssigmay[ivindk[[i]]]
    ivwres = ivw(sbetak,ssigmak,sbetayk,ssigmayk,i)
    ests = ivwres$est
    ses = ivwres$se
    ps = ivwres$p
    rejs = ivwres$rej
  }
  return(list(est=ests,se=ses,p=ps,rej=rejs))
}

