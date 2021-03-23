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
ivwout = function(sbetaBk,ssigmaBk,sbetay,ssigmay,ivindk,multiplebiom){
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
