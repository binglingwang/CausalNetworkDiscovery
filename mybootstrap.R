## switched part 1 and 2

## use total rss  myfunction12

## reduced bootstrap tau size
# test = DataTrim[[t]]
# x = test[, 1]
# y = test[, 2]

TestFunc = function(x, y, ntau = 100, tau.x = NULL, tau.y = NULL, stat = "pearson"){
  
  # if(missing(ntau)){ntau = 50}
  
  # function to get model output
  lmout = function(xp, xc, tau){
    #xp = x; xc = y
    #x = d$x; y = d$y; tau = max(xp); tau = min(xp); tau = tau.t;
    # tau=tau.null.x[5]
    
    
    if(min(xp) >= tau) {#mean(head(sort(xp)))
      xp.1 = xp[which(xp < tau)]; xp.2 = xp
      xc.1 = xc[which(xp < tau)]; xc.2 = xc
      
      lm.1 = NA; lm.2 = lm(xc ~ xp)
      
      coef.1 = c(NA, NA); coef.2 = coef(lm.2) # coefficients
      
      resids = resid(lm.2)
      
      rss.1 = 0; rss.2 = sum(resids^2) # residual sum of squares
      
      r.1 = 0; r.2 = cor(xp.2, xc.2)
      #rsq.1 = 0; rsq.2 = (cor(xp.2, xc.2))^2 # correlation pearson
      
      rsq.1.spearman = 0; rsq.2.spearman = (cor(xp.2, xc.2, method = "spearman"))^2 # correlation spearman
      
    }else if(max(xp) <= tau) { #mean(tail(sort(xp)))
      xp.1 = xp; xp.2 = xp[which(xp > tau)]
      xc.1 = xc; xc.2 = xc[which(xp > tau)]
      
      lm.1 = lm(xc ~ xp); lm.2 = NA
      
      coef.1 = coef(lm.1); coef.2 = c(NA, NA)
      
      resids = resid(lm.1);
      
      rss.1 = sum(resids^2); rss.2 = 0
      
      r.1 = cor(xp.1, xc.1); r.2 = 0
      #rsq.1 = (cor(xp.1, xc.1))^2; rsq.2 = 0
      
      rsq.1.spearman = (cor(xp.1, xc.1, method = "spearman"))^2; rsq.2.spearman = 0
      
    }else{
      xp.1 = xp[which(xp < tau)]; xp.2 = xp[which(xp >= tau)]
      xc.1 = xc[which(xp < tau)]; xc.2 = xc[which(xp >= tau)]
      
      lm.1 = lm(xc.1 ~ xp.1); lm.2 = lm(xc.2 ~ xp.2)
      
      coef.1 = coef(lm.1); coef.2 = coef(lm.2)
      
      resids = c(resid(lm.1), resid(lm.2))
      
      rss.1 = sum((resid(lm.1))^2); rss.2 = sum((resid(lm.2))^2)
      
      r.1 = cor(xp.1, xc.1); r.2 = cor(xp.2, xc.2)
      #rsq.1 = (cor(xp.1, xc.1))^2; rsq.2 = (cor(xp.2, xc.2))^2
      
      rsq.1.spearman = (cor(xp.1, xc.1, method = "spearman"))^2;
      rsq.2.spearman = (cor(xp.2, xc.2, method = "spearman"))^2
      
    }
    
    
    # total sum of residuals
    n1 = length(xp.1)
    n2 = length(xp.2)
    Trss = rss.1 + rss.2
    
    # average Rsquare statistic
    if(is.na(r.1)){r.1 = 0};
    if(is.na(r.2)){r.2 = 0};
    if(is.na(rsq.1.spearman)){rsq.1.spearman = 0};
    if(is.na(rsq.2.spearman)){rsq.2.spearman = 0};
    if(sum(is.na(coef.1))){coef.1[is.na(coef.1)] = 0}
    if(sum(is.na(coef.2))){coef.2[is.na(coef.2)] = 0}
    
    RsqAve = (n1 * r.1^2 + n2 * r.2^2) / (n1 + n2)
    RsqAve.spearman = (n1 * rsq.1.spearman + n2 * rsq.2.spearman) / (n1 + n2)
    
    param = NULL
    param$xp.1 = xp.1
    param$xp.2 = xp.2
    param$xc.1 = xc.1
    param$xc.2 = xc.2
    param$n1 = n1
    param$n2 = n2
    param$r.1 = r.1
    param$r.2 = r.2
    param$coef1 = coef.1 #coef1
    param$coef2 = coef.2 #coef1
    param$residuals = resids
    param$Trss = Trss
    param$RsqAve = RsqAve
    param$RsqAve.spearman = RsqAve.spearman
    
    return(param)
  }
  

  # function to select the best estimate among different tau
  findMax = function(x, y) {
    
    # tau.x = seq(from = min(x), to = max(x), length.out = ntau)
    # tau.y = seq(from = min(y), to = max(y), length.out = ntau)
    if(is.null(tau.x)) tau.x = sapply(list(x), quantile, probs = round((0:ntau)/ntau, 3)) #[-c(1:3, (ntau-3):ntau)]
    if(is.null(tau.y)) tau.y = sapply(list(y), quantile, probs = round((0:ntau)/ntau, 3)) #[-c(1:3, (ntau-3):ntau)]
    
    m.x = sapply(tau.x, function(t){lmout(xp = x, xc = y, tau = t)})
    m.y = sapply(tau.y, function(t){lmout(xp = y, xc = x, tau = t)})
    
    
    ########################### rss total
    Trss.x = unlist(m.x["Trss", ])
    Trss.y = unlist(m.y["Trss", ])
    
    ind.x = which.min(Trss.x)
    ind.y = which.min(Trss.y)
    
    r1.x = unlist(m.x["r.1", ])
    r2.x = unlist(m.x["r.2", ])
    r1.y = unlist(m.y["r.1", ])
    r2.y = unlist(m.y["r.2", ])
    
    r.x = unlist(m.x[c("r.1", "r.2"), ind.x])
    r.y = unlist(m.y[c("r.1", "r.2"), ind.y])
    
    n.x = unlist(m.x[c("n1", "n2"), ind.x])
    n.y = unlist(m.y[c("n1", "n2"), ind.y])
    
    ########################### Coefficients
    coef.x = unlist(m.x[c("coef1", "coef2"), ind.x])
    coef.y = unlist(m.y[c("coef1", "coef2"), ind.y])
    names(coef.x) = names(coef.y) = c("a1", "b1", "a2", "b2")
    
    ########################### Residuals
    resids.x = unlist(m.x[c("residuals"), ind.x])
    resids.y = unlist(m.y[c("residuals"), ind.y])
    
    out = NULL
    #out.pearson = out.spearman = out.Cd = list()
    
    if("pearson" %in% stat){
      ########################### Rsquared pearson
      RsqAve.x = unlist(m.x["RsqAve", ])
      RsqAve.y = unlist(m.y["RsqAve", ])
      
      maxRsq.x = RsqAve.x[ind.x] #max(RsqAve.x, na.rm = TRUE)
      maxRsq.y = RsqAve.y[ind.y] #max(RsqAve.y, na.rm = TRUE)
      
      testStat = maxRsq.x / maxRsq.y
      
      if ( testStat >= 1 ) {
        out$parent = 1; out$testStat = testStat; out$resids = resids.y
        out$maxStat.p = maxRsq.x; out$tau.p = tau.x[ind.x]; out$r.p = r.x; 
        out$n.p = n.x; out$coef.p = coef.x
        out$maxStat.c = maxRsq.y; out$tau.c = tau.y[ind.y]; out$r.c = r.y; 
        out$n.c = n.y; out$coef.c = coef.y
        out$rp.1 = r1.x; out$rp.2 = r2.x; out$rc.1 = r1.y; out$rc.2 = r2.y; 
      } else {
        out$parent = 2; out$testStat =  1 / testStat; out$resids = resids.x
        out$maxStat.p = maxRsq.y; out$tau.p = tau.y[ind.y]; out$r.p = r.y; 
        out$n.p = n.y; out$n.c = n.x; out$coef.p = coef.y
        out$maxStat.c = maxRsq.x; out$tau.c = tau.x[ind.x]; out$r.c = r.x; out$coef.c = coef.x
        out$rp.1 = r1.y; out$rp.2 = r2.y; out$rc.1 = r1.x; out$rc.2 = r2.x; 
      }
     #out.pearson = out 
    }

    
    if("spearman" %in% stat){
      ########################### Rsquared pearson
      RsqAve.x = unlist(m.x["RsqAve.spearman", ])
      RsqAve.y = unlist(m.y["RsqAve.spearman", ])
      
      maxRsq.x = RsqAve.x[ind.x] #max(RsqAve.x, na.rm = TRUE)
      maxRsq.y = RsqAve.y[ind.y] #max(RsqAve.y, na.rm = TRUE)
      
      testStat = maxRsq.x / maxRsq.y
      
      if ( testStat >= 1 ) {
        out$parent = 1; out$testStat =  testStat; out$resids = resids.y
        out$maxStat.p = maxRsq.x; out$tau.p = tau.x[ind.x];
        out$n.p = n.x; out$r.p = r.x; out$coef.p = coef.x;
        out$maxStat.c = maxRsq.y; out$tau.c = tau.y[ind.y]; 
        out$n.c = n.y; out$r.c = r.y; out$coef.c = coef.y
        out$rp.1 = r1.x; out$rp.2 = r2.x; out$rc.1 = r1.y; out$rc.2 = r2.y; 
      } else {
        out$parent = 2; out$testStat =  1 / testStat; out$resids = resids.x
        out$maxStat.p = maxRsq.y; out$tau.p = tau.y[ind.y]; 
        out$n.p = n.y; out$r.p = r.y; out$coef.p = coef.y
        out$maxStat.c = maxRsq.x; out$tau.c = tau.x[ind.x]; 
        out$n.c = n.x; out$r.c = r.x; out$coef.c = coef.x
        out$rp.1 = r1.y; out$rp.2 = r2.y; out$rc.1 = r1.x; out$rc.2 = r2.x;
      }
      #out.spearman = out 
    }
    return(out)
  }
  
  result = findMax(x, y)
  return(result)
}



R2func = function(xp, xc, tau){

  #x = d$x; y = d$y; tau = max(xp); tau = min(xp); tau = tau.t;
  # tau=tau.null.x[5]

  if(min(xp) >= tau) {
    xp.1 = xp[which(xp < tau)]; xp.2 = xp
    xc.1 = xc[which(xp < tau)]; xc.2 = xc

    r.1 = 0; r.2 = (cor(xp.2, xc.2)) # correlation pearson

    r.1.spearman = 0; r.2.spearman = (cor(xp.2, xc.2, method = "spearman")) # correlation spearman

  } else if(max(xp) <= tau) {
    xp.1 = xp; xp.2 = xp[which(xp > tau)]
    xc.1 = xc; xc.2 = xc[which(xp > tau)]

    r.1 = (cor(xp.1, xc.1)); r.2 = 0

    r.1.spearman = (cor(xp.1, xc.1, method = "spearman")); r.2.spearman = 0

  }else{
    xp.1 = xp[which(xp < tau)]; xp.2 = xp[which(xp >= tau)]
    xc.1 = xc[which(xp < tau)]; xc.2 = xc[which(xp >= tau)]

    r.1 = (cor(xp.1, xc.1)); r.2 = (cor(xp.2, xc.2))

    r.1.spearman = (cor(xp.1, xc.1, method = "spearman"))^2;
    r.2.spearman = (cor(xp.2, xc.2, method = "spearman"))^2

  }


  # total sum of residuals
  n1 = length(xp.1)
  n2 = length(xp.2)

  # average Rsquare statistic
  if(is.na(r.1)){r.1 = 0};
  if(is.na(r.2)){r.2 = 0};
  if(is.na(r.1.spearman)){r.1.spearman = 0};
  if(is.na(r.2.spearman)){r.2.spearman = 0};


  RsqAve = (n1 * r.1^2 + n2 * r.2^2) / (n1 + n2)
  RsqAve.spearman = (n1 * r.1.spearman^2 + n2 * r.2.spearman^2) / (n1 + n2)

  param = NULL
  param$xp.1 = xp.1
  param$xp.2 = xp.2
  param$xc.1 = xc.1
  param$xc.2 = xc.2
  param$n1 = n1
  param$n2 = n2
  param$r.1 = r.1
  param$r.2 = r.2
  param$RsqAve = RsqAve
  param$RsqAve.spearman = RsqAve.spearman

  return(param)
}

# RRfunc = function(x, y, tau.x, tau.y, stat="pearson") {
#   if(stat == "pearson"){
#     R2.x = R2func(xp = x, xc = y, tau = tau.x)$RsqAve
#     R2.y = R2func(xp = y, xc = x, tau = tau.y)$RsqAve
#   } else{
#     R2.x = R2func(xp = x, xc = y, tau = tau.x)$RsqAve.spearman
#     R2.y = R2func(xp = y, xc = x, tau = tau.y)$RsqAve.spearman
#   }
#   
#   RR = max(R2.x/R2.y, R2.y/R2.x)
#   
#   return(RR)
# }

RRfunc = function(x, y, ntau=100, tau.x, tau.y, stat="pearson") {
  obj.x = R2func(xp = x, xc = y, tau = tau.x)
  obj.y = R2func(xp = y, xc = x, tau = tau.y)

  if(stat == "pearson"){
    r.x = c(obj.x$r.1, obj.x$r.2)
    r.y = c(obj.y$r.1, obj.y$r.2)
    R2.x = obj.x$RsqAve
    R2.y = obj.y$RsqAve
  }else{
    R2.x = obj.x$RsqAve.spearman
    R2.y = obj.y$RsqAve.spearman
  }

  #RR = max(R2.x/R2.y, R2.y/R2.x)

  out = NULL
  out$x = list(xp.1 = obj.x$xp.1, xp.2 = obj.x$xp.2, xc.1 = obj.x$xc.1, xc.2 = obj.x$xc.2)
  out$y = list(yp.1 = obj.y$xp.1, yp.2 = obj.y$xp.2, yc.1 = obj.y$xc.1, yc.2 = obj.y$xc.2)
  out$n.x = c(obj.x$n1, obj.x$n2)
  out$n.y = c(obj.y$n1, obj.y$n2)
  out$r.x = r.x
  out$r.y = r.y
  out$RR = max(R2.x/R2.y, R2.y/R2.x)
  return(out)
}

bootstrapfunc = function(x, y, ntau = 100, tau.x, tau.y, func){
  n = length(x)
  s = sample(n, replace = T);
  res = func(x[s], y[s], ntau, tau.x, tau.y)
}
# x = res10N20E80P$dagList[[1]]$DAGdata[, 1]
# y = res10N20E80P$dagList[[1]]$DAGdata[, 2]

mybootstrap = function(x, y, k=100, ntau = 100, xlabel = "1", ylabel = "2", tau.x = NULL, tau.y = NULL, stat = "pearson",  #func = myfunction, 
                       tau.boot.step = NULL, bootstrap.parameter = FALSE, reest_bp = TRUE, approx = TRUE, testonly = TRUE,
                       PrintProcess = FALSE) {
  func = TestFunc
  ################################################################# Bootstrap
  # if(missing(k)){k = 100}
  # if(missing(ntau)){ntau = 100}
  n = length(x)
  object = func(x = x, y = y, ntau = ntau, tau.x = tau.x, tau.y = tau.y, stat = stat)
  
  if(object$parent == 1){
    xp = x; xc = y
    tx = object$tau.p; ty = object$tau.c
    parent.label = xlabel; child.label = ylabel
  }  else {
    xp = y; xc = x
    tx = object$tau.c; ty = object$tau.p
    parent.label = ylabel; child.label = xlabel
  }
  
  tau.p = object$tau.p
  tau.c = object$tau.c
  
  if(min(xp) >= tau.p) {
    xp.1 = xp[which(xp < tau.p)]; xp.2 = xp
    xc.1 = xc[which(xp < tau.p)]; xc.2 = xc
  }else if(max(xp) <= tau.p) {
    xp.1 = xp; xp.2 = xp[which(xp > tau.p)]
    xc.1 = xc; xc.2 = xc[which(xp > tau.p)]
  }else{
    xp.1 = xp[which(xp < tau.p)]; xp.2 = xp[which(xp >= tau.p)]
    xc.1 = xc[which(xp < tau.p)]; xc.2 = xc[which(xp >= tau.p)]
  }
  
  coefs = object$coef.p
  xc.1.hat = coefs[1] + coefs[2] * xp.1
  xc.2.hat = coefs[3] + coefs[4] * xp.2
  
  a = xc.1.hat[which.min(xp.1)]
  b = xc.1.hat[which.max(xp.1)]
  c = xc.2.hat[which.min(xp.2)]
  d = xc.2.hat[which.max(xp.2)]
  
  if( max(a, b) <= min(c, d) || min(a, b) >= max(c, d) ) {
    dist.c = 0
  }else if(abs(max(a, b) - min(c, d)) <= abs(min(a, b) - max(c, d))) {
    dist.c = max(a, b) - min(c, d)
  }else {
    dist.c = min(a, b) - max(c, d)
  }
  
  #xnull = xp; ynull = xc
  xc[which(xp < tau.p)] = xc.1 - dist.c
  
  # if (object$parent == 1) {
  #   x.null = xp.null; y.null = xc.null
  # } else {
  #   x.null = xc.null; y.null = xp.null
  # }
  
  # plot(xp, xc, cex = 0.1)
  # plot(xc, xp, cex = 0.1)
  # plotseg(x.null, y.null, obj.null)
  # plot(x.null, y.null, cex = .1)

  
  
  if(reest_bp){
    obj.null = func(x = xp, y = xc, ntau = ntau, stat = stat)
  }else{
    ###########################################
    # estimate null at same breakpoints 202005
    ###########################################
    obj.null = func(x = xp, y = xc, ntau = ntau, tau.x = tau.p, tau.y = tau.c, stat = stat)
  }

  
  # plotseg(xp, xc, obj.null)
  if(obj.null$parent == 1) {
    tx.null = obj.null$tau.p; ty.null = obj.null$tau.c;
    x.null = xp; y.null = xc
  } else {
    tx.null = obj.null$tau.c; ty.null = obj.null$tau.p; 
    x.null = xp; y.null = xc
  }
  
  dat.null = cbind(x.null, y.null)

  if(min(c(obj.null$n.p, obj.null$n.c))<=3) approx = F
  # source(file = "~/Bitbucket_Repos/chipseqproject/Functions/plotseg.R")
  # plotseg(x1 = x.null, x2 = y.null, object = obj.null)
  
  ntau.boot = 2 * tau.boot.step + 1
  
  res.null = NULL

  if(approx){
    
    n1.p = obj.null$n.p[1]; n2.p = obj.null$n.p[2]
    m1.p = atanh(obj.null$r.p[1]); m2.p = atanh(obj.null$r.p[2])
    s1.p = 1/sqrt(n1.p-3); s2.p = 1/sqrt(n2.p-3)
    
    n1.c = obj.null$n.c[1]; n2.c = obj.null$n.c[2]
    m1.c = atanh(obj.null$r.c[1]); m2.c = atanh(obj.null$r.c[2])
    s1.c = 1/sqrt(n1.c-3); s2.c = 1/sqrt(n2.c-3)
    
    phi1.p = rnorm(n = 100, mean = m1.p, sd = s1.p)
    phi2.p = rnorm(n = 100, mean = m2.p, sd = s2.p)
    
    phi1.c = rnorm(n = 100, mean = m1.c, sd = s1.c)
    phi2.c = rnorm(n = 100, mean = m2.c, sd = s2.c)
     
    # rsq1.x = tanh(phi1.x)^2
    # rsq2.x = tanh(phi2.x)^2
    # rsq1.y = tanh(phi1.y)^2
    # rsq2.y = tanh(phi2.y)^2
    # 
    # 
    # R2.x = ((r$n.x[1] * rsq1.x + r$n.x[2] * rsq2.x) / (r$n.x[1] + r$n.x[2]))
    # R2.y = ((r$n.y[1] * rsq1.y + r$n.y[2] * rsq2.y) / (r$n.y[1] + r$n.y[2]))
    
    testStat = object$testStat
    testStat.boot = (n1.p * tanh(phi1.p)^2 + n2.p * tanh(phi2.p)^2) / (n1.c * tanh(phi1.c)^2 + n2.c * tanh(phi2.c)^2)
    testStat.boot[testStat.boot<1] = 1/testStat.boot[testStat.boot<1]
    
    # testStat.boot = (n1.x * r1.x^2 + n2.x * r2.x^2) / (n1.y * r1.y^2 + n2.y * r2.y^2)
    # testStat.boot[testStat.boot<1] = 1/testStat.boot[testStat.boot<1]
    
    p.r2 = ifelse(testStat == 1, 1, sum(testStat.boot >= testStat) / k)
    
    
    ## boostrap original data for coefficients significance
    
    
    dat = cbind(x, y)
    colnames(dat) = c(xlabel, ylabel)
    
    outList = NULL
    
    outList$testStat = testStat
    outList$tau = tau.p
    outList$p.r2 = p.r2
    outList$parent = parent.label
    outList$child = child.label
    outList$obj = object
    outList$obj.null = obj.null
    outList$testStat.boot = testStat.boot
    outList$dat.null = dat.null
    outList$dat = dat
    
  }else{
    
  if(testonly){
    
    # for(i in 1:k){
    #   s = sample(n, replace = T);
    #   # func.s = suppressWarnings(func(x = x.null[s], y = y.null[s], tau.x = tx.null, tau.y = ty.null,
    #   #                                ntau = ntau, stat = stat))
    #   func.s = RRfunc(x = x.null[s], y = y.null[s], tau.x = tx.null, tau.y = ty.null)
    #   res.null = cbind(res.null, func.s);
    #   if(PrintProcess == TRUE) print(i)
    # }
    
    res.null = replicate(k, bootstrapfunc(x = x.null, y = y.null, tau.x = tx.null, tau.y = ty.null, func = RRfunc))
    
    ## r2 average pvalue
    testStat = object$testStat
    testStat.boot = as.vector(unlist(res.null["RR",]))
    # hist(testStat.boot)
    #parent.pct = sum(as.vector(unlist(res.null["parent",])) == object$parent) / k
    
    p.r2 = ifelse(testStat == 1, 1, sum(testStat.boot >= testStat) / k)
    
    ## boostrap original data for coefficients significance
    
    
    dat = cbind(x, y)
    colnames(dat) = c(xlabel, ylabel)
    # # p value of statistic under null
    # if(p.r2 > 0.5) p.r2 = (1-p.r2)
    
    
    outList = NULL
    
    outList$testStat = testStat
    outList$tau = tau.p
    outList$p.r2 = p.r2
    outList$parent = parent.label
    outList$child = child.label
    outList$obj = object
    outList$obj.null = obj.null
    outList$testStat.boot = testStat.boot
    outList$dat.null = dat.null
    outList$dat = dat
    
  }else {
    if(is.null(tau.boot.step)){
      
      #ptm <- proc.time()
      # for(i in 1:k){
      #   s = sample(1:n, n, replace = TRUE);
      #   func.s = suppressWarnings(func(x = x.null[s], y = y.null[s], #tau.x = tx.null, tau.y = ty.null,
      #                                  ntau = ntau, stat = stat))
      #   res.null = cbind(res.null, func.s);
      #   if(PrintProcess == TRUE) print(i)
      # }
      res.null = replicate(k, bootstrapfunc(x = x.null, y = y.null, ntau = 100, tau.x = NULL, tau.y = NULL, func = func))
      
      #proc.time() - ptm
    }else{
      
      tau.null.x = quantile( x.null, seq(max(0, round(sum(x.null <= tx.null)/n - tau.boot.step / ntau, 3)),
                                         min(1, round(sum(x.null <= tx.null)/n + tau.boot.step / ntau, 3)),
                                         length.out = ntau.boot) );
      tau.null.y = quantile( y.null, seq(max(0, round(sum(y.null <= ty.null)/n - tau.boot.step / ntau, 3)),
                                         min(1, round(sum(y.null <= ty.null)/n + tau.boot.step / ntau, 3)),
                                         length.out = ntau.boot) );
      
      for(i in 1:k){
        s = sample(1:n, n, replace = TRUE);
        func.s = suppressWarnings(func(x = x.null[s], y = y.null[s], ntau = ntau,
                                       tau.x = tau.null.x, tau.y = tau.null.y, stat = stat));
        
        
        if(func.s$parent == 1) {
          null.p = x.null; null.c = y.null
          tau.null.p = tau.null.x; tau.null.c = tau.null.y; 
        } else {
          null.p = y.null; null.c = x.null
          tau.null.p = tau.null.y; tau.null.c = tau.null.x;
        }
        
        it = 0
        while( (func.s$tau.p != min(null.p)) & (func.s$tau.p != max(null.p)) &
               (func.s$tau.p == min(tau.null.p) | func.s$tau.p == max(tau.null.p)) ){
          tau.null.p = quantile(null.p, seq(max(0, round(sum(null.p <= func.s$tau.p)/n - tau.boot.step / ntau, 3)),
                                            min(1, round(sum(null.p <= func.s$tau.p)/n + tau.boot.step / ntau, 3)), 
                                            length.out = ntau.boot))
          tau.null.c = quantile(null.c, seq(max(0, round(sum(null.c <= func.s$tau.c)/n - tau.boot.step / ntau, 3)),
                                            min(1, round(sum(null.c <= func.s$tau.c)/n + tau.boot.step / ntau, 3)), 
                                            length.out = ntau.boot))
          func.s = suppressWarnings(func(x = null.p[s], y = null.c[s], ntau = ntau,
                                         tau.x = tau.null.p, tau.y = tau.null.c, stat = stat))
          it = it + 1
          if(it > 10) break
        }
        
        res.null = cbind(res.null, func.s)
        if(PrintProcess == TRUE) print(i)
      }
    }
    
    
    ## r2 average pvalue
    testStat = object$testStat
    testStat.boot = as.vector(unlist(res.null["testStat",]))
    parent.pct = sum(as.vector(unlist(res.null["parent",])) == object$parent) / k
    
    p.r2 = ifelse(testStat == 1, 1, sum(testStat.boot >= testStat) / k)
    
    
    ## boostrap original data for coefficients significance
    
    if (bootstrap.parameter == TRUE) {
      
      res.orig = NULL
      if(is.null(tau.boot.step)){
        
        for(i in 1:k){
          s = sample(1:n, n, replace = TRUE);
          func.s = suppressWarnings(func(x = x[s], y = y[s], ntau = ntau, stat = stat));
          res.orig = cbind(res.orig, func.s);
          if(PrintProcess == TRUE) print(i)
        }
        
      }else{
        
        tau.orig.x = quantile( x, seq(max(0, round(sum(x <= tx)/n - tau.boot.step / ntau, 3)),
                                      min(1, round(sum(x <= tx)/n + tau.boot.step / ntau, 3)), 
                                      length.out = ntau.boot) );
        tau.orig.y = quantile( y.null, seq(max(0, round(sum(y <= ty)/n - tau.boot.step / ntau, 3)),
                                           min(1, round(sum(y <= ty)/n + tau.boot.step / ntau, 3)), 
                                           length.out = ntau.boot) );
        
        for(i in 1:k){
          s = sample(1:n, n, replace = TRUE);
          func.s = suppressWarnings(func(x = x[s], y = y[s], ntau = ntau, stat = stat,
                                         tau.x = tau.orig.x, tau.y = tau.orig.y));
          
          if (func.s$parent == 1) {
            orig.p = x; orig.c = y
            tau.orig.p = tau.orig.x; tau.orig.c = tau.orig.y; 
          } else {
            orig.p = y; orig.c = x
            tau.orig.p = tau.orig.y; tau.orig.c = tau.orig.x;
          }
          
          it = 0
          while( (func.s$tau.p != min(orig.p)) & (func.s$tau.p != max(orig.p)) &
                 (func.s$tau.p == min(tau.orig.p) | func.s$tau.p == max(tau.orig.p)) ){
            tau.orig.p = quantile(orig.p, seq(max(0, round(sum(orig.p <= func.s$tau.p)/n - tau.boot.step / ntau, 3)),
                                              min(1, round(sum(orig.p <= func.s$tau.p)/n + tau.boot.step / ntau, 3)), 
                                              length.out = ntau.boot))
            tau.orig.c = quantile(orig.c, seq(max(0, round(sum(orig.c <= func.s$tau.c)/n - tau.boot.step / ntau, 3)),
                                              min(1, round(sum(orig.c <= func.s$tau.c)/n + tau.boot.step / ntau, 3)), 
                                              length.out = ntau.boot))
            func.s = suppressWarnings(func(x = orig.p[s], y = orig.c[s], ntau = ntau, 
                                           tau.x = tau.orig.p, tau.y = tau.orig.c, stat = stat))
            it = it + 1
            if(it > 10) break
          }
          
          res.orig = cbind(res.orig, func.s)
          if(PrintProcess == TRUE) print(i)
        }
      }
      
      coef.boot = matrix(unlist(res.orig["coef.p",]), ncol = 4, byrow = TRUE)
      colnames(coef.boot) = c("a1", "b1", "a2", "b2")
      
      p.a1 = sum( coef.boot[, 1] >= 0 ) / k
      p.a1 = ifelse( p.a1 <= 0.5, 2 * p.a1, 2 * (1 - p.a1) )
      
      p.b1 = sum( coef.boot[, 2] >= 0 ) / k
      p.b1 = ifelse( p.b1 <= 0.5, 2 * p.b1, 2 * (1 - p.b1) )
      
      p.a2 = sum( coef.boot[, 3] >= 0 ) / k
      p.a2 = ifelse( p.a2 <= 0.5, 2 * p.a2, 2 * (1 - p.a2) )
      
      p.b2 = sum( coef.boot[, 4] >= 0 ) / k
      p.b2 = ifelse( p.b2 <= 0.5, 2 * p.b2, 2 * (1 - p.b2) )
      
      p.b12 = sum( (coef.boot[, 4] - coef.boot[, 2]) >= 0 ) / k
      p.b12 = ifelse( p.b12 <= 0.5, 2 * p.b12, 2 * (1 - p.b12) )
    }
    
    
    
    dat = cbind(x, y)
    colnames(dat) = c(xlabel, ylabel)
    # # p value of statistic under null
    # if(p.r2 > 0.5) p.r2 = (1-p.r2)
    
    
    outList = NULL
    
    outList$testStat = testStat
    outList$tau = tau.p
    outList$p.r2 = p.r2
    outList$parent = parent.label
    outList$child = child.label
    outList$parent.pct = parent.pct
    
    if (bootstrap.parameter == TRUE) {
      outList$coef.boot = coef.boot
      outList$p.a1 = p.a1
      outList$p.b1 = p.b1
      outList$p.a2 = p.a2
      outList$p.b2 = p.b2
      outList$p.b12 = p.b12
    }
    
    outList$obj = object
    outList$obj.null = obj.null
    outList$testStat.boot = testStat.boot
    outList$dat.null = dat.null
    outList$dat = dat
    
  }
}
  
  return(outList)
}



#save("mybootstrap", file="~/Bitbucket_Repos/chipseqproject/Functions/mybootstrap.Rdata")


