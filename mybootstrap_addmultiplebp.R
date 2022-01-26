
mybootstrap = function(x, y, k=100, ntau = 100, xlabel = "1", ylabel = "2", stat = "pearson",  #func = myfunction, 
                       tau.boot.step = NULL, bootstrap.parameter = FALSE, bootstrap.null = TRUE, approx = TRUE, testonly = TRUE,
                       PrintProcess = FALSE, bp_method = c("single")) {

  ################################################################# Bootstrap
  # if(missing(k)){k = 100}
  # if(missing(ntau)){ntau = 100}
  n = length(x)
  
  if(bp_method == "single"){
    func = TestFunc
    object = func(x = x, y = y, ntau = ntau, stat = stat)
  }else{
    func = test_multiple_bp
    object = func(x = x, y = y)
  }
  
  
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
  
  
  if(bp_method == "single"){
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
    
    obj.null = func(x = xp, y = xc, ntau = ntau, stat = stat)
    
  }else{
    intcpt = object$model.p$coefficients[1]
    coefs = sapply(1:length(object$model.p$coefficients[-1]), function(i) {
      sum(object$model.p$coefficients[-1][1:i]) 
    })
    
    id.group = object$model.p$id.group
    
    # dist.c = NULL
    for(i in 2:length(coefs)){
      # if(i==1){
      #   yhat_bp[i] = coefs[i] + coefs[i+1] * tau.p[i]
      # }else{
      #   yhat_bp[i] = yhat_bp[i-1] + sum(coefs[2:(i+1)]) * (tau.p[i] - tau.p[i-1])
      # }
      
      if(sign(coefs[i]) != sign(coefs[1])){
        xc[id.group==i] = xc[id.group==i] - 2*coefs[i]*(xp[id.group==i]-xp[id.group==i][1])
        xc[id.group>i] = xc[id.group>i] - 2*coefs[i]*(tail(xp[id.group==i],n=1)-xp[id.group==i][1])
        
        # coefs[i] = - coefs[i]
      }
      
    }
    
    obj.null = func(x = xp, y = xc)
  }

  plot(xp, xc, cex = 0.5)
  plot(object$model.p, add = T, col = "red")
  lines(obj$model.p, col='orange')
  
  
  # plotseg(xp, xc, obj.null)
  if(obj.null$parent == 1) {
    tx.null = obj.null$tau.p; ty.null = obj.null$tau.c;
    x.null = xp; y.null = xc
  } else {
    tx.null = obj.null$tau.c; ty.null = obj.null$tau.p; 
    x.null = xp; y.null = xc
  }
  
  dat.null = cbind(x.null, y.null)
  

  if(bp_method == "multiple"){
    
    
  }
  
  
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



