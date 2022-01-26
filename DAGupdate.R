##########################################
# Added method 2
# First estimate an empty graph to fix the nonlinear edges;
# Then apply PC or sbn to the graph to get a final estimated graph
##########################################

#########################################################
# Add conditional independence test
#########################################################

#setwd("C:/Users/bwang4/Desktop/Personal")
#setwd("Bitbucket_Repos/chipseqproject")
#source("/Functions/mydatagen.R")
#source("/Functions/myfunction.R")
source("Bitbucket_Repos/causal-network/FunctionsV2/mybootstrap.R")
#source("Functions/myplot.R")
#source("~/Bitbucket_Repos/chipseqproject/Functions/applyOrientationRules.R")


packages = function(x){
  #x = as.character(match.call()[[2]])
  if(!require(x,character.only=TRUE)){
    install.packages(pkgs=x, quiet = TRUE)
    require(x,character.only=TRUE)
  }else{
    require(x,character.only=TRUE)
  }
}

pkgBiocLite = function(x){
  #x = as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    source("http://bioconductor.org/biocLite.R")
    biocLite("BiocInstaller", suppressUpdates = TRUE)
    biocLite(x, suppressUpdates = TRUE)
    require(x,character.only=TRUE)
  }else{
    require(x,character.only=TRUE)
  }
}

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocInstaller", suppressUpdates = TRUE)
pkgBiocLite("graph")
pkgBiocLite("RBGL")
packages("pcalg")
packages("sparsebn")
packages("igraph")
packages("bnlearn")

###################################################################################################
DAGupdate = function(dat = NULL, labels = NULL, sbn_alpha = 0.2, pc_alpha = 0.01, skel_alpha = NULL,
                     bootstrap = mybootstrap, p_value = 0.01, stat = "pearson", iter = 100, 
                     ntau = 100, tau.boot.step = 3, pc.maxdegree = Inf, V1 = TRUE, V2 = TRUE,
                     UsePC = FALSE, UseSparsebn = FALSE, UseSkeleton = FALSE, UseEmptyG = TRUE, 
                     matchPC = FALSE, CheckAllPairs = TRUE, approx = TRUE, testonly = TRUE, 
                     conIndAlpha = 0.01, conIndTest = FALSE, test.method = "pearson", seed = 12345,
                     printProgress = FALSE) {
  
  set.seed(seed)
  
  #P = ncol(dat)
  # varible labels
  labels = if (is.null(labels)){
    paste("X", 1:ncol(dat), sep = "")
  } else {
    labels
  }
  colnames(dat) = labels
  
  func = function(dat, labels, AdjMat, EdgeList){  
    #print("Updating Graph")
    
    EdgeList.nl = edgeList(matrix(0, ncol(dat), ncol(dat)))
    AdjMat.nl = matrix(0, ncol(dat), ncol(dat))
    
    names(EdgeList) = names(EdgeList.nl) = labels
    colnames(AdjMat) = rownames(AdjMat) = colnames(AdjMat.nl) = rownames(AdjMat.nl) = labels
    # convert to adjacency matrix
    #AdjMat =  G.adjmat 
    #AdjMat = skeleton.adjmat#G.adjmat
    
    dat.resid = dat
    # regress each node on its parents to get residuals
    for(i in 1:ncol(dat)) {
      
      # parents of node 
      parent_node = which(AdjMat[, i] - AdjMat[i, ] == 1)
      
      # residual of node i if has parents
      dat.resid[, i] = if(sum(parent_node)) {
        #print(i)
        mod = paste(labels[i], "~", paste(labels[parent_node], collapse = " + "))
        resid(lm(mod, data = as.data.frame(dat)))
      }else {
        dat[, i]
      }
    }
    
    AdjMat.update = AdjMat
    # EdgeList = EdgeList
    
    # which edges are undirected
    undirectedInd = matrix(which((AdjMat == t(AdjMat) & (AdjMat != 0)), 
                                 arr.ind = TRUE), ncol = 2)
    undirectedPairs = matrix(undirectedInd[undirectedInd[, 1] < undirectedInd[, 2], ], ncol = 2)
    
    # undirectedPairs = undirectedPairs
    NpInd = NcInd =  NULL
    #dat.update = vector("list", ncol(AdjMat.update))
    # PairRes = NULL
    undirectedRes = vector("list", nrow(undirectedPairs))
    
    pairInd.update = 1:nrow(undirectedPairs)
    firstround = TRUE
    s=1
    while( s <= nrow(undirectedPairs) ) {
      
      for(k in pairInd.update) {
        # node i and node j
        i =  undirectedPairs[k, 1]; j = undirectedPairs[k, 2]
        
        Vi = dat.resid[, i]; Vj = dat.resid[, j]
        
        parent_i = which(AdjMat.update[, i] - AdjMat.update[i, ] == 1)
        parent_j = which(AdjMat.update[, j] - AdjMat.update[j, ] == 1)
        
        
        # determine reversible edge
        if( sum(NpInd %in% parent_i) | sum(NpInd %in% parent_j) | firstround ) {
          undirectedRes[[k]] = bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
                                         k = iter, ntau = ntau, tau.boot.step = tau.boot.step, 
                                         stat = stat, approx = approx, testonly = testonly)
        }
      }
      
      #colnames(PairRes) = pairnames
      
      # Sort by pvalue ascending, testStat descending
      #edge.update = order(unlist(PairRes[c("p.r2"), ]), -unlist(PairRes[c("testStat"), ]))[1]
      edge.update = order(sapply(undirectedRes, function(x) x$p.r2), 
                          -sapply(undirectedRes, function(x) x$testStat))[[1]]

      #undirectedRes[[2]]$p.r2
      if( undirectedRes[[edge.update]]$p.r2 <= p_value ) {
        # edge parent & child
        Np = which(labels %in% undirectedRes[[edge.update]]$parent)
        Nc = which(labels %in% undirectedRes[[edge.update]]$child)
        adjmat.valid.check = AdjMat.update
        # remove the edge from child to parent
        adjmat.valid.check[Nc, Np] = 0
        # remove undirected edges to check if there is directed cycle
        adjmat.valid.check[which(adjmat.valid.check == t(adjmat.valid.check))] = 0
        valid = isValidGraph(t(adjmat.valid.check), type = "dag")
       
        if( valid ) {
          #########################################################
          # Add conditional independence test
          #########################################################
          if(conIndTest) {
            cutoff = qnorm(1-conIndAlpha/2)
            
            # parent set of the child node
            Sind = which(AdjMat.update[, Nc]==1 & AdjMat.update[, Nc]!= t(AdjMat.update)[ , Nc])
            conIndDat = dat[, c(Np, Nc, Sind)] # test on original variable
            conIndDat.resid = dat.resid[, c(Np, Nc, Sind)] # estimates come from residuals
            if(length(Sind)!=0) {
              S = 1:length(Sind) +2
            }else {
              S = NULL
            }
            
            test.tau = undirectedRes[[edge.update]]$tau
            conIndDat1 = conIndDat[which(conIndDat.resid[, 1] <= test.tau), ]
            conIndDat2 = conIndDat[which(conIndDat.resid[, 1] > test.tau), ]
            
            test1 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat1, method=test.method), n=nrow(conIndDat1), cutoff=cutoff)
            test2 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat2, method=test.method), n=nrow(conIndDat2), cutoff=cutoff)
            test = test1 + test2
          }else{
            test = FALSE
          }
          
          
          if(!test) {
            ## if conditional test rejected for both parts
            
            NpInd = cbind(NpInd, Np)
            NcInd = cbind(NcInd, Nc)
            AdjMat.update[Nc, Np] = 0
            # update child node to residuals
            #dat.update[[Np]][[Nc]] = undirectedRes[[edge.update]]$obj$resids
            dat.resid[, Nc] = undirectedRes[[edge.update]]$obj$resids
            
            EdgeList.nl[[Nc]] = c(EdgeList.nl[[Nc]], Np)
            AdjMat.nl[Np, Nc] = 1
            
            # undirectedRes.update = undirectedRes[[-edge.update]]
            pairInd.update = pairInd.update[-which(pairInd.update == edge.update)]

          }
        }
      }
      
  
      firstround = FALSE
      s=s+1
      
      if( undirectedRes[[edge.update]]$p.r2 <= p_value ) {
        undirectedRes[[edge.update]]$p.r2 = Inf
      } else break
      # if(is.na(sum(undirectedPairs))) break
      # print(s)
      
    }
    
    
    # update Graph edgeList
    EdgeList.update = edgeList(AdjMat.update)
    
    whitelist = matrix(labels[which(AdjMat.nl==1, arr.ind = TRUE)], ncol = 2)
    #whitelist = matrix(labels[which(AdjMat.nl==1, arr.ind = TRUE)], ncol = 2)
    G.bn = to_bn(to_graphNEL(EdgeList.update))
    G.bn$learning$whitelist = whitelist
    cpdag.bn = cpdag(G.bn, moral = FALSE, wlbl = TRUE)
    AdjMat.update = amat(cpdag.bn)
    EdgeList.update = edgeList(AdjMat.update)
    
    
    AdjMat.update.all = AdjMat.update
    AdjMat.nl.all = AdjMat.nl
    #EdgeList.nl.all = EdgeList.nl
    
    noEdgeInd = which(AdjMat.update.all + t(AdjMat.update.all) == 0, arr.ind = TRUE)
    noEdgePairs =  matrix(noEdgeInd[noEdgeInd[, 1] < noEdgeInd[, 2], ], ncol = 2)
    
    if( CheckAllPairs & length(noEdgePairs)>0 ) {

      # dat_resid = dat
      res_noEdgePairs = list()
      for(k in 1:nrow(noEdgePairs)) {
        # node i and node j
        i =  noEdgePairs[k, 1]; j = noEdgePairs[k, 2]
        
        Vi = dat.resid[, i]; Vj = dat.resid[, j]
        
        # determine reversible edge
        
        res_noEdgePairs[[k]] = bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
                                         k = iter, ntau = ntau, tau.boot.step = tau.boot.step, 
                                         stat = stat, testonly = testonly)
        # PairRes = cbind(PairRes, bootstrap(x = Vi, y = Vj, xlabel = labels[i], ylabel = labels[j],
        #                          k = iter, ntau = ntau, tau.boot.step = tau.boot.step))
        
      }
      
      p_noEdgePairs = sapply(res_noEdgePairs, function(x) x$p.r2)
      # fdr_noEdgePairs = round(sapply(1:length(p_noEdgePairs), FUN = function(i) {
      #   length(p_noEdgePairs) * p_noEdgePairs[i] / length(which(p_noEdgePairs <= p_noEdgePairs[i]))}), 3)
      
      
      for(k in which(p_noEdgePairs <= p_value)){
        Np = which(labels %in% res_noEdgePairs[[k]]$parent)
        Nc = which(labels %in% res_noEdgePairs[[k]]$child)
        
        adjmat.valid.check = AdjMat.update.all
        adjmat.valid.check[Np, Nc] = 1
        adjmat.valid.check[which(adjmat.valid.check == t(adjmat.valid.check))] = 0
        valid = isValidGraph(t(adjmat.valid.check), type = "dag")
        
        if( valid ) {
          #########################################################
          # Add conditional independence test
          #########################################################
          if(conIndTest) {
            cutoff = qnorm(1-conIndAlpha/2)
            
            # parent set of the child node
            Sind = which(AdjMat.update.all[, Nc]==1 & AdjMat.update.all[, Nc]!= t(AdjMat.update.all)[ , Nc])
            conIndDat = dat[, c(Np, Nc, Sind)]
            conIndDat.resid = dat.resid[, c(Np, Nc, Sind)]
            
            if(length(Sind)!=0) {
              S = 1:length(Sind) +2
            }else {
              S = NULL
            }

            test.tau = res_noEdgePairs[[k]]$tau
            conIndDat1 = conIndDat[which(conIndDat.resid[, 1] <= test.tau), ]
            conIndDat2 = conIndDat[which(conIndDat.resid[, 1] > test.tau), ]
            
            test1 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat1, method=test.method), n=nrow(conIndDat1), cutoff=cutoff)
            test2 = condIndFisherZ(x=1, y=2, S=S, C=cor(conIndDat2, method=test.method), n=nrow(conIndDat2), cutoff=cutoff)
            test = test1 + test2
          }else{
            test = FALSE
          }

          
          if(!test) {
            ## if conditional test rejected for both parts
            AdjMat.update.all[Np, Nc] = 1
            AdjMat.nl.all[Np, Nc] = 1
          }
        }
      }
      
      
    } else {
      AdjMat.update.all = AdjMat.nl.all = matrix(0, ncol(dat), ncol(dat))
      EdgeList.update.all = EdgeList.nl.all = edgeList(AdjMat.update.all)
    }
    

    EdgeList.update.all = edgeList(AdjMat.update.all)
    
    
    whitelist.all = matrix(labels[which(AdjMat.nl.all==1, arr.ind = TRUE)], ncol = 2)
    # whitelist.all = matrix(labels[which(AdjMat.nl.all==1, arr.ind = TRUE)], ncol = 2)
    G.bn.all = to_bn(to_graphNEL(EdgeList.update.all))
    G.bn.all$learning$whitelist = whitelist.all
    cpdag.bn.all = cpdag(G.bn.all, moral = FALSE, wlbl = TRUE)
    # this step may not result a true cpdag since some edges are restricted
    
    AdjMat.update.all = amat(cpdag.bn.all)
    EdgeList.update.all = edgeList(AdjMat.update.all)
    EdgeList.nl.all = edgeList(AdjMat.nl.all)
    
    # out = list()
    # out$G_update = EdgeList
    # out$adjacency.matrix = AdjMat.update
    # return(out)
    # return(list(AdjMat.update = AdjMat.update, 
    #             EdgeList.update = EdgeList.update))
    return(list(cpdagBN = cpdag.bn, cpdagAddBN = cpdag.bn.all,
                EdgeList=EdgeList.update, AdjMatrix=AdjMat.update, 
                EdgeListNL = EdgeList.nl, AdjMatrixNL = AdjMat.nl,
                EdgeListAdd = EdgeList.update.all, AdjMatrixAdd = AdjMat.update.all,
                EdgeListAddNL = EdgeList.nl.all, AdjMatrixAddNL = AdjMat.nl.all
                ))
  }
  
  out = res1 = res2 = res_sbn = res_pc = res_skeleton = res_emptyG= NULL
  
  
  
  if(V1) {
    # print("Method 1")
    
    if(UseEmptyG) {
      if(printProgress) cat("NL only...")      
      
      G.amat = matrix(0, ncol(dat), ncol(dat))
      G.edgelist = edgeList(G.amat)
      res_emptyG = func(dat = dat, labels = labels, 
                        AdjMat = G.amat, EdgeList = G.edgelist)
      #plot(res_nobase$EdgeListAdd, layout = test.layout)
    }
    
    
    # Use pcalg graph
    if(UsePC) {
      # pc algorithm Graph
      #print("updating PC")
      if(printProgress) cat("PC + NL...")

      ### retry until find an extendable DAG or randomly generate one
      pc.learn = pc(suffStat = list(C = cor(dat), n = nrow(dat)),
                    indepTest = gaussCItest, m.max = pc.maxdegree,## indep.test: partial correlations
                    alpha = pc_alpha, labels = labels, u2pd = "retry")
      
      # isValidGraph(amat = as(pc.learn, "matrix"), type = "cpdag")
      G.edgelist = edgeList(t(as(pc.learn, "matrix")))
      G.learn = to_bn(G.edgelist)
      
      #G.learn = pc.stable(as.data.frame(dat), alpha = pc_alpha)
      cpdag.learn = cpdag(G.learn, moral = FALSE)
      cpdag.amat = amat(cpdag.learn)
      cpdag.edgelist = edgeList(cpdag.amat)
      #isValidGraph(t(cpdag.amat), "cpdag")
      
      skeleton.learn= bnlearn::skeleton(cpdag.learn)
      skeleton.amat = amat(skeleton.learn)
      skeleton.edgelist = edgeList(skeleton.amat)
       
      res_pc$G_est = list(G.learn = G.learn, cpdag.learn = cpdag.learn, EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
                          skeleton.EdgeList = skeleton.edgelist, skeleton.AdjMatrix = skeleton.amat)
      res_pc$G_update = func(dat = dat, labels = labels, 
                             AdjMat = cpdag.amat, EdgeList = cpdag.edgelist)
      if(UseSkeleton) {
        res_pc$skeleton_update = func(dat = dat, labels = labels, 
                                      AdjMat = skeleton.amat, EdgeList = skeleton.edgelist)
      }
      # pc.learn = pc(suffStat = list(C = cor(dat), n = nrow(dat)),
      #               indepTest = gaussCItest, ## indep.test: partial correlations
      #               alpha = pc_alpha, labels = labels, u2pd = "retry")
      # G.adjmat = t(as(pc.learn, "matrix"))
      # if(!isValidGraph(t(G.adjmat), type = "cpdag")) {
      #   G.adjmat = as(dag2cpdag(as(G.adjmat, "graphNEL")), "matrix")
      # }
      # G.edgelist = edgeList(G.adjmat)
      

      
      # pc algorithm skeleton
      # if(is.null(skel_alpha)) skel_alpha = pc_alpha
      # skeleton.learn = pcalg::skeleton(suffStat = list(C = cor(dat), n = nrow(dat)),
      #                           indepTest = gaussCItest,#method = "original", ## indep.test: partial correlations
      #                           alpha=skel_alpha, labels = labels)
      # skeleton.adjmat = t(as(skeleton.learn, "matrix"))
      # skeleton.edgelist = edgeList(skeleton.adjmat)

    }
    
    
    # Use sparsebn graph
    if(UseSparsebn) {
      #detach(package:igraph)
      if(printProgress) cat("SBN + NL...") 
      
      packages("graph")
      dat.sbn = sparsebnData(x = dat, type = "continuous")
      
      if(matchPC) {
        sbn.learn = estimate.dag(data = dat.sbn, edge.threshold = sum(as(pc.learn, "matrix")))
      } else {
        sbn.learn = estimate.dag(data = dat.sbn)
      }
      if(is.null(sbn_alpha)) sbn_alpha = 0.1
      G.learn = sbn.learn[[select.parameter(x = sbn.learn, data = dat.sbn, alpha = sbn_alpha)]]$edges
      G.learn = to_bn(to_graphNEL(G.learn))
      cpdag.learn = cpdag(G.learn, moral = FALSE)
      cpdag.amat = amat(cpdag.learn)
      cpdag.edgelist = edgeList(cpdag.amat)
      
      skeleton.learn= bnlearn::skeleton(G.learn)
      skeleton.amat = amat(skeleton.learn)
      skeleton.edgelist = edgeList(skeleton.amat)
      
      res_sbn$G_est = list(G.learn = G.learn, cpdag.learn = cpdag.learn, EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
                          skeleton.EdgeList = skeleton.edgelist, skeleton.AdjMatrix = skeleton.amat)
      res_sbn$G_update = func(dat = dat, labels = labels, 
                             AdjMat = cpdag.amat, EdgeList = cpdag.edgelist)
      if(UseSkeleton) {
        res_sbn$skeleton_update = func(dat = dat, labels = labels, 
                                      AdjMat = skeleton.amat, EdgeList = skeleton.edgelist)
      }
      
      # G.amat = as(G.edgelist, "matrix")
      # # convert sparsebn graph to cpDAG using pcalg dag2cpDAG
      # G.amat = as(dag2cpdag(as(G.amat, "graphNEL")), "matrix")
      # G.edgelist = edgeList(G.amat)
      # # sparsebn skeleton
      # 
      
      # skeleton.amat = G.amat
      # for(i in 1:nrow(G.amat)) {
      #   for(j in 1:ncol(G.amat))
      #     if(G.amat[i, j] == 1) {
      #       skeleton.amat[j, i] = 1
      #     }
      # }
      # skeleton.edgelist = edgeList(skeleton.amat)
      # 
      # res_sbn$G_est = list(G.learn = sbn.learn, EdgeList = G.edgelist, AdjMatrix = G.amat, 
      #                      skeleton.EdgeList = skeleton.edgelist, skeleton.AdjMatrix = skeleton.amat)
      # res_sbn$G_update = func(dat = dat, labels = labels, 
      #                         AdjMat = G.amat, EdgeList = G.edgelist)
      # if(UseSkeleton) {
      #   res_sbn$skeleton_update = func(dat = dat, labels = labels, 
      #                                  AdjMat = skeleton.amat, EdgeList = skeleton.edgelist)
      # }
    }

    
    #res1 = list()
    res1$res_pc = res_pc
    res1$res_sbn = res_sbn
    res1$res_emptyG = res_emptyG
    out$res1 = res1
  }
  
  
  
  
  if(V2) {
    # print("Method 2")
    
    # estimate nonlinear edges from empty graph
    #print("updating empty graph")
    
    if(V1){
      res_emptyG = res1$res_emptyG 
    }else{
      G.amat = matrix(0, ncol(dat), ncol(dat))
      edgelist = edgeList(G.amat)
      res_emptyG = func(dat = dat, labels = labels, 
                        AdjMat = G.amat, EdgeList = G.edgelist)
    }

    G.amat.nl = res_emptyG$AdjMatrixAddNL
    G.edgeList.nl = edgeList(G.amat.nl)
    
    
    # apply PC algorithm
    if(UsePC) {
    
      if(printProgress) cat("NL + PC...")      
      
      # add nonlinear to fixed edge set
      whitelist = matrix(labels[which(G.amat.nl==1, arr.ind = TRUE)], ncol = 2)
      #colnames(whitelist) = c("from", "to")
      G.learn = pc.stable(as.data.frame(dat), alpha = pc_alpha, 
                          whitelist = whitelist, max.sx = pc.maxdegree)
      cpdag.learn = cpdag(G.learn, moral = FALSE, wlbl = TRUE)
      # fixedEdges = G.amat.nl
      # fixedEdges[t(fixedEdges)==1] = 1
      # 
      
      
      # pc.learn = pc(suffStat = list(C = cor(dat), n = nrow(dat)),
      #               indepTest = gaussCItest, ## indep.test: partial correlations
      #               alpha = pc_alpha, labels = labels, u2pd = "retry", fixedEdges = fixedEdges)
      
      cpdag.amat = amat(cpdag.learn)
      cpdag.edgelist = edgeList(cpdag.amat)
      
      res_pc$G_update = list(cpdagBN = cpdag.learn, cpdagAddBN = cpdag.learn,
                             EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
                             EdgeListNL = G.edgeList.nl, AdjMatrixNL = G.amat.nl,
                             EdgeListAdd = cpdag.edgelist, AdjMatrixAdd = cpdag.amat, 
                             EdgeListAddNL = G.edgeList.nl, AdjMatrixAddNL = G.amat.nl)
    }
    
    # apply sbn algorithm
    if(UseSparsebn) {
      #detach(package:igraph)
      if(printProgress) cat("NL + SBN...Done! \n")      
      packages("graph")
      
      # add nonlinear to whitelist
      dat.sbn = sparsebnData(x = dat, type = "continuous")
      whitelist = matrix(labels[which(G.amat.nl==1, arr.ind = TRUE)], ncol = 2)
      if(length(whitelist)==0) whitelist = NULL
      
      if(matchPC) {
        sbn.learn = estimate.dag(data = dat.sbn, edge.threshold = sum(as(pc.learn, "matrix")), whitelist = whitelist)
      } else {
        sbn.learn = estimate.dag(data = dat.sbn, whitelist = whitelist)
      }
      
      if(is.null(sbn_alpha)) sbn_alpha = 0.1
      
      G.edgelist = sbn.learn[[select.parameter(x = sbn.learn, data = dat.sbn, alpha = sbn_alpha)]]$edges
      #G.adjmat = as(G.edgelist, "matrix")
      G.learn = to_bn(to_graphNEL(G.edgelist))
      G.learn$learning$whitelist = whitelist
      cpdag.learn = cpdag(G.learn, moral = FALSE, wlbl = TRUE)
      cpdag.amat = amat(cpdag.learn)
      cpdag.edgelist = edgeList(cpdag.amat)
      
      
     
      #G.adjmat = G.adjmat
      #G.edgelist = edgeList(G.adjmat)
      
      res_sbn$G_update = list(cpdagBN = cpdag.learn, cpdagAddBN = cpdag.learn,
                              EdgeList = cpdag.edgelist, AdjMatrix = cpdag.amat, 
                              EdgeListNL = G.edgeList.nl, AdjMatrixNL = G.amat.nl,
                              EdgeListAdd = cpdag.edgelist, AdjMatrixAdd = cpdag.amat, 
                              EdgeListAddNL = G.edgeList.nl, AdjMatrixAddNL = G.amat.nl)
    }
    
    #res2 = list()
    res2$res_pc = res_pc
    res2$res_sbn = res_sbn
    res2$res_emptyG = res_emptyG
    out$res2 = res2
  }

  return(out)
  
}


