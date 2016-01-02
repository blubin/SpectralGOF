##################################
##################################
#                                # 
#      FUNCTIONS FOR SGOF        #
#                                #
##################################
##################################


require(sna)
SGOF<-function(fittedERGM=NULL,networkModel=NULL, observedNetwork=NULL, observedSpectrum=NULL, nullDistribution=NULL,nsim=250){
  ##
  # preliminaries
  ##
  if(!is.null(fittedERGM)) require(ergm)
  if(!is.null(networkModel)){
    if(!is.function(networkModel)) stop('networkModel must be supplied as a function')
  }
  if(is.null(fittedERGM)&is.null(networkModel)){
    stop('exactly one of fittedERGM or networkModel must be supplied as an argument')
  }
  if(!is.null(networkModel)){
    if(is.null(observedNetwork)&is.null(observedSpectrum)){
      stop('with networkModel, exactly one of observedNetwork or observedSpectrum must be supplied as an argument')
    }
  }
  if(is.null(observedSpectrum)){
    if(!is.null(fittedERGM)){
      n<-dim(fittedERGM$network[,])[1]
      S<-eigen(diag(rowSums(fittedERGM$network[,]))-fittedERGM$network[,], only.values=TRUE)$values
      S<-S/sum(S)
    }else{
      if(dim(observedNetwork)[1]!=dim(observedNetwork)[2]) stop('observedNetwork must be a square adjacency matrix')
      if(!all(observedNetwork==t(observedNetwork))) warning('observedNetwork must be symmetric (undirected). What follows may be meaningless')
      if(any(diag(observedNetwork)!=0)) warning('observedNetwork has loops. What follows may be meaningless')
      S<-eigen(diag(rowSums(observedNetwork))-observedNetwork, only.values=TRUE)$values
      S<-S/sum(S)
      n<-length(S)
    }
  }else{
    S<-observedSpectrum
    if(any(is.complex(S))) warning('supplied spectrum contains complex eigenvalues. What follows may be meaningless.')
    n<-length(S)
  }
  if(!is.null(nullDistribution)){
    if(dim(nullDistribution)[2]!=n) stop('the null distribution supplied implies the wrong number of nodes')
    if(dim(nullDistribution)[1]<100) warning('fewer than 100 simulated spectra from the null distribution may underestimate the true dispersion of the null')
    NDSS<-nullDistribution
    NDSS<-NDSS/rowSums(NDSS)
  }
  if(!is.null(fittedERGM)&is.null(nullDistribution)){
    print("fitting null model")
    fit<-ergm(fittedERGM$network~edges)
    print("simulating from null model")
    simfit<-list()
    i<-1
    while(length(simfit)<nsim){
      simfit[[i]]<-simulate(fit)
      i<-i+1
    }
    NDSS<-lapply(simfit[], FUN=function(net){
      sim<-net[,]
      eigen(diag(rowSums(sim))-sim, only.values=TRUE)$values
    })
    NDSS<-t(array(unlist(NDSS),c(n,length(simfit))))
    NDSS<-NDSS/rowSums(NDSS)
  }
  if(is.null(fittedERGM)&is.null(nullDistribution)){
    e<-sum(observedNetwork)/2
    p<-e/(n*(n-1)/2)
    print("simulating from null model")
    simfit<-list()
    for(i in 1:nsim){
      simfit[[i]]<-rgraph(n=n, m=1, tprob=p, mode="graph", diag=FALSE, replace=FALSE,tielist=NULL, return.as.edgelist=FALSE)
    }
    NDSS<-lapply(simfit[], FUN=function(net){
      sim<-net[,]
      eigen(diag(rowSums(sim))-sim, only.values=TRUE)$values
    })
    NDSS<-t(array(unlist(NDSS),c(n,length(simfit))))
    NDSS<-NDSS/rowSums(NDSS)
  }
  
  ##
  # simulate from fitted model
  ##
  if(!is.null(fittedERGM)){
    print("simulating from fitted ergm")
    simfit<-list()
    for(i in 1:nsim){
      simfit[[i]]<-simulate(fittedERGM)
    }
    
    DSS<-lapply(simfit[], FUN=function(net){
      sim<-net[,]
      eigen(diag(rowSums(sim))-sim, only.values=TRUE)$values
    })
    DSS<-t(array(unlist(DSS),c(n,length(simfit))))
    DSS<-DSS/rowSums(DSS)
  }
  if(!is.null(networkModel)){
    print("simulating from network model")
    simfit<-list()
    for(i in 1:nsim){
      simfit[[i]]<-networkModel()[,]
    }
    DSS<-lapply(simfit[], FUN=function(net){
      sim<-net[,]
      eigen(diag(rowSums(sim))-sim, only.values=TRUE)$values
    })
    DSS<-t(array(unlist(DSS),c(n,length(simfit))))
    DSS<-DSS/rowSums(DSS)
    
  }
    
  eucNum<-mean(apply(DSS, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) ) 
  
  eucCI<-quantile(apply(DSS, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }),probs=c(.95,.05) ) 
  
  eucDen<-mean(apply(NDSS, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) ) 
  
  SGOFeuc<-1-(eucNum/eucDen)
  print(paste(round(SGOFeuc,digits=3)," (5th, 95th quantiles:",round(1-(eucCI[1]/eucDen), digits=2),", ",round(1-(eucCI[2]/eucDen), digits=2),")",sep=""))  
  output<-list(observedSpectrum=S, 
               allSimSpectra=DSS, 
               allNullSpectra=NDSS, 
               SGOF=SGOFeuc, 
               SGOF5=1-(eucCI[1]/eucDen),
               SGOF95=1-(eucCI[2]/eucDen))
  class(output)<-"SGOF"
  return(output)
}

print.SGOF<-function(sgofObject){
  print(paste(round(sgofObject$SGOF,digits=3)," (",round(sgofObject$SGOF5,digits=3),", ",round(sgofObject$SGOF95,digits=3),")", sep=""))
}

SGOFse<-function(SGOFobject){
  S<-SGOFobject$observedSpectrum
  eucNum<-(apply(SGOFobject$allSimSpectra, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) ) 
  
  eucDen<-mean(apply(SGOFobject$allNullSpectra, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) )
  
  return(sd(1-(eucNum/eucDen))/sqrt(length(eucNum)))
}

SGOFdens<-function(SGOFobject){
  S<-SGOFobject$observedSpectrum
  eucNum<-(apply(SGOFobject$allSimSpectra, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) ) 
  
  eucDen<-mean(apply(SGOFobject$allNullSpectra, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) )
  
  plot(density(1-(eucNum/eucDen)), main="Distribution of SGOF")
  abline(v=SGOFobject$SGOF, col="red")
  
  text(SGOFobject$SGOF, 0.05, label="mean SGOF", col="red", pos=4, cex=.7)
}

summarySGOF<-function(sgofobject){
  cat("****************************************\n")
  cat("Spectral goodness of fit (SGOF) analysis\n")
  cat("****************************************\n\n")
  cat("mean SGOF of fitted model:", round(sgofobject$SGOF, digits=3),"\n") 
  cat("     5th percentile:", round(sgofobject$SGOF5, digits=3),"\n")
  cat("     95th percentile:", round(sgofobject$SGOF95, digits=3),"\n\n")
  cat("standard error of the mean:", round(SGOFse(sgofobject),digits=4),"\n\n")
  cat("Null model SGOF percentiles:\n")
  NMP<-quantile(apply(sgofobject$allNullSpectra, 1, FUN=function(vec,S=sgofobject$observedSpectrum){   
    sqrt(sum((vec-S)^2))  
  }),probs=c(.95,.05) )
  eucDen<-mean(apply(sgofobject$allNullSpectra, 1, FUN=function(vec,S=sgofobject$observedSpectrum){   
    sqrt(sum((vec-S)^2))  
  }) ) 
  cat("     5th percentile:", round(1-(NMP[1]/eucDen), digits=3),"\n")
  cat("     95th percentile:", round(1-(NMP[2]/eucDen), digits=3),"\n\n")
  cat("\n")
  cat("****************************************\n\n")
}


plotSGOFerrors<-function(sgofobject,includeLegend=TRUE, style="flat", 
                         logepsilon=1e-6, ptcex=1, insetForFlat=TRUE, insetLoc=c(0.3,0.7,0.025,0.425)){
  DSS<-sgofobject$allSimSpectra
  NDSS<-sgofobject$allNullSpectra
  S<-sgofobject$observedSpectrum
  n<-length(S)
  mse<-apply( apply(DSS,1,FUN=function(vec){(vec-S)^2}) ,1,mean)
  nmse<-apply( apply(NDSS,1,FUN=function(vec){(vec-S)^2}) ,1,mean)
  FMD<-apply(DSS,1,FUN=function(vec){sqrt(sum(vec-S)^2)})
  FMS<-DSS[ which(abs(FMD-mean(FMD))==min(abs(FMD-mean(FMD)))),]
  NMD<-apply(NDSS,1,FUN=function(vec){sqrt(sum(vec-S)^2)})
  NMS<-NDSS[ which(abs(NMD-mean(NMD))==min(abs(NMD-mean(NMD)))),]
  
  samesign<-which(sign(S-FMS)==sign(S-NMS))
  diffsign<-setdiff(1:n,samesign)
  
  better<-which(abs(NMS-S)[samesign]>abs(FMS-S)[samesign])
  worse<-which(abs(NMS-S)[samesign]<abs(FMS-S)[samesign])
  if(style=="original"){
    plot(S, ylim=c(0,max(c(S,FMS,NMS) )), cex=.1 , pch=20, ylab="Normalized Laplacian Eigenvalue", axes=F)
    box()
    axis(side=1,at=c(0,n), labels=expression(italic(lambda[n]),italic(lambda[0])))
    axis(side=2,labels=TRUE)
    segments(x0=c(1:n)[samesign][better], y0=FMS[samesign][better],y1=S[samesign][better], col="#254CA8",lwd=.7)
    segments(x0=c(1:n)[samesign][worse], y0=NMS[samesign][worse],y1=S[samesign][worse], col="#254CA8",lwd=.7)
    segments(x0=c(1:n)[samesign], y0=NMS[samesign],y1=FMS[samesign], col=ifelse(abs(NMS-S)[samesign]>abs(FMS-S)[samesign], "#23CB1A","#F92025"),lwd=.7) #
    
    segments(x0=c(1:n)[diffsign], y0=NMS[diffsign],y1=S[diffsign], col="#23CB1A",lwd=.7)
    segments(x0=c(1:n)[diffsign], y0=FMS[diffsign],y1=S[diffsign], col="#F92025",lwd=.7)
    points(NMS,cex=.3*ptcex,pch=20, col=rgb(.2,.2,.2,.4))
    points(S,cex=.4*ptcex,pch=20,col=rgb(0,0,0,.6))
    points(FMS,cex=.4*ptcex,pch=20,col="#FBB220")
    if(includeLegend){
      location<-ifelse(sum(abs(diff(S[1:floor(n/3)])))>=sum(abs(diff(S[ceiling(n/3):n]))), "topright","bottomleft")
      op <- par(family = "serif")
      legend(location, legend=c("Remaining Errors", "Explained Errors", "New Errors", "Observed Spectrum","Fitted Spectrum", "Null Spectrum"), 
             pt.bg=c("#254CA8","#23CB1A","#F92025", NA, NA,NA), pch=c(20), fill=c("#254CA8","#23CB1A","#F92025", NA, NA,NA),
             col=c("#254CA8","#23CB1A","#F92025",rgb(0,0,0),"#FBB220", rgb(.2,.2,.2,.4)), pt.cex=c(0,0,0,.8, .8, .6), 
             border=c("white", "white","white", "white", "white" ,"white"), bty="n" )
      par(op)
    }
  }
  if(style=="flat"){
    plot(S-S, ylim=c(min(c(FMS-S,NMS-S)),max(c(FMS-S,NMS-S) )), cex=.1 , pch=20, ylab="Normalized Laplacian Eigenvalue minus observed spectrum", axes=F)
    box()
    axis(side=1,at=c(0,n), labels=expression(italic(lambda[n]),italic(lambda[0])))
    axis(side=2,labels=TRUE)
    segments(x0=c(1:n)[samesign][better], y0=(FMS-S)[samesign][better],y1=(S-S)[samesign][better], col="#254CA8",lwd=.7)
    segments(x0=c(1:n)[samesign][worse], y0=(NMS-S)[samesign][worse],y1=(S-S)[samesign][worse], col="#254CA8",lwd=.7)
    segments(x0=c(1:n)[samesign], y0=(NMS-S)[samesign],y1=(FMS-S)[samesign], col=ifelse(abs(NMS-S)[samesign]>abs(FMS-S)[samesign], "#23CB1A","#F92025"),lwd=.7) #
    
    segments(x0=c(1:n)[diffsign], y0=(NMS-S)[diffsign],y1=(S-S)[diffsign], col="#23CB1A",lwd=.7)
    segments(x0=c(1:n)[diffsign], y0=(FMS-S)[diffsign],y1=(S-S)[diffsign], col="#F92025",lwd=.7)
    points(NMS-S,cex=.3*ptcex,pch=20, col=rgb(.2,.2,.2,.4))
    points(S-S,cex=.4*ptcex,pch=20,col=rgb(0,0,0,.6))
    points(FMS-S,cex=.4*ptcex,pch=20,col="#FBB220")
    if(includeLegend){
      location<-"bottomright"
      op <- par(family = "serif")
      legend(location, legend=c("Remaining Errors", "Explained Errors", "New Errors", "Observed Spectrum","Fitted Spectrum", "Null Spectrum"), 
             pt.bg=c("#254CA8","#23CB1A","#F92025", NA, NA,NA), pch=c(20), fill=c("#254CA8","#23CB1A","#F92025", NA, NA,NA),
             col=c("#254CA8","#23CB1A","#F92025",rgb(0,0,0),"#FBB220", rgb(.2,.2,.2,.4)), pt.cex=c(0,0,0,.8*ptcex, .8*ptcex, .6*ptcex), 
             border=c("white", "white","white", "white", "white" ,"white"), bty="n" )
      par(op)
    }
  }
  if(style=="logflat"){
    logify <- function(x) {
      return(sign(x)*ifelse(abs(x)<logepsilon,0,log(abs(x))-log(logepsilon)))
    }
    plot(logify(S-S), ylim=c(min(c(logify(FMS-S),logify(NMS-S))),max(c(logify(FMS-S),logify(NMS-S)))), cex=.1 , pch=20, ylab=expression(paste(sgn %.% log, " of Absolute Normalized Laplacian Eigenvalue less observed spectrum")), axes=F)
    box()
    axis(side=1,at=c(0,n), labels=expression(italic(lambda[n]),italic(lambda[0])))
    axis(side=2,labels=TRUE)
    segments(x0=c(1:n)[samesign][better], y0=logify((FMS-S)[samesign][better]),y1=logify((S-S)[samesign][better]), col="#254CA8",lwd=.7)
    segments(x0=c(1:n)[samesign][worse], y0=logify((NMS-S)[samesign][worse]),y1=logify((S-S)[samesign][worse]), col="#254CA8",lwd=.7)
    segments(x0=c(1:n)[samesign], y0=logify((NMS-S)[samesign]),y1=logify((FMS-S)[samesign]), col=ifelse(abs(NMS-S)[samesign]>abs(FMS-S)[samesign], "#23CB1A","#F92025"),lwd=.7) #
    
    segments(x0=c(1:n)[diffsign], y0=logify((NMS-S)[diffsign]),y1=logify((S-S)[diffsign]), col="#23CB1A",lwd=.7)
    segments(x0=c(1:n)[diffsign], y0=logify((FMS-S)[diffsign]),y1=logify((S-S)[diffsign]), col="#F92025",lwd=.7)
    points(logify(NMS-S),cex=.3*ptcex,pch=20, col=rgb(.2,.2,.2,.4))
    points(logify(S-S),cex=.4*ptcex,pch=20,col=rgb(0,0,0,.6))
    points(logify(FMS-S),cex=.4*ptcex,pch=20,col="#FBB220")
    if(includeLegend){
      location<-"bottomright"
      op <- par(family = "serif")
      legend(location, legend=c("Remaining Errors", "Explained Errors", "New Errors", "Observed Spectrum","Fitted Spectrum", "Null Spectrum"), 
             pt.bg=c("#254CA8","#23CB1A","#F92025", NA, NA,NA), pch=c(20), fill=c("#254CA8","#23CB1A","#F92025", NA, NA,NA),
             col=c("#254CA8","#23CB1A","#F92025",rgb(0,0,0),"#FBB220", rgb(.2,.2,.2,.4)), pt.cex=c(0,0,0,.8*ptcex, .8*ptcex, .6*ptcex), 
             border=c("white", "white","white", "white", "white" ,"white"), bty="n" )
      par(op)
    }
  }
  if((style=="flat" || style=="logflat") && insetForFlat) {
    par(fig=insetLoc, new=T)
    plot(S, ylim=c(0,max(c(S,FMS,NMS) )), cex=.1 , pch=20, axes=F, xlab=NA, ylab=NA )
    box()
    mtext("Observed Spectrum",3,line=-1, las=0)
    labs = seq(0, max(c(S,FMS,NMS)), length.out=3)
    axis(side=2, tick=TRUE, las=2, at=labs, labels=round(labs,2))
    par(fig=c(0,1,0,1))
  }
}

SGOFnull<-function(observedNetwork, nsim=250){
  if(dim(observedNetwork)[1]!=dim(observedNetwork)[2]) stop('Square adjacency matrix required')
  n<-dim(observedNetwork)[1]
  e<-sum(observedNetwork)/2
  p<-e/(n*(n-1)/2)
  print("simulating from null model")
  simfit<-list()
  for(i in 1:nsim){
    simfit[[i]]<-rgraph(n=n, m=1, tprob=p, mode="graph", diag=FALSE, replace=FALSE,tielist=NULL, return.as.edgelist=FALSE)
  }
  NDSS<-lapply(simfit[], FUN=function(net){
    sim<-net[,]
    eigen(diag(rowSums(sim))-sim, only.values=TRUE)$values
  })
  NDSS<-t(array(unlist(NDSS),c(n,length(simfit))))
  NDSS<-NDSS/rowSums(NDSS)
  S<-eigen(diag(rowSums(observedNetwork))-observedNetwork, only.values=TRUE)$values
  S<-S/sum(S)
  
  eucCI<-quantile(apply(NDSS, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }),probs=c(.05,.95) ) 
  
  eucDen<-mean(apply(NDSS, 1, FUN=function(vec){   
    sqrt(sum((vec-S)^2))  
  }) ) 
  
  output<-list(
    nullSGOF=0,
    allNullSpectra=NDSS, 
    nullSGOF5=as.numeric(1-(eucCI[2]/eucDen)),
    nullSGOF95=as.numeric(1-(eucCI[1]/eucDen))
     )
  class(output)<-"SGOFnull"
  return(output
    )
 }
print.SGOFnull<-function(sgofObject){
  print(paste(round(sgofObject$nullSGOF,digits=3)," (",round(sgofObject$nullSGOF5,digits=3),", ",round(sgofObject$nullSGOF95,digits=3),")", sep=""))
}
