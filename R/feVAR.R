#####  PREPROCESSING  #####

# apply box-cox transformation (auxiliary)
makeBoxCox <- function(x, par, a=NULL) {
  if(is.null(a)) {
    if(sum(x<=0,na.rm=T)==0) a <- 0 else a <- -min(x,na.rm=T)
    if(sum(x+a==0,na.rm=T)>0) a <- a+min((x+a)[which(x+a>0)],na.rm=T)/2
    }
  if(is.null(par)) {
    #f <- function(lam) {
    #  z <- ((x+a)^lam-1)/lam
    #  sum(dnorm(z,mean(z),sd(z),log=T),na.rm=T)
    #  }
    #par <- optimize(f, interval=c(-2,2), maximum=T)$maximum
    bc <- MASS::boxcox(lm(I(x+a)~1, y=T), lambda=seq(-2,2,0.001), plotit=F)
    par <- bc$x[which.max(bc$y)]
    z <- ((x+a)^par-1)/par
    } else if(par==1) {
    a <- 0
    z <- x
    } else if(par==0) {
    z <- log(x+a)  
    } else {
    z <- ((x+a)^par-1)/par
    }
  attr(z,"offset") <- a
  attr(z,"lambda") <- par
  z
  }

# invert box-cox transformation (auxiliary)
invBoxCox <- function(z, par, a) {
  if(par==1) {
    z
    } else {
    if(par==0) {
      exp(z)-a  
      } else {
      (z*par+1)^(1/par)-a
      }
    }
  }

# linear interpolation (auxiliary)
linInterp <- function(x, xpre=F, xpost=F, xdef=NA) {
  if(sum(!is.na(x))>=2) {
    x <- approx(x, xout=1:length(x))$y
    }
  if(sum(!is.na(x))==0) {
    x <- rep(xdef,length(x))
    } else {
    if(xpre) {
      xsx <- min(which(!is.na(x)))
      x[1:xsx] <- x[xsx]
      }
    if(xpost) {
      xdx <- max(which(!is.na(x)))
      x[xdx:length(x)] <- x[xdx]
      }
    }
  x
  }

# lines function with linear interpolation (auxiliary)
lines2 <- function(x, y=NULL, ...) {
  lines(x=linInterp(x), y=linInterp(y), ...)
  }

# fill missing with most recent value (auxiliary)
fillNAs <- function(x, xpre=F, xpost=F) {
  if(sum(is.na(x))>0) {
    n <- length(x)
    if(sum(!is.na(x))>0) {
      xsx <- min(which(!is.na(x)))
      if(xpost) xdx <- n else xdx <- max(which(!is.na(x)))
      for(i in (xsx+1):xdx) {
        if(is.na(x[i])) x[i] <- x[i-1]
        }
      if(xpre) x[1:xsx] <- x[xsx]
      }
    }
  x
  }

# unit root test for one variabile (auxiliary)
oneTest <- function(x, unit, max.lag, trend, check) {
  if(is.null(unit)) {
    x <- linInterp(x, xpre=F, xpost=F, xdef=NA)
    n <- length(na.omit(x))
    } else {
    gr <- levels(factor(unit))
    nvet <- c()
    for(w in gr) {
      ind <- which(unit==w)
      x[ind] <- linInterp(x[ind], xpre=F, xpost=F, xdef=NA)
      nvet[w] <- length(na.omit(x[ind]))
      }
    isOK <- which(!is.na(x))
    x <- x[isOK]
    unit <- unit[isOK]
    n <- min(nvet)
    }
  if(is.null(max.lag)) {
    max.lag <- round(sqrt(n))
    } else if(!is.numeric(max.lag)) {
    #max.lag <- min(n-3,trunc((n-1)^(1/3)))
    max.lag <- round(sqrt(n))
    if(check) warning("Argument 'max.nlags' was determined automatically: ",max.lag,sep="",call.=F)
    } else {
    max.lag <- round(max(max.lag,na.rm=T))
    }
  if(max.lag<0) {
    max.lag <- 0
    if(check) warning("Argument 'max.nlags' was set to the minimum possible value: 0",call.=F)
    }
  if(max.lag>n-5) {
    max.lag <- max(0,n-5)
    if(check) warning("Argument 'max.nlags' was set to the maximum possible value: ",max.lag,sep="",call.=F)
    }
  if(is.null(unit)) {
    res <- res1 <- adfFun(x=x, max.lag=max.lag, trend=trend)
    res2 <- kpssFun(x=x, max.lag=max.lag, trend=trend)
    for(i in 1:3) {
      res[[i]] <- c(adf=res1[[i]],kpss=res2[[i]])
      }
    } else {
    gr <- levels(factor(unit))
    res1 <- res2 <- vector("list",length=3)
    for(w in gr) {
      ind <- which(unit==w)
      iadf <- adfFun(x=x[ind], max.lag=max.lag, trend=trend)
      ikpss <- kpssFun(x=x[ind], max.lag=max.lag, trend=trend)
      for(j in 1:length(res1)) {
        res1[[j]] <- c(res1[[j]],iadf[[j]])
        res2[[j]] <- c(res2[[j]],ikpss[[j]])
        }
      }
    #
    adfComb <- function(p) {
      p[which(p<=0)] <- 1e-4
      p[which(p>=1)] <- 1-1e-4
      ind <- which(!is.na(p))
      if(length(ind)>0) {
        m <- length(ind)
        logp <- qnorm(p[ind])
        rhat <- 1-var(logp,na.rm=T)
        rstar <- max(rhat,-1/(m-1))
        auxz <- sum(logp,na.rm=T)/sqrt(m*(1+(m-1)*(rstar+0.2*sqrt(2/(m+1))*(1-rstar))))
        #auxz <- sum(logp,na.rm=T)/sqrt(m)
        #c(p,'(combined)'=2*pnorm(-abs(auxz)))
        2*pnorm(-abs(auxz))
        } else {
        NaN
        #c(p,'(combined)'=NaN)
        }
      }
    #
    kpssComb <- function(stat) {
      if(trend) {
        mu <- 1/15
        s2 <- 11/6300
        } else {
        mu <- 1/6
        s2 <- 1/45  
        }
      z <- (mean(stat,na.rm=T)-mu)/sqrt(s2/length(stat))
      c(2*pnorm(-abs(z)))
      }
    #
    res1 <- lapply(res1, function(z){names(z)<-gr; z})
    names(res1) <- names(iadf)
    res1$p.value <- adfComb(res1$p.value)
    res2 <- lapply(res2, function(z){names(z)<-gr; z})
    names(res2) <- names(ikpss)
    res2$p.value <- kpssComb(res2$statistic)
    res <- list()
    res$statistic <- cbind(adf=res1$statistic,kpss=res2$statistic)
    res$lag.selected <- cbind(adf=res1$lag.selected,kpss=res2$lag.selected)
    res$p.value <- c(adf=res1$p.value,kpss=res2$p.value)
    }
  res
  }

# function for adf test (auxiliary)  
adfFun <- function(x, max.lag, trend) {
  #
  doADF <- function(k) {
    y <- diff(x)
    n <- length(y)
    k <- k+1
    z <- embed(y,k)
    yt <- z[,1]
    xt1 <- x[k:n]
    tt <- k:n
    if(k>1) {
      yt1 <- z[,2:k,drop=F]
      if(trend) res <- lm(yt~xt1+yt1+tt) else res <- lm(yt~xt1+yt1)
      } else {
      if(trend) res <- lm(yt~xt1+tt) else res <- lm(yt~xt1)
      }
    suppressWarnings(
      res.sum <- summary.lm(res)$coefficients
      )
    if(nrow(res.sum)>=2) {
      STAT <- res.sum[2,1]/res.sum[2,2]
      if(trend) {
        table <- rbind(c(-4.38, -3.95, -3.60, -3.24, -2.14, -1.14, -0.81, -0.50, -0.15),
                       c(-4.16, -3.80, -3.50, -3.18, -2.16, -1.19, -0.87, -0.58, -0.24),
                       c(-4.05, -3.73, -3.45, -3.15, -2.17, -1.22, -0.90, -0.62, -0.28), 
                       c(-3.98, -3.69, -3.42, -3.13, -2.18, -1.23, -0.92, -0.64, -0.31),
                       c(-3.97, -3.67, -3.42, -3.13, -2.18, -1.24, -0.93, -0.65, -0.32),
                       c(-3.96, -3.67, -3.41, -3.13, -2.18, -1.25, -0.94, -0.66, -0.32))
        } else {
        table <- rbind(c(-3.75, -3.33, -2.99, -2.64, -1.53, -0.37,  0.00, 0.34, 0.71),
                       c(-3.59, -3.23, -2.93, -2.60, -1.55, -0.41, -0.04, 0.28, 0.66),
                       c(-3.50, -3.17, -2.90, -2.59, -1.56, -0.42, -0.06, 0.26, 0.63),
                       c(-3.45, -3.14, -2.88, -2.58, -1.56, -0.42, -0.07, 0.24, 0.62),
                       c(-3.44, -3.13, -2.87, -2.57, -1.57, -0.44, -0.07, 0.24, 0.61),
                       c(-3.42, -3.12, -2.86, -2.57, -1.57, -0.44, -0.08, 0.23, 0.6))
        }
      tablen <- dim(table)[2]
      tableT <- c(25, 50, 100, 250, 500, 1e+05)
      tablep <- c(0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99)
      tableipl <- numeric(tablen)
      for(i in (1:tablen)) {
        tableipl[i] <- approx(tableT,table[,i],n,rule=2)$y
        }
      PVAL <- approx(tableipl,tablep,STAT,rule=2)$y
      } else {
      STAT <- PVAL <- NA
      }
    c(STAT,PVAL)
    }
  #
  x <- na.omit(x)
  k <- 0
  if(length(x)>=5 & var(x)>0) {  ## <--
    if(max.lag>0) k <- ar(x,order.max=max.lag)$order
    res <- doADF(k)
    } else {
    res <- c(NaN,NaN)
    }
  list(statistic=res[1], lag.selected=k, p.value=res[2])
  }

# function for kpss test (internal use only)
kpssFun <- function(x, max.lag, trend) {
  #
  doKPSS <- function(lag) {
    n <- length(x)
    if(trend) {
      xt <- 1:n
      e <- residuals.lm(lm(x ~ xt))
      table <- c(0.216, 0.176, 0.146, 0.119, 0)
      } else {
      e <- residuals.lm(lm(x ~ 1))
      table <- c(0.739, 0.574, 0.463, 0.347, 0)
      }
    tablep <- c(0.01, 0.025, 0.05, 0.1, 1)
    s <- cumsum(e)
    eta <- sum(s^2)/(n^2)
    s2 <- sum(e^2)/n
    k <- 0
    for(i in 1:lag) {
      ik <- 0
      for(j in (i+1):n) {
        ik <- ik+e[j]*e[j-i]
        }
      k <- k+(1-i/(lag+1))*ik
      }
    stat <- eta/(s2+2*k/n)
    pval <- approx(table,tablep,stat,rule=2)$y
    c(statistic=stat, p.value=pval)
    }
  #
  x <- na.omit(x)
  k <- 0
  if(length(x)>=5 & var(x)>0) {  ## <--
    if(max.lag>0) k <- ar(x,order.max=max.lag)$order
    res <- doKPSS(k)
    } else {
    res <- c(NaN,NaN)
    }
  list(statistic=unname(res[1]), lag.selected=k, p.value=unname(res[2]))
  }

# perform unit root test
unirootTest <- function(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.nlags=NULL) {
  if(!is.numeric(box.cox)) {
    box.cox <- 1
    warning("Argument 'box.cox' was set to the default value: 1",call.=F)
    }
  if(!is.numeric(ndiff)) {
    ndiff <- 0
    warning("Argument 'ndiff' was set to the default value: 0",call.=F)
    }
  var.names <- var.names[which(!is.na(var.names))]
  dataD <- preProcess(var.names=var.names, unit=unit, time=time, data=data, box.cox=box.cox, ndiff=ndiff, check=T, quiet=T)
  if(is.null(unit)) gr <- NULL else gr <- dataD[,unit]
  if(ndiff>0) trend <- F else trend <- T
  testList <- list()
  for(i in 1:length(var.names)) {
    iadf <- oneTest(x=dataD[,var.names[i]], unit=gr, max.lag=max.nlags, trend=trend, check=(i==1))
    iadf$box.cox <- attr(dataD,"box.cox")[var.names[i],]
    iadf$ndiff <- unname(attr(dataD,"ndiff")[var.names[i]])
    testList[[i]] <- iadf
    }
  names(testList) <- var.names
  class(testList) <- "unirootTest"
  testList
  }

# print method for class 'unirootTest'
print.unirootTest <- function(x, ...) {
  cat("p-values","\n")
  cat("  ADF:  null hypothesis is 'unit root'","\n")
  cat("  KPSS: null hypothesis is 'no unit roots'","\n")
  tab <- sapply(x, function(z){
    pval <- round(z$p.value,4)
    #if(is.matrix(pval)) pval["(combined)",] else pval
    pval
    })
  print(t(tab))
  }

# function for differencing (auxiliary)
diffFun <- function(var.names, data, ndiff, xdef.vet) {
  newdat <- data
  n <- nrow(data)
  for(i in 1:length(var.names)) {
    if(isQuant(data[,var.names[i]])) {
      idat <- linInterp(data[,var.names[i]], xpre=T, xpost=T, xdef=xdef.vet[i])
      if(ndiff[var.names[i]]>0) {
        idat_lag <- c(rep(NA,ndiff[i]),idat)[1:length(idat)]
        newdat[,var.names[i]] <- idat-idat_lag
        } else {
        newdat[,var.names[i]] <- idat
        }
      } else {
      newdat[,var.names[i]] <- fillNAs(data[,var.names[i]], xpre=T, xpost=T)
      }
    }
  attr(newdat,"ndiff") <- ndiff
  newdat
  }

# recognize quantitative variable (auxiliary)
isQuant <- function(x) {
  is.numeric(x)&!setequal(unique(na.omit(x)),c(0,1))
  }

# pre-processing (auxiliary)
preProcess <- function(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.ndiff=2, check=TRUE, quiet=TRUE) {
  if(check) {
    #
    if(missing(data)) stop("Argument 'data' is missing",call.=F)
    if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame",call.=F)
    if(missing(var.names)) stop("Argument 'var.names' is missing",call.=F)
    if(!is.character(var.names)) {
      stop("Argument 'var.names' must be a character vector",call.=F)
      } else {
      var.names <- var.names[which(!is.na(var.names))]
      if(length(var.names)<1) stop("Argument 'var.names' must be a character vector of length 1 or greater",call.=F)
      }
    auxchk <- setdiff(var.names,colnames(data))  
    if(length(auxchk)>0) stop("Unknown variable '",auxchk[1],"' in argument 'var.names'",call.=F)
    #
    unit <- unit[1]
    if(!is.null(unit)&&is.na(unit)) unit <- NULL
    if(!is.null(unit)) {
      if(length(setdiff(unit,colnames(data)))>0) stop("Unknown variable '",unit,"' provided to argument 'unit'",call.=F)
      if(length(intersect(unit,var.names))>0) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'",call.=F)
      if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values",call.=F)
      data[,unit] <- factor(data[,unit])
      }
    #
    time <- time[1]
    if(!is.null(time)&&is.na(time)) time <- NULL
    if(!is.null(time)) {
      if(length(setdiff(time,colnames(data)))>0) stop("Unknown variable '",time,"' provided to argument 'time'",call.=F)
      if(length(intersect(time,var.names))>0) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'",call.=F)
      if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'",call.=F)
      if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'",call.=F)
      if(sum(is.na(data[,time]))>0) stop("Variable '",time,"' provided to argument 'time' contains missing values",call.=F)
      if(is.null(unit)) {
        if(sum(duplicated(data[,time]))>0) stop("Variable '",time,"' contains duplicated values",call.=F)
        } else {
        if(sum(sapply(split(data,data[,unit]),function(x){sum(duplicated(x[,time]))}))>0) stop("Variable '",time,"' contains duplicated values",call.=F)
        }
      }
    }
  #
  nomiQ <- c()
  for(i in 1:length(var.names)) {
    if(isQuant(data[,var.names[i]])) nomiQ <- c(nomiQ,var.names[i])
    }
  if(is.null(box.cox)) {
    if(quiet==F) cat("Automated Box-Cox transformation ... ")
    bcvet <- rep(NA,length(var.names))
    names(bcvet) <- var.names
    for(i in 1:length(nomiQ)) {
      bcvet[nomiQ[i]] <- attr(makeBoxCox(data[,nomiQ[i]], NULL), "lambda")
      }
    box.cox <- bcvet
    if(quiet==F) cat("Done","\n")
    } else {
    if(!is.numeric(box.cox)) {
      box.cox <- 1
      warning("Argument 'box.cox' was set to the default value: 1",call.=F)
      }
    if(length(box.cox)==1&is.null(names(box.cox))) {
      box.cox <- rep(box.cox,length(var.names))
      names(box.cox) <- var.names
      } else if(is.null(names(box.cox))) {
      box.cox <- box.cox[1:length(var.names)]
      names(box.cox) <- var.names
      box.cox[which(is.na(box.cox))] <- 1
      } else {
      box.cox <- box.cox[var.names]
      names(box.cox) <- var.names
      box.cox[which(is.na(box.cox))] <- 1
      }
    if(sum(box.cox< -2)>0) {  
      box.cox[which(box.cox< -2)] <- -2  ## <--
      warning("Argument 'box.cox' was set to the minimum possible value: -2",call.=F)
      }
    if(sum(box.cox>2)>0) {  
      box.cox[which(box.cox>2)] <- 2  ## <--
      warning("Argument 'box.cox' was set to the maximum possible value: 2",call.=F)
      }
    }
  #
  if(is.null(ndiff)) {
    if(quiet==F) cat("Automated differencing ... ")
    maxd <- lagCalc(var.names=var.names, unit=unit, data=data)
    if(is.null(max.ndiff)) {
      max.ndiff <- maxd
      } else {
      max.ndiff <- max(min(round(max.ndiff),maxd),0)
      }
    dvet <- rep(0,length(var.names))
    names(dvet) <- var.names
    for(i in 1:length(nomiQ)) {
      fine <- 0
      while(fine==0) {
        iurt <- unirootTest(var.names=nomiQ[i],unit=unit,time=time,data=data,box.cox=box.cox[nomiQ[i]],ndiff=dvet[i],max.nlags=NULL)
        ipval <- iurt[[1]]$p.value
        #if(is.matrix(ipval)) ipvalOK <- ipval[nrow(ipval),] else ipvalOK <- ipval
        #if(ipvalOK[1]>0.05&ipvalOK[2]<0.05) {
        if(ipval[1]>0.05&ipval[2]<0.05) {
          dvet[i] <- dvet[i]+1
          } else {
          fine <- 1
          }
        if(dvet[i]>=max.ndiff) fine <- 1
        }
      }
    ndiff <- dvet
    if(quiet==F) cat("Done","\n")
    } else {
    if(!is.numeric(ndiff)) {
      ndiff <- 0
      warning("Argument 'ndiff' was set to the default value: 0",call.=F)
      }
    if(length(ndiff)==1&is.null(names(ndiff))) {
      ndiff <- rep(ndiff,length(var.names))
      names(ndiff) <- var.names
      } else if(is.null(names(ndiff))) {
      ndiff <- ndiff[1:length(var.names)]
      names(ndiff) <- var.names
      ndiff[which(is.na(ndiff))] <- 0
      } else {
      ndiff <- ndiff[var.names]
      names(ndiff) <- var.names
      ndiff[which(is.na(ndiff))] <- 0
      }
    ndiff <- round(ndiff)
    if(sum(ndiff<0|ndiff==Inf)>0) {  
      ndiff[which(ndiff<0|ndiff==Inf)] <- 0
      warning("Argument 'ndiff' was set to the minimum possible value: 0",call.=F)
      }
    }
  # sort by time
  if(!is.null(time)) {
    if(is.null(unit)) {
      data <- data[order(data[,time]),]
      } else {
      data <- data[order(data[,unit],data[,time]),]
      }
    }
  # box-cox transformation
  dataL <- data
  aval <- rep(NA,length=length(var.names))
  names(aval) <- var.names
  for(i in 1:length(nomiQ)) {
    ilam <- box.cox[nomiQ[i]]
    if(ilam!=1 & sum(data[,nomiQ[i]]<=0,na.rm=T)>0) {
      #warning("Variable '",nomiQ[i],"' contains 0 or negative values: an offset was added",call.=F)
      }
    iz <- makeBoxCox(data[,nomiQ[i]], ilam)
    aval[nomiQ[i]] <- attr(iz,"offset")
    dataL[,nomiQ[i]] <- iz
    }
  bcOK <- cbind(offset=aval, lambda=box.cox)
  # differencing
  if(is.null(unit)) {
    n <- nrow(data)
    for(i in 1:length(var.names)) {
      if(ndiff[var.names[i]]>n-5) {
        ndiff[var.names[i]] <- 0
        warning("Differencing not applied to variable '",var.names[i],"'",call.=F)
        }
      }
    xdef.vet <- rep(NA,length(var.names))
    names(xdef.vet) <- var.names
    xdef.vet[nomiQ] <- colMeans(dataL[,nomiQ,drop=F],na.rm=T)
    dataD <- diffFun(var.names=var.names, data=dataL, ndiff=ndiff, xdef.vet=xdef.vet)
    attr(dataD,"box.cox") <- bcOK
    attr(dataD,"ndiff") <- ndiff
    obj <- dataD
    #if(max(ndiff)>0) {
    #  obj <- dataD[setdiff(1:nrow(data),1:max(ndiff)),,drop=F]
    #  } else {
    #  obj <- dataD
    #  }
    } else {
    dataD <- dataL
    isNA <- c()
    gr <- levels(data[,unit])
    val0 <- matrix(nrow=length(gr),ncol=length(var.names))
    rownames(val0) <- gr
    colnames(val0) <- var.names
    n_gr <- c()
    xdef.vet <- rep(NA,length(var.names))
    names(xdef.vet) <- var.names
    xdef.vet[nomiQ] <- colMeans(dataL[,nomiQ,drop=F],na.rm=T)
    for(w in 1:length(gr)) {
      ind <- which(data[,unit]==gr[w])
      n_gr[w] <- length(ind)
      if(max(ndiff)>0) isNA <- c(isNA, ind[1]:ind[max(ndiff)])
      dataD[ind,] <- diffFun(var.names=var.names, data=dataL[ind,], ndiff=ndiff, xdef.vet=xdef.vet)
      val0[w,] <- as.numeric(data[ind[1],var.names])
      }
    n <- max(n_gr)
    for(i in 1:length(var.names)) {
      if(ndiff[var.names[i]]>n-5) {
        ndiff[var.names[i]] <- 0
        warning("Differencing not applied to variable '",var.names[i],"'",call.=F)
        }
      }
    attr(dataD,"box.cox") <- bcOK
    attr(dataD,"ndiff") <- ndiff
    obj <- dataD
    #obj <- dataD[setdiff(1:nrow(data),isNA),,drop=F]
    }
  #attr(obj,"box.cox") <- bcOK
  #attr(obj,"ndiff") <- ndiff
  obj
  }


#####  VAR MODELING  #####

# lag function
LAG <- function(x, lag, unit=NULL, del=0) {
  #
  lfun <- function(v) {
    if(lag>0) {
      z <- c(rep(NA,lag),v)[1:length(v)]
      } else {
      z <- v  
      }
    if(del>0) z[1:del] <- NA
    z
    }
  #
  lag <- max(round(lag),0)
  del <- max(round(del),0)
  if(is.null(unit)) {
    lfun(x)
    } else {
    gr <- levels(factor(unit))
    lx <- c()
    for(i in 1:length(gr)) {
      ind <- which(unit==gr[i])
      lx[ind] <- lfun(x[ind])
      }
    lx
    }
  }

# generate matrix of lags (auxiliary)
lagMat <- function(x, nlags, unit=NULL, del=0) {
  nlags <- max(round(nlags),0)
  del <- max(round(del),0)
  if(nlags>0) {
    #
    lfun <- function(v) {
      n <- length(v)
      M <- matrix(nrow=n,ncol=nlags)
      for(i in 1:nlags) {
        M[,i] <- c(rep(NA,i),v)[1:n]
        }
      colnames(M) <- 1:nlags
      if(del>0) M[1:del,] <- NA
      M
      }
    #
    if(is.null(unit)) {
      lfun(x)
      } else {
      M <- matrix(nrow=length(x),ncol=nlags)
      gr <- levels(factor(unit))
      for(i in 1:length(gr)) {
        ind <- which(unit==gr[i])
        M[ind,] <- lfun(x[ind])
        }
      colnames(M) <- 1:nlags
      M
      }
    } else {
    x
    }
  }

# fit a single equation (auxiliary)
oneRegFit <- function(y.name, var.names, nlags, exogenous, data, unit, interc, trnd, del) {
  nomi <- c(y.name,var.names)
  names(nlags) <- nomi
  if(!is.null(exogenous)) {
    estr <- paste0(paste(exogenous,collapse="+"),"+")
    } else {
    estr <- ""
    }
  if(trnd[y.name]==1) {
    tstr <- "`(trend)`+"
    } else if(trnd[y.name]>1) {
    tstr <- paste0(unit,":`(trend)`+")
    } else {
    tstr <- ""
    }
  if(is.null(unit)) {
    xstr <- c()
    for(i in 1:length(nomi)) {
      if(nlags[i]>=1) {
        if(del>0) {
          xstr[i] <- paste("LAG(",nomi[i],",",1:nlags[i],",del=",del,")", sep="",collapse="+")
          } else {
          xstr[i] <- paste("LAG(",nomi[i],",",1:nlags[i],")", sep="",collapse="+")
          }
        } else {
        xstr[i] <- NA
        }
      }
    xstr <- na.omit(xstr)
    if(length(xstr)>0) {
      xstrOK <- paste(xstr, collapse="+")
      } else {
      xstrOK <- "1"
      }
    if(interc[y.name]==0) {
      form <- formula(paste(y.name,"~-1+",tstr,estr,xstrOK,sep=""))
      } else if(interc[y.name]==1) {
      form <- formula(paste(y.name,"~",tstr,estr,xstrOK,sep=""))
      } else {
      form <- formula(paste(y.name,"~-1+",unit,"+",tstr,estr,xstrOK,sep=""))
      }
    } else {
    xstr <- c()
    for(i in 1:length(nomi)) {
      if(nlags[i]>=1) {
        if(del>0) {
          xstr[i] <- paste("LAG(",nomi[i],",",1:nlags[i],",",unit,",del=",del,")", sep="", collapse="+")
          } else {
          xstr[i] <- paste("LAG(",nomi[i],",",1:nlags[i],",",unit,")", sep="",collapse="+")
          }
        } else {
        xstr[i] <- NA
        }
      }
    xstr <- na.omit(xstr)
    if(length(xstr)>0) {
      xstrOK <- paste(xstr, collapse="+")
      } else {
      xstrOK <- "1"  
      }
    if(interc[y.name]==0) {
      form <- formula(paste(y.name,"~-1+",tstr,estr,xstrOK,sep=""))        
      } else if(interc[y.name]==1) {
      form <- formula(paste(y.name,"~",tstr,estr,xstrOK,sep=""))        
      } else {
      form <- formula(paste(y.name,"~-1+",unit,"+",tstr,estr,xstrOK, sep=""))
      }
    }
  mod <- lm(form, data=data, na.action=na.exclude)
  mod$call$formula <- form
  mod$call$data.orig <- NULL
  mod$call$data.used <- NULL
  mod
  }

# fit all equations (auxiliary)
fitEquations <- function(var.names, unit, exogenous, data, nlags, interc, trnd, penaltyFun, del, quiet=quiet) {
  mod <- list()
  for(i in 1:length(var.names)) {
    iy <- var.names[i]
    ix <- setdiff(var.names, iy)
    mod[[i]] <- oneRegFit(y.name=iy, var.names=ix, nlags=rep(nlags,length(var.names)), exogenous=exogenous, data=data, unit=unit, interc=interc, trnd=trnd, del=del)
    }
  names(mod) <- var.names
  mod
  }

# compute maximum lag length (auxiliary)
lagMaxCalc <- function(var.names, unit, exogenous, data, interc, trnd) {
  if(is.null(exogenous)) estr <- "1" else estr <- paste(exogenous,collapse="+")
  form1 <- formula(paste("~",estr,sep=""))
  nex <- ncol(model.matrix(form1, data=data))
  form2 <- formula(paste0("~",paste(c(var.names,exogenous),collapse="+"),sep=""))
  n <- nrow(model.matrix(form2, data=data))
  if(is.null(unit)) ng <- 1 else ng <- length(split(data,data[,unit]))
  floor((n-max(interc)-max(trnd)-nex-5)/(ng+length(var.names)))
  }
 
# compute default lag length (auxiliary)
lagCalc <- function(var.names, unit, data) {
  if(is.null(unit)) {
    nt <- nrow(data)
    min(nt-3,round(sqrt(nt)))
    } else {
    nt <- sapply(split(data, data[,unit]), nrow)
    min(nt-3,round(sqrt(nt)))
    }
  }

# get number of lags per variable (auxiliary)
getLags <- function(model) {
  bvet <- model$coefficients
  xstr <- names(bvet)[grep("^LAG\\(",names(bvet))]
  if(length(xstr)>0) {
    xnam <- xlag <- c()
    for(i in 1:length(xstr)) {
      istr1 <- strsplit(xstr[i],"^LAG\\(")[[1]][2]
      istr2 <- strsplit(istr1,",")[[1]]
      xnam[i] <- istr2[1]
      xlag[i] <- as.numeric(gsub("\\)","",istr2[2]))
      }
    mat <- data.frame(x=xnam,lag=xlag)
    rownames(mat) <- xstr
    mat
    }
  }

# create beta matrices (auxiliary)
createBeta <- function(modList) {
  xnam <- names(modList)
  lagList <- lapply(modList, getLags)
  suppressWarnings(
    nlags <- max(sapply(lagList,function(x){max(x$lag)}))
    )
  if(nlags>0) {
    betaList <- vector("list",length=nlags)
    names(betaList) <- 1:nlags
    for(i in 1:nlags) {
      imat <- matrix(0,nrow=length(xnam),ncol=length(xnam))
      rownames(imat) <- colnames(imat) <- xnam
      betaList[[i]] <- imat
      }
    for(i in 1:length(xnam)) {
      imod <- modList[[xnam[i]]]
      ilag <- lagList[[xnam[i]]]
      if(!is.null(ilag)) {
        istr <- rownames(ilag)
        icoef <- imod$coefficients
        for(j in 1:length(istr)) {
          ijl <- ilag[j,"lag"]
          ijx <- ilag[j,"x"]
          betaList[[ijl]][ijx,xnam[i]] <- icoef[istr[j]]
          }
        }
      }
    } else {
    betaList <- list()
    imat <- matrix(0,nrow=length(xnam),ncol=length(xnam))
    rownames(imat) <- colnames(imat) <- xnam
    betaList[[1]] <- imat
    }
  betaList
  }

# compute fixed effects (auxiliary)
fixeffCalc <- function(regList, data, call) {
  nomi <- names(regList)
  if(is.null(call$unit)) {
    alpha <- vector(length=length(nomi))
    names(alpha) <- nomi
    for(i in 1:length(nomi)) {
      ireg <- regList[[nomi[i]]]
      if(call$ndiff[nomi[i]]==0) {
        alpha[i] <- ireg$coefficients["(Intercept)"]
        } else {
        iB <- ireg$coefficients
        iB <- iB[setdiff(names(iB),"(Intercept)")]
        istr <- colnames(ireg$model)
        iym <- mean(data[,istr[1]],na.rm=T)
        if(length(istr)>1) {
          iX <- model.matrix(ireg)
          iX <- iX[,setdiff(colnames(iX),"(Intercept)")]
          ixm <- colMeans(iX, na.rm=T)
          alpha[i] <- iym-ixm%*%iB
          } else {
          alpha[i] <- iym  
          }
        }
      }
    } else {
    unam <- call$unit
    gr <- levels(data[,unam])
    alpha <- matrix(nrow=length(gr), ncol=length(nomi))
    rownames(alpha) <- gr
    colnames(alpha) <- nomi
    for(i in 1:length(nomi)) {
      ireg <- regList[[nomi[i]]]
      if(call$ndiff[nomi[i]]==0) {
        alpha[,i] <- ireg$coefficients[paste0(unam,gr)]
        } else {
        iB <- ireg$coefficients
        iB <- iB[setdiff(names(iB),"(Intercept)")]
        istr <- colnames(ireg$model)
        iym <- sapply(split(data[,istr[1]], data[,unam]), mean, na.rm=T)
        if(length(istr)>1) {
          iX <- model.matrix(ireg)
          iX <- iX[,setdiff(colnames(iX),"(Intercept)")]
          ixm <- lapply(split(data.frame(iX), data[rownames(iX),unam]), colMeans, na.rm=T)
          ialph <- c()
          for(j in 1:length(iym)) {
            ialph[j] <- iym[j]-ixm[[j]]%*%iB
            }
          alpha[,i] <- ialph
          } else {
          alpha[,i] <- iym
          }
        }
      }
    }
  alpha
  }

# MASTER FUNCTION
feVAR <- function(var.names, unit=NULL, time=NULL, exogenous=NULL, data, max.nlags=NULL, nlags=NULL, trend=c("none","global","unit"), ic=c("bic","aic","aicc","hqic"), box.cox=1, ndiff=0, max.ndiff=2, auto.restrict=FALSE, imputation=TRUE, em.tol=1e-4, em.maxiter=100, quiet=FALSE) {
  #
  if(missing(data)) stop("Argument 'data' is missing",call.=F)
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame",call.=F)
  if(missing(var.names)) stop("Argument 'var.names' is missing",call.=F)
  if(!is.character(var.names)) {
    stop("Argument 'var.names' must be a character vector of length 2 or greater",call.=F)
    } else {
    var.names <- var.names[which(!is.na(var.names))]
    if(length(var.names)<2) stop("Argument 'var.names' must be a character vector of length 2 or greater",call.=F)
    for(i in 1:length(var.names)) {
      if((var.names[i]%in%colnames(data))==F) stop("Unknown variable '",var.names[i],"' in argument 'var.names'",call.=F)
      if(isQuant(data[,var.names[i]])==F) stop("Variable '",var.names[i],"' in argument 'var.names' is not quantitative",call.=F)
      }
    }
  #
  if(!is.null(unit)) {
    unit <- unit[which(!is.na(unit))]
    if(!is.character(unit)|length(unit)!=1) stop("Argument 'unit' must be either NULL or a character vector of length 1",call.=F)
    if(length(unit)>0) {
      if(length(setdiff(unit,colnames(data)))>0) stop("Unknown variable '",unit,"' provided to argument 'unit'",call.=F)
      if(length(intersect(unit,var.names))>0) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'",call.=F)
      if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values",call.=F)
      data[,unit] <- factor(data[,unit])
      } else {
      unit <- NULL
      }
    }
  #
  if(!is.null(time)) {
    time <- time[which(!is.na(time))]
    if(!is.character(time)|length(time)!=1) stop("Argument 'time' must be either NULL or a character vector of length 1",call.=F)
    if(length(time)>0) {
      if(length(setdiff(time,colnames(data)))>0) stop("Unknown variable '",time,"' provided to argument 'time'",call.=F)
      if(length(intersect(time,var.names))>0) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'",call.=F)
      if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'",call.=F)
      if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'",call.=F)
      if(sum(is.na(data[,time]))>0) stop("Variable '",time,"' provided to argument 'time' contains missing values",call.=F)
      if(is.null(unit)) {
        if(sum(duplicated(data[,time]))>0) stop("Variable '",time,"' contains duplicated values",call.=F)
        } else {
        if(sum(sapply(split(data,data[,unit]),function(x){sum(duplicated(x[,time]))}))>0) stop("Variable '",time,"' contains duplicated values",call.=F)
        }
      } else {
      time <- NULL
      }
    }
  #
  if(!is.null(exogenous)) {
    exogenous <- exogenous[which(!is.na(exogenous))]
    if(!is.character(exogenous)|length(exogenous)==0) stop("Argument 'exogenous' must be either a character vector or NULL",call.=F)
    auxch1 <- setdiff(exogenous,colnames(data))
    if(length(auxch1)>0) stop("Unknown variable '",auxch1[1],"' provided to argument 'exogenous'",call.=F)
    auxch2 <- intersect(exogenous,var.names)
    if(length(auxch2)>0) stop("Variable '",auxch2[1],"' appears in both arguments 'var.names' and 'exogenous'",call.=F)
    auxch3 <- intersect(exogenous,unit)
    if(length(auxch3)>0) stop("Variable '",auxch3[1],"' appears in both arguments 'unit' and 'exogenous'",call.=F)
    if(!is.null(time)) {
      auxch4 <- intersect(exogenous,time)
      if(length(auxch4)>0) stop("Variable '",auxch4[1],"' appears in both arguments 'time' and 'exogenous'",call.=F)
      }
    }
  auto.restrict <- auto.restrict[1]
  if(is.na(auto.restrict)||(!is.logical(auto.restrict)|is.null(auto.restrict))) auto.restrict <- FALSE
  imputation <- imputation[1]
  if(is.na(imputation)||(!is.logical(imputation)|is.null(imputation))) imputation <- FALSE
  if(!is.numeric(em.tol)) em.tol <- 1e-4 else em.tol <- round(max(c(0,em.tol),na.rm=T))
  if(is.na(em.tol)|em.tol=="Inf") em.tol <- 1e-4
  if(!is.numeric(em.maxiter)) em.maxiter <- 100 else em.maxiter <- round(max(c(1,em.maxiter),na.rm=T))
  if(is.na(em.maxiter)|em.maxiter=="Inf") em.maxiter <- 100
  quiet <- quiet[1]
  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  #
  if(is.null(unit)) {
    xtrend <- 0:(nrow(data)-1)
    } else {
    gr <- levels(data[,unit])
    xtrend <- c()
    for(i in 1:length(gr)) {
      ind <- which(data[,unit]==gr[i])
      xtrend[ind] <- 0:(length(ind)-1)
      }
    }
  if(is.null(time)) {
    data[,"(time)"] <- xtrend
    time <- "(time)"
    }
  dataOrig <- data
  data <- preProcess(var.names=c(var.names,exogenous), unit=unit, time=time, data=data, box.cox=box.cox, ndiff=ndiff, max.ndiff=max.ndiff, check=F, quiet=quiet)
  interc <- rep(1,length(var.names))
  trend <- intersect(trend,c("none","global","unit"))
  if(length(trend)==0) trend <- "none" else trend <- trend[1]
  if(trend!="none"&is.null(unit)) trend <- "global"
  if(trend=="none") {
    trnd <- rep(0,length(var.names))
    } else if(trend=="unit"&!is.null(unit)) {
    trnd <- rep(nlevels(data[,unit]),length(var.names))
    } else {
    trnd <- rep(1,length(var.names))
    }
  names(interc) <- names(trnd) <- var.names
  if(!is.null(unit)) {
    interc[names(which(attr(data,"ndiff")==0))] <- nlevels(data[,unit])
    if(trend=="unit") interc[names(which(attr(data,"ndiff")==1))] <- nlevels(data[,unit])
    }
  if(trend=="none") {
    interc[names(which(attr(data,"ndiff")>=1))] <- 0
    } else {
    interc[names(which(attr(data,"ndiff")>1))] <- 0
    trnd[names(which(attr(data,"ndiff")>=1))] <- 0
    }
  if(trend!="none") data[,"(trend)"] <- xtrend
  laglim <- lagMaxCalc(var.names=var.names, unit=unit, exogenous=exogenous, data=data, interc=interc, trnd=trnd)
  lag0 <- min(laglim, lagCalc(var.names=c(var.names,exogenous), unit=unit, data=data))
  if(is.null(nlags)) {
    if(is.null(max.nlags)) {
      max.nlags <- lag0
      } else {
      if(!is.numeric(max.nlags)) {
        warning("Argument 'max.nlags' was determined automatically: ",lag0,sep="",call.=F)
        max.nlags <- lag0
        } else {
        max.nlags <- round(max(max.nlags,na.rm=T))
        }
      if(max.nlags>laglim) {
        warning("Argument 'max.nlags' was set to the maximum possible value: ",laglim,sep="",call.=F)
        max.nlags <- laglim
        }
      if(max.nlags<1) {
        if(max.nlags<0) warning("Argument 'nlags' was set to the minimum possible value: 0",call.=F)
        max.nlags <- NULL
        nlags <- 0
        }
      }
    } else {
    if(!is.numeric(nlags)) {
      warning("Argument 'nlags' was set to the maximum possible value: ",laglim,sep="",call.=F)
      nlags <- laglim
      } else {
      nlags <- round(max(nlags,na.rm=T))
      }
    if(nlags>laglim) warning("Argument 'nlags' was set to the maximum possible value: ",laglim,sep="",call.=F)
    if(nlags<0) warning("Argument 'nlags' was set to the minimum possible value: 0",call.=F)
    nlags <- max(0,min(nlags,laglim))
    }
  p <- length(var.names)
  ic <- intersect(tolower(ic),c("bic","aic","aicc","hqic"))
  if(length(ic)==0) ic <- "bic" else ic <- ic[1]
  if(ic=="aic") {
    penaltyFun <- function(npar,n){2*npar}
    } else if(ic=="aicc") {
    penaltyFun <- function(npar,n){2*npar+(n>npar+1)*(2*npar^2+2*npar)/(n-npar-1)}
    } else if(ic=="hqic") {
    penaltyFun <- function(npar,n){2*npar*log(log(n))}
    } else {
    penaltyFun <- function(npar,n){npar*log(n)}
    }
  if(imputation) {
    emImp <- sum(is.na(dataOrig[,var.names]))>0  ## <-- imputare anche le esogene
    } else {
    emImp <- F
    }
  if(is.null(nlags)&(auto.restrict|emImp)) nlags <- max.nlags
  if(emImp==F) {
    if(is.null(nlags)) {
      icval <- c()
      for(i in 1:max.nlags) {
        auxmod <- fitEquations(var.names=var.names, unit=unit, exogenous=exogenous, data=data, nlags=i, interc=interc, trnd=trnd, del=max.nlags, quiet=T)
        ll <- sum(sapply(auxmod, function(x){logLik(x)[1]}))
        npar <- sum(sapply(auxmod, function(x){attr(logLik(x),"df")}))
        icval[i] <- -2*ll[1]+penaltyFun(npar, sum(sapply(auxmod,nobs)))
        if(quiet==F) cat("VAR with ",i," lags. ",toupper(ic),": ",icval[i],sep="","\n")
        }
      lagOK <- which.min(icval)
      if(quiet==F) cat("Number of lags selected: ",lagOK,sep="","\n")
      mod <- fitEquations(var.names=var.names, unit=unit, exogenous=exogenous, data=data, nlags=lagOK, interc=interc, trnd=trnd, del=0, quiet=T)
      } else {
      mod <- fitEquations(var.names=var.names, unit=unit, exogenous=exogenous, data=data, nlags=nlags, interc=interc, trnd=trnd, penaltyFun=penaltyFun, del=0, quiet=quiet)
      }
    modOK <- list(equations=mod)
    modOK$call <- list(var.names=var.names, unit=unit, time=time, exogenous=exogenous, trend=trend, box.cox=attr(data,"box.cox"), ndiff=attr(data,"ndiff"))
    modOK$data.orig <- dataOrig[,c(unit,time,var.names,exogenous)]
    if(sum(trnd)>0) {
      modOK$data.used <- data[,c(unit,time,var.names,exogenous,"(trend)")]
      } else {
      modOK$data.used <- data[,c(unit,time,var.names,exogenous)]
      }
    class(modOK) <- "feVAR"
    } else {
    isNA <- 1*apply(dataOrig[,var.names,drop=F],2,is.na)
    filldat <- dataOrig
    filldat[,var.names] <- apply(dataOrig[,var.names,drop=F],2,function(x){x[which(is.na(x))]<-mean(x,na.rm=T); x})
    ll <- -Inf
    fine <- count <- 0
    while(fine==0) {
      count <- count+1
      auxdat <- preProcess(var.names=c(var.names,exogenous), unit=unit, time=time,
        data=filldat, box.cox=box.cox, ndiff=ndiff, max.ndiff=max.ndiff, check=F, quiet=T)
      auxeq <- fitEquations(var.names=var.names, unit=unit, exogenous=exogenous,
        data=auxdat, nlags=nlags, interc=interc, trnd=trnd, penaltyFun=penaltyFun, del=0, quiet=T)
      auxmod <- list(equations=auxeq)
      auxmod$call <- list(var.names=var.names, unit=unit, time=time, exogenous=exogenous, trend=trend, box.cox=attr(data,"box.cox"), ndiff=attr(data,"ndiff"))
      auxmod$data.orig <- filldat[,c(unit,time,var.names,exogenous)]
      if(sum(trnd)>0) {
        auxmod$data.used <- auxdat[,c(unit,time,var.names,exogenous,"(trend)")]
        } else {
        auxmod$data.used <- auxdat[,c(unit,time,var.names,exogenous)]
        }
      resMat <- na.omit(do.call(cbind, lapply(auxmod$equations, residuals)))
      auxmod$Sigma <- cov(resMat, use="complete.obs")
      #S <- matrix(0,nrow=length(var.names),ncol=length(var.names))
      #rownames(S) <- colnames(S) <- var.names
      #for(i in 1:nrow(resMat)) S <- S+resMat[i,]%*%t(resMat[i,])
      #auxmod$Sigma <- S/nrow(na.omit(resMat))
      class(auxmod) <- "feVAR"
      auxll <- logLik(auxmod)[1]
      if(auxll<ll) {
        fine <- 1
        } else {
        modOK <- auxmod
        if(auxll-ll<em.tol|count>=em.maxiter) fine <- 1
        }
      if(quiet==F) {
        nch <- nchar(ll)-nchar(auxll)
        if(nch>0) catsep <- rep(" ",nch) else catsep <- ""
        cat("\r","EM iteration ",count,". Log likelihood: ",logLik(modOK)[1],catsep,sep="")
        flush.console()
        }
      if(auxll>ll) {
        auxpr <- predict(modOK)$predicted
        nback <- max(modOK$call$ndiff[var.names]+nlags)
        if(nback>0) {
          tval <- sort(unique(modOK$data.orig[,time]))
          tr_ind <- rev(which(modOK$data.orig[,time]%in%tval[setdiff(1:length(tval),1:nback)]))
          auxback <- predict(modOK, n.ahead=nback, subset=tr_ind)$predicted
          auxbackOK <- lapply(auxback,function(x){x[,time]<-tval[nback:1]; x})
          for(i in 1:length(var.names)) {
            iauxpr <- auxpr[[var.names[i]]]
            iauxpr[which(modOK$data.orig[,time]%in%tval[1:nback]),] <- auxbackOK[[var.names[i]]]
            if(is.null(unit)) {
              iauxpr <- iauxpr[order(iauxpr[,time]),]
              } else {
              iauxpr <- iauxpr[order(iauxpr[,unit],iauxpr[,time]),]
              }
            auxpr[[var.names[i]]] <- iauxpr
            }
          }
        for(i in 1:length(var.names)) {
          ind <- which(isNA[,var.names[i]]==1)
          filldat[ind,var.names[i]] <- auxpr[[var.names[i]]][ind,"mean"]
          }
        ll <- auxll
        }
      }
    if(quiet==F) {
      cat("\n")
      if(count<em.maxiter) {
        cat("EM converged after ",count," iterations",sep="","\n")
        } else {
        cat("EM reached the maximum number of iterations","\n")      
        }
      }
    }
  if(auto.restrict) modOK <- restrictFit(modOK, interc=interc, ic=ic, quiet=quiet)
  modOK$intercepts <- fixeffCalc(modOK$equations, data=modOK$data.orig, call=modOK$call)
  modOK$Beta <- createBeta(modOK$equations)
  resMat <- na.omit(do.call(cbind, lapply(modOK$equations, residuals)))
  modOK$Sigma <- cov(resMat, use="complete.obs")
  #S <- matrix(0,nrow=length(var.names),ncol=length(var.names))
  #rownames(S) <- colnames(S) <- var.names
  #for(i in 1:nrow(resMat)) S <- S+resMat[i,]%*%t(resMat[i,])
  #modOK$Sigma <- S/nrow(na.omit(resMat))
  modOK$companion <- compMat(modOK$Beta)
  modOK$eigen.module <- Mod(eigen(modOK$companion)$values)
  modOK
  }
  
# function for automated restriction (auxiliary)
restrictFit <- function(model, interc, ic, quiet) {
  nomi <- model$call$var.names
  for(i in 1:length(nomi)) {
    if(quiet==F) {
      cat("\r","Performing restriction for equation ",i,"/",length(nomi),sep="")
      flush.console()
      }
    imod <- model$equations[[nomi[i]]]
    idat <- imod$model
    if(model$call$trend=="unit"&model$call$ndiff[nomi[i]]==0) {
    #if(trnd[nomi[i]]>1)
      itr <- paste0("+",model$call$unit,":`(trend)`")
      } else {
      itr <- ""
      }
    if(interc[nomi[i]]==0) {
      iform <- formula(paste0(nomi[i],"~-1+.",itr))
      iform0 <- paste0(nomi[i],"~-1")
      } else if(interc[nomi[i]]==1) {
      iform <- formula(paste0(nomi[i],"~.",itr))
      iform0 <- paste0(nomi[i],"~1")
      } else {
      if(is.null(model$call$unit)) {
        iform <- formula(paste0(nomi[i],"~.",itr))
        iform0 <- paste0(nomi[i],"~1")
        } else {
        iform <- formula(paste0(nomi[i],"~-1+.",itr))
        iform0 <- paste0(nomi[i],"~-1+",model$call$unit)
        }
      }
    if(ic=="aic") {
      ipen <- 2
      } else if(ic=="aicc") {
      inP <- attr(logLik(imod),"df")
      iN <- nobs(imod)
      ipen <- 2+(iN>inP+1)*(2*inP+2)/(iN-inP-1)
      } else if(ic=="hqic") {
      ipen <- 2*log(log(nobs(imod)))
      } else {
      ipen <- log(nobs(imod))
      }
    imod <- stepAIC(lm(iform, data=idat), scope=list(lower=formula(iform0)), k=ipen, trace=F, direction="both")
    istr0 <- setdiff(all.vars(imod$call$formula),model$call$unit)[-1]
    if(length(istr0)>0) {
      istr <- gsub("\\(trend\\)","`\\(trend\\)`",istr0)
      if(paste0(model$call$unit,":`(trend)`")%in%rownames(anova(imod))) istr <- c(istr,paste0(model$call$unit,":`(trend)`"))
      iformOK <- formula(paste0(iform0,"+",paste(istr,collapse="+")))
      } else {
      iformOK <- iform0
      }
    imodOK <- lm(iformOK, data=model$data.used, na.action=na.exclude)
    imodOK$call$formula <- iformOK
    model$equations[[i]] <- imodOK
    }
  if(quiet==F) cat("\n")
  model
  }

# print method for class 'feVAR'
print.feVAR <- function(x, ...) {
  cat("Fixed effects VAR on ",length(x$call$var.names)," variables",sep="","\n")
  if(is.null(x$call$unit)) {
    ng <- 1
    nt <- nrow(x$data.used)
    } else {
    ng <- nlevels(x$data.used[,x$call$unit])
    nt <- range(sapply(split(x$data.used,x$data.used[,x$call$unit]),nrow))
    }
  cat("  number of units: ",ng,sep="","\n")
  cat("  number of time points: ",ifelse(diff(nt)==0,nt,paste0(nt,collapse="-")),sep="","\n")
  cat("  maximum number of lags: ",length(x$Beta),sep="","\n")
  }

# coef method for class 'feVAR'
coef.feVAR <- function(object, ...) {
  lapply(object$equations, function(x){x$coefficients})
  }

# summary method for class 'feVAR'
summary.feVAR <- function(object, ...) {
  lapply(object$equations, summary, ...)
  }

# confint method for class 'feVAR'
confint.feVAR <- function(object, parm, level=0.95, ...) {
  lapply(object$equations, function(x){confint.lm(x, level=level)})
  }

# fitted method for class 'feVAR'
fitted.feVAR <- function(object, ...) {
  fit <- lapply(object$equations, function(x){x$fitted.values})
  #nomiObs <- rownames(object$data.used)
  #for(i in 1:length(fit)) {
  #  ifit <- fit[[i]][nomiObs]
  #  names(ifit) <- nomiObs
  #  fit[[i]] <- ifit
  #  }
  fit <- lapply(object$equations, fitted)
  tab <- do.call(cbind, fit)
  nomi <- colnames(tab)
  tabOK <- data.frame(object$data.used[,object$call$unit],
    object$data.used[,object$call$time],tab)
  colnames(tabOK) <- c(object$call$unit,object$call$time,nomi)
  tabOK
  }

# residuals method for class 'feVAR'
residuals.feVAR <- function(object, ...) {
  #res <- lapply(object$equations, function(x){x$residuals})
  #nomiObs <- rownames(object$data.used)
  #for(i in 1:length(res)) {
  #  ires <- res[[i]][nomiObs]
  #  names(ires) <- nomiObs
  #  res[[i]] <- ires
  #  }
  res <- lapply(object$equations, residuals)
  tab <- do.call(cbind, res)
  nomi <- colnames(tab)
  tabOK <- data.frame(object$data.used[,object$call$unit],
    object$data.used[,object$call$time],tab)
  colnames(tabOK) <- c(object$call$unit,object$call$time,nomi)
  tabOK
  }

# cooks.distance method for class 'feVAR'
cooks.distance.feVAR <- function(model, ...) {
  tab <- do.call(cbind, lapply(model$equations, cooks.distance, ...))
  tab <- cbind(tab, '(maximum)'=apply(tab,1,function(x){
    ifelse(sum(!is.na(x))>0,max(x,na.rm=T),NA)
  }))
  tabOK <- data.frame(tab)
  colnames(tabOK) <- c(colnames(tab))
  tabOK
  }

# logLik method for class 'feVAR'
logLik.feVAR <- function(object, ...) {
  obj <- na.omit(do.call(cbind, lapply(object$equations, residuals)))
  k <- ncol(obj)
  S <- object$Sigma
  iS <- solve(S)
  ll <- -k/2*log(2*pi)-1/2*log(det(S))
  for(i in 1:nrow(obj)) ll <- ll-1/2*(t(obj[i,])%*%iS%*%obj[i,])
  npar <- sum(sapply(object$equations,function(x){attr(logLik(x,...),"df")-1}))+k*(k+1)/2
  attr(ll,"df") <- npar
  class(ll) <- "logLik"
  ll
  }

# nobs method for class 'feVAR'
nobs.feVAR <- function(object, ...) {
  min(sapply(object$equations, nobs, ...))
  }

# extractAIC method for class 'feVAR'
extractAIC.feVAR <- function(fit, scale, k=2, ...) {
  ll <- logLik.feVAR(fit, ...)
  npar <- attr(ll,"df")
  c(npar,-2*ll[1]+k*npar)
  }

# AIC method for class 'feVAR'
AIC.feVAR <- function(object, ...) {
  ll <- logLik.feVAR(object, ...)
  -2*ll[1]+2*attr(ll,"df")
  }

# BIC method for class 'feVAR'
BIC.feVAR <- function(object, ...) {
  ll <- logLik.feVAR(object, ...)
  n <- sum(sapply(object$equations,nobs))
  -2*ll[1]+attr(ll,"df")*log(n)
  }

# in-sample prediction (auxiliary)
insamPred <- function(object, unit.id, level) {
  suppressWarnings(
    pred <- lapply(object$equations, function(x){
      pr <- predict(x,interval="prediction",level=level)
      dd <- data.frame(object$data.used[,c(object$call$unit,object$call$time)],pr)
      colnames(dd) <- c(object$call$unit,object$call$time,"mean","lower","upper")
      dd
      })
    )
  if(!is.null(object$call$unit)&!is.null(unit.id)) {
    predList <- lapply(pred,function(x){
      x[which(x[,object$call$unit]%in%unit.id),]
      })
    } else {
    predList <- pred
    }
  names(predList) <- names(pred)
  predList
  }

# predict method for class 'feVAR'
predict.feVAR <- function(object, n.ahead=0, unit.id=NULL, newdata=NULL, level=0.95, subset=NULL, ...) {
  if(!is.numeric(level)) level <- 0.95 else level <- min(c(1,max(c(level,0),na.rm=T)),na.rm=T)
  if(!is.numeric(n.ahead)) n.ahead <- 0 else n.ahead <- round(max(n.ahead,na.rm=T))
  if(is.na(n.ahead)|n.ahead<0|abs(n.ahead)=="Inf") n.ahead <- 0
  if(!is.null(object$call$unit)) {
    glev <- levels(object$data.used[,object$call$unit])
    if(is.null(unit.id)) {
      unit.id <- glev
      obsdat <- object$data.orig[,c(object$call$unit,object$call$time,object$call$var.names)]
      } else {
      if(is.numeric(unit.id)) {
        unit.id <- intersect(unit.id,1:length(glev))
        if(length(unit.id)==0) {
          unit.id <- 1:length(glev)
          warning("Invalid argument 'unit.id': all units have been considered",call.=F)
          }
        unit.id <- glev[unit.id]
        } else {
        unit.id <- intersect(unit.id, glev)
        if(length(unit.id)==0) {
          unit.id <- glev
          warning("Invalid argument 'unit.id': all units have been considered",call.=F)
          }
        }
      obsdat <- object$data.orig[which(object$data.orig[,object$call$unit]%in%unit.id),c(object$call$unit,object$call$time,object$call$var.names)]
      }
    } else {
    unit.id <- NULL
    obsdat <- object$data.orig[,c(object$call$time,object$call$var.names)]
    }
  if(!is.numeric(subset)|is.null(subset)) subset <- 1:nrow(object$data.used)
  if(n.ahead==0) {
    predList <- insamPred(object, unit.id=unit.id, level=level)
    } else {
    nomi <- object$call$var.names
    if(is.null(object$call$unit)) {
      dat0 <- object$data.used[subset,]
      t_ax <- max(dat0[,object$call$time])+1:n.ahead
      obs0 <- dat0[nrow(dat0),]
      obs0[,nomi] <- NA
      fit <- sx <- dx <- data.frame(matrix(nrow=n.ahead,ncol=length(nomi)))
      colnames(fit) <- colnames(sx) <- colnames(dx) <- nomi
      for(j in 1:n.ahead) {
        dat0 <- rbind(dat0, obs0)
        if("(trend)"%in%colnames(dat0)) dat0[nrow(dat0),"(trend)"] <- dat0[nrow(dat0)-1,"(trend)"]+1
        #
        ## <-- prendere variabili esogene da 'newdata'
        #
        pred <- lapply(object$equations, predict, newdata=dat0, interval="prediction", level=level)
        fit[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),1]})
        sx[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),2]})
        dx[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),3]})
        dat0[nrow(dat0),nomi] <- fit[j,nomi]
        }
      predList <- list()
      for(i in 1:length(nomi)) {
        itab <- data.frame(t_ax,fit[,nomi[i]],sx[,nomi[i]],dx[,nomi[i]])
        colnames(itab) <- c(object$call$time,"mean","lower","upper")
        rownames(itab) <- NULL
        predList[[i]] <- itab
        }
      names(predList) <- nomi
      } else {
      D0 <- object$data.used[subset,]
      dataList <- split(D0, D0[,object$call$unit])[unit.id]
      t_ax <- lapply(dataList, function(x){
        max(x[,object$call$time])+1:n.ahead
        })
      fitList <- sxList <- dxList <- list()
      for(i in 1:length(dataList)) {
        dat0 <- dataList[[i]]
        obs0 <- dat0[nrow(dat0),]
        obs0[,nomi] <- NA
        fit <- sx <- dx <- data.frame(matrix(nrow=n.ahead,ncol=length(nomi)))
        colnames(fit) <- colnames(sx) <- colnames(dx) <- nomi
        for(j in 1:n.ahead) {
          dat0 <- rbind(dat0, obs0)
          if("(trend)"%in%colnames(dat0)) dat0[nrow(dat0),"(trend)"] <- dat0[nrow(dat0)-1,"(trend)"]+1
          #
          ## <-- prendere variabili esogene da 'newdata'
          #
          pred <- lapply(object$equations, predict, newdata=dat0, interval="prediction", level=level)
          fit[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),1]})
          sx[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),2]})
          dx[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),3]})
          dat0[nrow(dat0),nomi] <- fit[j,nomi]
          }
        fitList[[i]] <- fit
        sxList[[i]] <- sx
        dxList[[i]] <- dx
        }
      fitOK <- do.call(rbind,fitList)
      sxOK <- do.call(rbind,sxList)
      dxOK <- do.call(rbind,dxList)
      colnames(fitOK) <- colnames(sxOK) <- colnames(dxOK) <- nomi
      predList <- list()
      for(i in 1:length(nomi)) {
        itab <- data.frame(
          rep(names(dataList),each=n.ahead),
          do.call(c,t_ax),
          fitOK[,nomi[i]],sxOK[,nomi[i]],dxOK[,nomi[i]]
          )
        itab[,1] <- factor(itab[,1],levels=glev)
        colnames(itab) <- c(object$call$unit,object$call$time,"mean","lower","upper")
        rownames(itab) <- NULL
        predList[[i]] <- itab
        }
      names(predList) <- nomi
      }
    }
  #
  backFun <- function(tab, diff, bc, quan, dat) {
    pr <- tab[,"mean"]
    se <- (tab[,"mean"]-tab[,"lower"])/quan
    if(diff>0) {
      x0 <- makeBoxCox(dat, bc[2], a=bc[1])
      if(n.ahead>0) {
        val <- diffinv(pr, differences=diff, xi=x0[length(x0)+1-c(diff:1)])
        pr <- val[(diff+1):length(val)]
        } else {
        for(j in 1:length(pr)) {
          if(!is.na(pr[j])) {
            jval <- diffinv(pr[j], differences=diff, xi=x0[j-c(diff:1)])
            pr[j] <- jval[length(jval)]
            }
          }
        }
      }
    tab[,"mean"] <- invBoxCox(pr, bc[2], a=bc[1])
    tab[,"lower"] <- invBoxCox(pr-quan*se, bc[2], a=bc[1])
    tab[,"upper"] <- invBoxCox(pr+quan*se, bc[2], a=bc[1])
    tab
    }
  #
  if(!is.null(object$call$unit)) {
    for(i in 1:length(predList)) {
      inam <- names(predList)[i]
      ibc <- object$call$box.cox[inam,]
      idiff <- object$call$ndiff[inam]
      if(idiff>0|ibc[2]!=1) {
        iquan <- qt((1+level)/2,object$equations[[inam]]$df.residual)
        idatList <- split(object$data.orig[subset,inam], object$data.orig[subset,object$call$unit], drop=T)
        iprList <- split(predList[[i]], predList[[i]][,object$call$unit], drop=T)
        igr <- names(iprList)
        iprOK <- list()
        for(j in 1:length(igr)) {
          iprOK[[j]] <- backFun(tab=iprList[[igr[j]]], diff=idiff, bc=ibc, quan=iquan, dat=idatList[[igr[j]]])
          }
        predList[[i]] <- do.call(rbind, iprOK)
        }
      }
    } else {
    for(i in 1:length(predList)) {
      inam <- names(predList)[i]
      ibc <- object$call$box.cox[inam,]
      idiff <- object$call$ndiff[inam]
      if(idiff>0|ibc[2]!=1) {
        iquan <- qt((1+level)/2,object$equations[[inam]]$df.residual)
        idat <- object$data.orig[subset,inam]
        ipr <- predList[[i]]
        predList[[i]] <- backFun(tab=ipr, diff=idiff, bc=ibc, quan=iquan, dat=idat)
        }
      }
    }
  #
  attr(predList,"level") <- level
  obj <- list(predicted=predList, observed=obsdat, n.ahead=n.ahead, call=object$call)
  class(obj) <- "predict.feVAR"
  obj
  }

# print method for class 'predict.feVAR'
print.predict.feVAR <- function(x, ...) {
  print(x$predicted)
  }

# compute ACF (auxiliary)
acfCalc <- function(x, unit.id, max.lag, bptest) {
  nomi <- colnames(x)
  if(is.null(unit.id)) {
    lagM <- max(1,min(apply(x,2,function(z){1.5+round(sqrt(sum(!is.na(z))))})))
    if(is.null(max.lag)|!is.numeric(max.lag)) {
      max.lag <- lagM
      } else {
      max.lag <- min(max(1,max.lag[1]),lagM)
      }
    mat <- nt <- matrix(nrow=max.lag+1,ncol=length(nomi))
    nt[1,] <- apply(x,2,function(z){sum(!is.na(z))})
    mat[1,] <- 1
    for(i in 1:length(nomi)) {
      ires <- x[,nomi[i]]
      for(j in 1:max.lag) {
        ires_lag <- c(rep(NA,j),ires)[1:length(ires)]
        ijD <- na.omit(cbind(ires,ires_lag))
        mat[j+1,i] <- cor(ijD[,1],ijD[,2])
        nt[j+1,i] <- nrow(ijD)
        }
      } 
    #mat <- nt <- matrix(nrow=max.lag+1,ncol=length(nomi))
    #for(i in 1:length(nomi)) {
    #  ires <- res[,nomi[i]]
    #  imat <- cbind(ires, LAG(ires,nlags=max.lag))
    #  mat[,i] <- cor(imat,use="complete.obs")[1,]
    #  nt[,i] <- rep(nrow(na.omit(imat)),nrow(nt))
    #  }
    } else {
    resList <- split(data.frame(x), unit.id)
    lagM <- max(1,min(sapply(resList, function(z){
      min(apply(z,2,function(k){round(1.5+sqrt(sum(!is.na(k))))}))
    })))
    if(is.null(max.lag)|!is.numeric(max.lag)) {
      max.lag <- lagM
      } else {
      max.lag <- min(max(1,max.lag[1]),lagM)
      }
    mat <- nt <- matrix(nrow=max.lag+1,ncol=length(nomi))
    mat[1,] <- 1
    nt[1,] <- apply(x,2,function(z){sum(!is.na(z))})
    for(i in 1:length(nomi)) {
      iresList <- lapply(resList,function(z){z[,nomi[i]]})
      ires <- do.call(c,iresList)
      for(j in 1:max.lag) {
        ijres_lag <- do.call(c,lapply(iresList,function(z){
          c(rep(NA,j),z)[1:length(z)]
        }))
        ijD <- na.omit(cbind(ires,ijres_lag))
        mat[j+1,i] <- cor(ijD[,1],ijD[,2])
        nt[j+1,i] <- nrow(ijD)
        }
      }
    #mat <- nt <- matrix(nrow=max.lag+1,ncol=length(nomi))
    #for(i in 1:length(nomi)) {
    #  iresList <- lapply(resList,function(z){z[,nomi[i]]})
    #  imat <- do.call(rbind, lapply(iresList, function(x){cbind(x,LAG(x,nlags=max.lag))}))
    #  mat[,i] <- cor(imat,use="complete.obs")[1,]
    #  nt[,i] <- rep(nrow(na.omit(imat)),nrow(nt))
    #  }
    }
  colnames(mat) <- colnames(nt) <- nomi
  rownames(mat) <- rownames(nt) <- 0:max.lag
  matOK <- na.omit(mat)
  attributes(matOK)$na.action <- NULL
  lagOK <- nrow(matOK)-1
  if(bptest) {
    bpt <- c()
    for(i in 1:ncol(matOK)) {
      ir <- matOK[-1,i]
      iT <- nt[-1,i]
      bpt <- cbind(bpt, iT*(iT+2)*cumsum(ir^2/(iT-(1:length(ir)))))
      }
    rownames(bpt) <- 1:lagOK
    colnames(bpt) <- nomi
    bpt_p <- round(apply(bpt, 2, function(z){1-pchisq(z,1:nrow(bpt))}),4)
    list(acf=matOK, acf.se=1/sqrt(nt), statistic=bpt, p.value=bpt_p)
    } else {
    list(acf=matOK, acf.se=1/sqrt(nt))
    }
  }

# autocorrelation test on residuals
autocorTest <- function(model, max.lag=NULL) {
  res <- residuals(model)
  if(!is.null(model$call$unit)) {
    acf1 <- acfCalc(res[,model$call$var.names,drop=F], unit.id=res[,model$call$unit], max.lag=max.lag, bptest=T)
    acf2 <- acfCalc(res[,model$call$var.names,drop=F]^2, unit.id=res[,model$call$unit], max.lag=max.lag, bptest=T)
    } else {
    acf1 <- acfCalc(res[,model$call$var.names,drop=F], unit.id=NULL, max.lag=max.lag, bptest=T)
    acf2 <- acfCalc(res[,model$call$var.names,drop=F]^2, unit.id=NULL, max.lag=max.lag, bptest=T)
    }
  obj <- list(residual=acf1[3:4], sq.residual=acf2[3:4])
  class(obj) <- "autocorTest.feVAR"
  obj
  }

# print method for class 'autocorTest.feVAR'
print.autocorTest.feVAR <- function(x, ...) {
  p1 <- x$residual$p.value
  p2 <- x$sq.residual$p.value
  mat <- cbind(residual=p1[nrow(p1),],sq.residual=p2[nrow(p2),])
  cat("p-values of the Ljung-Box test on residuals (lag ",nrow(p1),")",sep="","\n")
  print(mat)
  }

# companion matrix (auxiliary)
compMat <- function(Beta) {
  Cmat_up <- do.call(cbind,Beta)
  nc <- ncol(Cmat_up)-nrow(Cmat_up)
  if(nc>0) {
    Cmat_low <- matrix(0,nrow=nc,ncol=ncol(Cmat_up))
    for(i in 1:nc) Cmat_low[i,i] <- 1
    } else {
    Cmat_low <- NULL
    }
  Cmat <- rbind(Cmat_up,Cmat_low)
  rownames(Cmat) <- colnames(Cmat)
  Cmat
  }

# ACF of residuals (auxiliary)
residualACF <- function(x, nomi, max.lag=NULL, squared=F, signif=0.05, ylim, xlab, ylab, titles, mfrow, mar, mgp, ...) {
  res <- residuals(x)
  if(squared) res[,nomi] <- res[,nomi]^2
  if(is.null(titles)) titles <- nomi
  if(!is.null(x$call$unit)) {
    obj <- acfCalc(res[,nomi], unit.id=res[,x$call$unit], max.lag=max.lag, bptest=F)
    } else {
    obj <- acfCalc(res[,nomi], unit.id=NULL, max.lag=max.lag, bptest=F)
    }
  mat <- obj$acf
  serr <- obj$acf.se
  lagOK <- nrow(mat)-1
  zval <- qnorm(1-signif/2)
  opar <- par(no.readonly=T)
  if(is.null(mfrow)) mfrow <- n2mfrow(length(nomi))
  on.exit(par(opar))
  par(mfrow=mfrow, mar=mar, mgp=mgp)
  for(i in 1:length(nomi)) {
    ise <- serr[,i]
    if(is.null(ylim)) iylim <- range(c(zval*ise,-zval*ise,mat[,nomi[i]]))
    plot(0:lagOK, mat[,nomi[i]], type="h", xlab=xlab, ylab=ylab, ylim=iylim, main=titles[i], ...)
    abline(h=0)
    lines(0:lagOK, -zval*ise, lty=2, col=2)
    lines(0:lagOK, zval*ise, lty=2, col=2)
    }  
  }

# residual qqnorm (auxiliary)
residualQQ <- function(x, nomi, xlab, ylab, cex, titles, mfrow, mar, mgp, ...) {
  res <- residuals(x)
  if(is.null(titles)) titles <- nomi
  opar <- par(no.readonly=T)
  if(is.null(mfrow)) mfrow <- n2mfrow(length(nomi))
  on.exit(par(opar))
  par(mfrow=mfrow, mar=mar, mgp=mgp)
  for(i in 1:length(nomi)) {
    qqnorm(res[,nomi[i]], xlab=xlab, ylab=ylab, cex=cex, main=titles[i], ...)
    qqline(res[,nomi[i]])
    }
  }

# fitted vs residuals (auxiliary)
fittedVSresiduals <- function(x, nomi, squared, xlab, ylab, cex, titles, add.grid, mfrow, mar, mgp, ...) {
  res <- residuals(x)
  if(squared) res[,nomi] <- res[,nomi]^2
  fit <- fitted.values(x)
  if(is.null(titles)) titles <- nomi
  opar <- par(no.readonly=T)
  if(is.null(mfrow)) mfrow <- n2mfrow(length(nomi))
  on.exit(par(opar))
  par(mfrow=mfrow, mar=mar, mgp=mgp)
  for(i in 1:length(nomi)) {
    plot(fit[,nomi[i]], res[,nomi[i]], type="n", xlab=xlab, ylab=ylab, main=titles[i], ...)
    if(add.grid) grid()
    points(fit[,nomi[i]], res[,nomi[i]], cex=cex)
    abline(h=0, lty=2)
    box()
    }
  }

# spaghetti plot of residuals (auxiliary)
spagResid <- function(x, nomi, xlab, ylab, ylim, titles, add.grid, mfrow, mar, mgp, ...) {
  res <- residuals(x)
  nt <- length(unique(x$data.used[,x$call$time]))
  if(is.null(titles)) titles <- nomi
  if(is.null(x$call$unit)) {
    opar <- par(no.readonly=T)
    if(is.null(mfrow)) mfrow <- n2mfrow(length(nomi))
    on.exit(par(opar))
    par(mfrow=mfrow, mar=mar, mgp=mgp)
    for(i in 1:length(nomi)) {
      if(is.null(ylim)) iylim <- range(res[,nomi[i]],na.rm=T)
      plot(1:nt, res[,nomi[i]], ylim=iylim, type="n", xlab=xlab, ylab=ylab, main=titles[i], ...)
      if(add.grid) grid()
      lines2(1:nt, res[,nomi[i]], col="grey70")
      abline(h=0, lty=2)
      box()
      }
    } else {
    resList <- split(res,x$data.used[,x$call$unit])
    opar <- par(no.readonly=T)
    if(is.null(mfrow)) mfrow <- n2mfrow(length(nomi))
    on.exit(par(opar))
    par(mfrow=mfrow, mar=mar, mgp=mgp)
    for(i in 1:length(nomi)) {
      if(is.null(ylim)) iylim <- range(res[,nomi[i]],na.rm=T)
      plot(1:nt, resList[[1]][,nomi[i]], ylim=iylim, type="n", xlab=xlab, ylab=ylab, main=titles[i], ...)
      if(add.grid) grid()
      for(j in 1:length(resList)) {
        lines2(1:nt, resList[[j]][,nomi[i]], col="grey70")
        }
      abline(h=0, lty=2)
      box()
      }
    }
  }

# plot method for class 'feVAR'
plot.feVAR <- function(x, type=c("ts","acf","acf2","fitVSres","qqnorm"), var.names=NULL, max.nlags=NULL, signif=0.05, ylim=NULL, cex.points=0.6,
  add.grid=TRUE, xlab=NULL, ylab=NULL, titles=NULL, las=0, mfrow=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  #
  nomi <- x$call$var.names
  if(!is.null(var.names)) {
    if(is.numeric(var.names)) {
      nomi <- nomi[intersect(var.names,1:length(nomi))]
      } else {
      nomi <- intersect(var.names,nomi)
      }
    if(length(nomi)==0) nomi <- x$call$var.names
    }
  add.grid <- add.grid[1]
  if(is.na(add.grid)||(!is.logical(add.grid)|is.null(add.grid))) add.grid <- TRUE 
  if(is.numeric(type)) {
    type <- intersect(type,1:6)[1]    
    } else {
    type <- intersect(type,c("ts","acf","acf2","qqnorm","fitVSres"))[1]
    }
  if(is.null(type)||is.na(type)) {
    type <- 1
    warning("Argument 'type' contains no valid values: it is set to 'acf'",call.=F)
    }
  if(type=="ts"|type==1) {
    if(is.null(xlab)) xlab <- "time point"
    if(is.null(ylab)) ylab <- "residuals"
    spagResid(x, nomi=nomi, xlab=xlab, ylab=ylab, ylim=ylim, titles=titles, add.grid=add.grid, mfrow=mfrow, mar=mar, mgp=mgp, ...)
    } else if(type%in%c("acf",2,"acf2",3)) {
    if(!is.numeric(signif)) signif <- 0.05 else signif <- min(c(1,max(c(signif,0),na.rm=T)),na.rm=T)
    if(is.null(xlab)) xlab <- "lag"
    if(is.null(ylab)) ylab <- "ACF"
    residualACF(x, nomi=nomi, max.lag=max.nlags, squared=ifelse(type%in%c("acf2",3),T,F), signif=signif, ylim=ylim, xlab=xlab, ylab=ylab, titles=titles, mfrow=mfrow, mar=mar, mgp=mgp, ...)
    } else if(type=="qqnorm"|type==5) {
    if(is.null(xlab)) xlab <- "expected"
    if(is.null(ylab)) ylab <- "observed"
    residualQQ(x, nomi=nomi, xlab=xlab, ylab=ylab, cex=cex.points, titles=titles, mfrow=mfrow, mar=mar, mgp=mgp, ...)
    } else if(type=="fitVSres"|type==4) {
    if(is.null(xlab)) xlab <- "fitted values"
    if(is.null(ylab)) ylab <- "residuals"
    fittedVSresiduals(x, squared=F, nomi=nomi, xlab=xlab, ylab=ylab, cex=cex.points, titles=titles, add.grid=add.grid, mfrow=mfrow, mar=mar, mgp=mgp, ...)
    #} else if(type=="fitVSres2"|type==5) {
    #if(is.null(xlab)) xlab <- "fitted values"
    #if(is.null(ylab)) ylab <- expression(paste(residuals^2))
    #fittedVSresiduals(x, squared=T, nomi=nomi, xlab=xlab, ylab=ylab, cex=cex.points, titles=titles, add.grid=add.grid, mfrow=mfrow, mar=mar, mgp=mgp, ...)
    } else {
    stop("Argument 'type' must be one among 'ts' (1), 'acf' (2), 'acf2' (3), 'fitVSres' (4), and 'qqnorm' (5)",call.=F)
    }
  }

# plot method for class 'predict.feVAR'
plot.predict.feVAR <- function(x, var.names=NULL, unit.id, newdata=NULL, start=NULL, ylim=NULL, add.grid=TRUE, las=0, cex.axis=c(1,1), cex.lab=c(1,1), xlab=NULL, ylab=NULL, titles=NULL,
  obs.col="grey40", fit.col="dodgerblue", out.col="red", new.col="grey40", interval.col="grey70", obs.lty=1, fit.lty=5, out.lty=1, new.lty=2,
  obs.cex=0.4, new.cex=0.4, mfrow=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  #
  add.grid <- add.grid[1]
  if(is.na(add.grid)||(!is.logical(add.grid)|is.null(add.grid))) add.grid <- TRUE
  nomi <- x$call$var.names
  if(!is.null(var.names)) {
    if(is.numeric(var.names)) {
      nomi <- nomi[intersect(1:length(nomi),var.names)]
      } else {
      suppressWarnings(
        auxnam <- intersect(1:length(nomi),na.omit(as.numeric(var.names)))
        )
      if(length(auxnam)>0) var.names <- c(var.names,nomi[auxnam])
      nomi <- intersect(nomi,var.names)
      }
    if(length(nomi)==0) nomi <- x$call$var.names
    }
  if(!is.null(x$call$unit)) {
    glev <- levels(factor(x$observed[,x$call$unit]))
    if(missing(unit.id)||sum(!is.na(unit.id))==0) {
      unit.id <- glev[1]
      } else {
      wrn <- F
      if(is.numeric(unit.id)) {
        auxid <- intersect(1:length(glev), unit.id)
        if(length(auxid)!=length(unit.id)) wrn <- T
        if(length(auxid)>0) {
          unit.id <- glev[auxid[1]]
          } else {
          unit.id <- glev[1]  
          }
        } else {
        auxid <- intersect(glev, unit.id)
        if(length(auxid)!=length(unit.id)) wrn <- T
        if(length(auxid)>0) {
          unit.id <- auxid[1]
          } else {
          unit.id <- glev[1]  
          }
        }
      if(wrn) warning("Only unit '",unit.id,"' is displayed",call.=F)
      }
    }
  if(!is.null(newdata)) {
    if(!identical(class(newdata),"data.frame")) {
      stop("Argument 'newdata' must be an object of class 'data.frame'",call.=F)
      }
    auxchk <- setdiff(nomi,colnames(newdata))
    if(length(auxchk)>0) stop("Variable '",auxchk[1],"' not found in 'newdata",call.=F)
    #if(!is.null(x$call$time) && (x$call$time %in% colnames(newdata))==F) {
    if(!identical(x$call$time,"(time)") && (x$call$time %in% colnames(newdata))==F) {
      stop("Variable '",x$call$time,"' not found in 'newdata'",call.=F)
      }
    if(!is.null(x$call$unit)) {
      if(x$call$unit %in% colnames(newdata)) {
        datnew <- newdata[which(newdata[,x$call$unit] %in% unit.id),]
        } else {
        stop("Variable '",x$call$unit,"' not found in 'newdata'",call.=F)
        }
      } else {
      datnew <- newdata  
      }
    } else {
    datnew <- NULL
    }
  if(is.null(x$call$unit)) {
    dat <- x$observed
    if(x$n.ahead==0) {
      insam <- x$predicted
      outsam <- c()
      } else {
      outsam <- x$predicted
      insam <- outsam
      for(i in 1:length(insam)) {
        insam[[i]] <- insam[[i]][1:nrow(dat),]
        insam[[i]][] <- NA
        insam[[i]][,x$call$time] <- dat[,x$call$time]
        }
      }
    } else {
    dat <- x$observed[which(x$observed[,x$call$unit]%in%unit.id),]
    if(x$n.ahead==0) {
      insam <- lapply(x$predicted, function(z){z[which(z[,x$call$unit]%in%unit.id),]})
      outsam <- c()
      } else {
      outsam <- lapply(x$predicted, function(z){z[which(z[,x$call$unit]%in%unit.id),]})
      insam <- outsam
      for(i in 1:length(insam)) {
        insam[[i]] <- insam[[i]][1:nrow(dat),]
        insam[[i]][] <- NA
        insam[[i]][,x$call$time] <- dat[,x$call$time]
        }
      }
    }
  if(is.null(titles)) titles <- nomi
  if(is.null(xlab)) xlab <- ""
  if(is.null(ylab)) ylab <- "" 
  if(is.null(las)|!is.numeric(las)) las <- 0 else las <- max(0,min(3,round(las)))
  if(is.null(cex.axis)|!is.numeric(cex.axis)) cex.axis <- c(1,1)
  if(length(cex.axis)<2) cex.axis <- rep(cex.axis,2)
  if(is.null(cex.lab)|!is.numeric(cex.lab)) cex.lab <- c(1,1)
  if(length(cex.lab)<2) cex.lab <- rep(cex.lab,2)
  if(is.null(mfrow)) mfrow <- n2mfrow(length(nomi))
  tvar <- dat[,x$call$time]
  if(!is.null(start)) {
    ind0 <- which(tvar>=start)
    if(length(ind0)>0) {
      tvar <- tvar[ind0]
      dat <- dat[which(dat[,x$call$time]%in%tvar),]
      for(i in 1:length(insam)) {
        insam[[i]] <- insam[[i]][which(insam[[i]][,x$call$time]%in%tvar),]
        }
      }
    }
  if(x$n.ahead>0) tvar <- c(tvar, outsam[[1]][,x$call$time])
  n <- nrow(dat)
  N <- n+x$n.ahead
  opar <- par(no.readonly=T)
  on.exit(opar)
  par(mfrow=mfrow, mar=mar, mgp=mgp)
  for(i in 1:length(nomi)) {
    iobs <- dat[,nomi[i]]
    idat <- rbind(insam[[nomi[i]]][,c("mean","lower","upper")],
                  outsam[[nomi[i]]][,c("mean","lower","upper")])
    if(is.null(ylim)) {
      ilimy <- range(c(iobs,idat,datnew[,nomi[i]]),na.rm=T)
      } else {
      ilimy <- ylim
      }
    plot(tvar, idat[,1], type="n", xlab=xlab, ylab=ylab, ylim=ilimy, las=las, main=titles[i], cex.axis=cex.axis[1], cex.lab=cex.lab[1], yaxt="n", ...) 
    axis(2, cex.axis=cex.axis[2], cex.lab=cex.lab[2], las=las)
    if(add.grid) grid()
    if(x$n.ahead>0) {
      idat2 <- idat
      idat2[1:n,2:3] <- iobs[1:n]
      polygon(c(tvar[n:N],tvar[N:n]), c(idat2[n:N,2],rev(idat2[n:N,3])),
              border=NA, col=adjustcolor(interval.col,alpha.f=0.5))
      lines2(tvar[n:N], c(iobs[n],idat[(n+1):N,1]), col=out.col, lty=out.lty)
      lines2(tvar[1:n], iobs, col=obs.col, lty=obs.lty)
      points(tvar[1:n], iobs, cex=obs.cex, col=obs.col)
      } else {
      polygon(c(tvar,rev(tvar)), c(idat[,2],rev(idat[,3])),
              border=NA, col=adjustcolor(interval.col,alpha.f=0.5))
      lines2(tvar, idat[,1], col=fit.col, lty=fit.lty)
      lines2(tvar, iobs, col=obs.col, lty=obs.lty)
      points(tvar, iobs, cex=obs.cex, col=obs.col)
      }
    if(x$n.ahead>0&!is.null(datnew)) {
      if(identical(x$call$time,"(time)")) {
        ind <- intersect(1:x$n.ahead,1:nrow(datnew))
        } else {
        ind <- which(datnew[,x$call$time]%in%(tvar[(n+1):N]))
        }
      if(length(ind)>0) {
        lines2(tvar[n:(n+length(ind))], c(iobs[n],datnew[ind,nomi[i]]), col=new.col, lty=new.lty)
        points(tvar[(n+1):(n+length(ind))], datnew[ind,nomi[i]], cex=new.cex, col=new.col)
        }
      }
    box()
    }
  }

# plot method for class 'IRF.feVAR'
plot.IRF.feVAR <- function(x, from=NULL, to=NULL, n.ahead=NULL, labels=NULL, ylim=NULL, add.grid=TRUE, add.titles=FALSE, cex.points=0, mfrow=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  nomi <- dimnames(x[[1]])[[2]]
  if(is.null(from)) {
    from <- nomi
    } else {
    if(is.numeric(from)) {
      from <- nomi[intersect(1:length(nomi),from)] 
      } else {
      suppressWarnings(
        auxfrom <- intersect(1:length(nomi),na.omit(as.numeric(from)))
        )
      if(length(auxfrom)>0) from <- c(from,nomi[auxfrom])
      from <- intersect(nomi,from)
      }
    }
  if(length(from)==0) from <- nomi
  if(is.null(to)) {
    to <- nomi
    } else {
    if(is.numeric(to)) {
      to <- nomi[intersect(1:length(nomi),to)] 
      } else {
      suppressWarnings(
        auxto <- intersect(1:length(nomi),na.omit(as.numeric(to)))
        )
      if(length(auxto)>0) to <- c(to,nomi[auxto])
      to <- intersect(nomi,to)
      }
    }
  if(length(to)==0) to <- nomi
  lagM <- dim(x$irf)[1]-1
  if(is.null(n.ahead)|!is.numeric(n.ahead)) {
    n.ahead <- lagM
    } else {
    n.ahead <- min(round(max(c(1,n.ahead),na.rm=T)),lagM)
    }
  add.grid <- add.grid[1]
  if(is.na(add.grid)||(!is.logical(add.grid)|is.null(add.grid))) add.grid <- TRUE
  add.titles <- add.titles[1]
  if(is.na(add.titles)||(!is.logical(add.titles)|is.null(add.titles))) add.titles <- FALSE
  #
  onePlot <- function(x1, x2, ...) {
    irfOK <- x$irf[1:(n.ahead+1),x1,x2]
    irfOK_sx <- x$irf_sx[1:(n.ahead+1),x1,x2]
    irfOK_dx <- x$irf_dx[1:(n.ahead+1),x1,x2]
    if(is.null(ylim)) ylim <- range(c(irfOK_sx,irfOK_dx))
    lseq <- 0:n.ahead
    plot(lseq, irfOK, type="n", ylim=ylim, ...)
    if(add.grid) grid()
    abline(h=0, col=2, lty=2)
    if(!is.null(irfOK_sx)&!is.null(irfOK_dx)) {
      polygon(c(lseq,rev(lseq)), c(irfOK_sx,rev(irfOK_dx)),
              border=NA, col=adjustcolor("grey75",alpha.f=0.5))
      }
    lines(lseq, irfOK)
    points(lseq, irfOK, cex=cex.points)
    }
  #
  opar <- par(no.readonly=T)
  #if(is.null(mfrow)) mfrow <- c(length(from),length(to))
  if(is.null(mfrow)) mfrow <- n2mfrow(length(from)*length(to))
  on.exit(par(opar))
  par(mfrow=mfrow, mar=mar, mgp=mgp)
  for(i in 1:length(from)) {
    for(j in 1:length(to)) {
      if(from[i] %in% names(labels)) ijlab1 <- labels[from[i]] else ijlab1 <- from[i]
      if(to[j] %in% names(labels)) ijlab2 <- labels[to[j]] else ijlab2 <- to[j]
      if(add.titles) ijtit <- paste0("From: ",ijlab1," | To: ",ijlab2) else ijtit <- NULL
      onePlot(from[i], to[j], xlab=paste0(ijlab1," (lag)"), ylab=ijlab2, main=ijtit, ...)
      }
    }
  }

# simulate multivariate normal (auxiliary)
mvnsim <- function(n, mu, S) {
  p <- length(mu)
  z <- matrix(rnorm(n*p),nrow=n,ncol=p)
  W <- chol(S)
  z%*%W+matrix(rep(mu,n),ncol=p,byrow=T)
  }

# simulate beta matrices (auxiliary)
simulBeta <- function(model) {
  modList <- model$equations
  nlags <- length(model$Beta)
  lagList <- lapply(modList, getLags)
  xnam <- names(modList)
  betaList <- vector("list",length=nlags)
  names(betaList) <- 1:nlags
  for(i in 1:nlags) {
    imat <- matrix(0,nrow=length(xnam),ncol=length(xnam))
    rownames(imat) <- colnames(imat) <- xnam
    betaList[[i]] <- imat
    }
  for(i in 1:length(xnam)) {
    imod <- modList[[xnam[i]]]
    ilag <- lagList[[xnam[i]]]
    if(!is.null(ilag)) {
      istr <- rownames(ilag)
      isim <- mvnsim(1,imod$coefficients,vcov(imod))[1,]
      for(j in 1:length(istr)) {
        ijl <- ilag[j,"lag"]
        ijx <- ilag[j,"x"]
        betaList[[ijl]][ijx,xnam[i]] <- isim[istr[j]]
        }
      }
    }
  betaList
  }

# compute impulse response functions
IRF <- function(model, n.ahead=10, nboot=100, level=0.95, quiet=FALSE) {
  if(!identical(class(model),"feVAR")) stop("Argument 'model' must be an object of class 'feVAR'",call.=F)
  if(!is.numeric(n.ahead)) n.ahead <- 10 else n.ahead <- round(max(c(1,n.ahead),na.rm=T))
  if(is.na(n.ahead)|n.ahead=="Inf") n.ahead <- 10
  if(!is.numeric(nboot)) nboot <- 0 else nboot <- round(max(c(0,nboot),na.rm=T))
  if(!is.numeric(level)) level <- 0.05 else level <- max(0,min(level,1,na.rm=T),na.rm=T)
  quiet <- quiet[1]
  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  #
  Func <- function(betaList) {
    p <- length(betaList)
    Psi <- array(dim=c(n.ahead+1,m,m))
    dimnames(Psi) <- list(0:n.ahead,xnam,xnam)
    Psi[1,,] <- matrix(0,nrow=m,ncol=m)
    diag(Psi[1,,]) <- 1
    for(j in 1:n.ahead) {
      psi_j <- 0
      for(s in 1:min(j,p)) {
        psi_j <- psi_j+Psi[j+1-s,,]%*%betaList[[s]]
        }
      Psi[j+1,,] <- psi_j
      }
    Psi
    }
  #
  m <- length(model$equations)
  xnam <- names(model$equations)
  irf <- Func(model$Beta)
  if(nboot>0) {
    irf_sim <- array(dim=c(nboot,n.ahead+1,m,m))
    dimnames(irf_sim) <- list(NULL,0:n.ahead,xnam,xnam)
    for(i in 1:nboot) {
      if(quiet==F) {
        cat("\r","Bootstrap resample ",i,"/",nboot,sep="")
        flush.console()
        }
      irf_sim[i,,,] <- Func(simulBeta(model))
      }
    irf_sx <- apply(irf_sim,2:4,quantile,prob=(1-level)/2)
    irf_dx <- apply(irf_sim,2:4,quantile,prob=(1+level)/2)
    if(quiet==F) cat("\n")
    obj <- list(irf=irf, irf_sx=irf_sx, irf_dx=irf_dx)
    attr(obj,"level") <- level
    } else {
    obj <- list(irf=irf)
    }
  class(obj) <- "IRF.feVAR"
  obj
  }

# print method for class 'IRF.feVAR'
print.IRF.feVAR <- function(x, ...) {
  cat("Impulse response functions","\n")
  cat("  number of variables: ",dim(x$irf)[2],"\n",sep="")
  cat("  number of steps ahead: ",dim(x$irf)[1]-1,"\n",sep="")
  if(length(x)>1) {
    cat("  confidence intervals: yes (level ",attr(x,"level"),")","\n",sep="")
    } else {
    cat("  confidence intervals: no","\n")
    }
  }

# estimate h-step ahead prediction error
forecastError <- function(model, n.ahead=1, quiet=FALSE) {
  if(!identical(class(model),"feVAR")) stop("Argument 'model' must be an object of class 'feVAR'",call.=F)
  if(!is.numeric(n.ahead)) n.ahead <- 1 else n.ahead <- round(max(n.ahead,na.rm=T))
  if(is.na(n.ahead)|n.ahead<1|abs(n.ahead)=="Inf") n.ahead <- 1
  quiet <- quiet[1]
  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  obsdat <- model$data.orig
  trdat <- model$data.used
  if(is.null(model$call$unit)) {
    maxh <- max(nrow(na.omit(trdat))-5,1)
    } else {
    maxh <- min(sapply(split(trdat, trdat[,model$call$unit]), function(x){max(nrow(na.omit(x))-5,1)}))
    }
  if(n.ahead>maxh) {
    n.ahead <- maxh
    warning("Argument 'n.ahead' was set to the maximum possible value: ",maxh,call.=F)
    }
  #
  errFun <- function(pred, obs, h) {
    nomi <- colnames(pred)
    if(is.null(model$call$unit)) uval <- NULL else uval <- obs[,model$call$unit]
    rmse <- rmse0 <- mae <- mae0 <- rsq <- mape <- c()
    for(i in 1:length(nomi)) {
      iobs <- obs[,nomi[i]]
      ipred <- pred[,nomi[i]]
      rmse[i] <- sqrt(mean((iobs-ipred)^2, na.rm=T))
      mae[i] <- mean(abs(iobs-ipred), na.rm=T)
      ipred0 <- LAG(iobs, h, uval)
      rmse0[i] <- sqrt(mean((iobs-ipred0)^2, na.rm=T))
      mae0[i] <- mean(abs(iobs-ipred0), na.rm=T)
      #rsq[i] <- cor(iobs,ipred,use="pairwise.complete")^2
      if(sum(iobs<=0,na.rm=T)>0) mape[i] <- NA else mape[i] <- 100*mean(abs((iobs-ipred)/iobs), na.rm=T)
      }
    tab <- cbind(rmse=rmse, mae=mae, rmsse=rmse/rmse0, mase=mae/mae0, mape=mape)
    rownames(tab) <- nomi
    tab
    }
  #
  tval <- sort(unique(trdat[,model$call$time]))
  nT <- length(tval)
  pList <- vector("list",length=n.ahead)
  for(i in 1:n.ahead) {
    pList[[i]] <- matrix(nrow=nrow(trdat),ncol=length(model$call$var.names))
    colnames(pList[[i]]) <- model$call$var.names
    }
  for(i in 1:(nT-n.ahead)) {
    if(quiet==F) {
      cat("\r","Performing pseudo-prediction ",i,"/",nT-n.ahead,sep="")
      flush.console()
      }
    isub <- which(trdat[,model$call$time]%in%tval[1:i])
    ipr <- predict.feVAR(model, n.ahead=n.ahead, subset=isub)$predicted
    for(j in 1:n.ahead) {
      ijpr <- do.call(cbind, lapply(ipr, function(x){
        tt <- unique(x[,model$call$time])
        x[which(x[,model$call$time]==tt[j]),"mean"]
      }))
      pList[[j]][which(obsdat[,model$call$time]==tval[i+j]),] <- ijpr
      }
    }
  if(quiet==F) cat("\n")
  eList <- list()
  #pListOK <- list()
  for(i in 1:n.ahead) {
    eList[[i]] <- errFun(pred=pList[[i]], obs=model$data.orig, h=i)
    #pListOK[[i]] <- cbind(trdat[,c(model$call$unit,model$call$time)],pList[[i]])
    }
  names(eList) <- 1:n.ahead
  #names(pListOK) <- 1:n.ahead
  #res <- list(metrics=eList, predicted=pListOK)
  #class(res) <- "predError.feVAR"
  #res
  eList
  }