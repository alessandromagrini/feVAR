## DA FARE
#  - IRF + grafico
#  - error variance decomposition
#  - spaghetti plot


## FUNZIONI DA ESPORTARE
#  - preProcess
#  - unirootTest
#  - LAG
#  - feVAR
#  - residualPlot
#  - predictPlot
#
## METODI S3: print, summary, coef, fitted, residuals,
#             logLik, extractAIC, AIC, BIC, predict,
#             print.unirootTest


###  PREPROCESSING  ###

# log function (auxiliary)
logFun <- function(x) {
  ind <- which(x<=0)
  if(length(ind)>0) {
    x[ind] <- min(x[setdiff(1:length(x),ind)])/2
    }
  log(x)
  }

# apply box-cox transformation (auxiliary)
makeBoxCox <- function(x, par) {
  if(par==1) {
    x
    } else if(par==0) {
    logFun(x)  
    } else {
    (x^par-1)/par  
    }
  }

# linear interpolation (auxiliary)
linInterp <- function(x) {
  if(sum(!is.na(x))>=2) {
    approx(x, xout=1:length(x))$y
    } else {
    x
    }
  }

# unit root test for one variabile (auxiliary)
oneTest <- function(x, unit=NULL, max.lag=NULL) {
  if(is.null(unit)) {
    x <- na.omit(linInterp(x))
    n <- length(x)
    #if(n<5) stop("At least 5 observations are required",call.=F)
    } else {
    gr <- levels(factor(unit))
    nvet <- c()
    for(w in gr) {
      ind <- which(unit==w)
      x[ind] <- linInterp(x[ind])
      nvet[w] <- length(na.omit(x[ind]))
      #if(nvet[w]<5) stop("At least 5 observations are required for each unit",call.=F)
      }
    isOK <- which(!is.na(x))
    x <- x[isOK]
    unit <- unit[isOK]
    n <- min(nvet)
    }
  if(is.null(max.lag)) {
    #max.lag <- min(n-3,trunc((n-1)^(1/3)))
    max.lag <- round(sqrt(n))
    } else {
    if(length(max.lag)>1) max.lag <- max.lag[1]
    if(!is.numeric(max.lag) || max.lag!=round(max.lag) || max.lag<0) stop("Argument 'max.lag' must be a non-negative integer value")
    if(max.lag>n-3) {
      max.lag <- n-3
      #warning("Argument 'max.lag' was set to ",n-3)
      }
    }
  if(is.null(unit)) {
    res <- res1 <- adfFun(x=x, max.lag=max.lag)
    res2 <- kpssFun(x=x, max.lag=max.lag)
    for(i in 1:3) {
      res[[i]] <- c(adf=res1[[i]],kpss=res2[[i]])
      }
    } else {
    gr <- levels(factor(unit))
    res1 <- res2 <- vector("list",length=3)
    for(w in gr) {
      ind <- which(unit==w)
      iadf <- adfFun(x=x[ind], max.lag=max.lag)
      ikpss <- kpssFun(x=x[ind], max.lag=max.lag)
      for(j in 1:length(res1)) {
        res1[[j]] <- c(res1[[j]],iadf[[j]])
        res2[[j]] <- c(res2[[j]],ikpss[[j]])
        }
      }
    #
    pvalComb <- function(x) {
      x[which(x<=0)] <- 1e-8
      x[which(x>=1)] <- 1-1e-8
      ind <- which(!is.na(x))
      if(length(ind)>0) {
        m <- length(ind)
        logp <- qnorm(x[ind])
        rhat <- 1-var(logp)
        rstar <- max(rhat,-1/(m-1))
        auxz <- sum(logp)/sqrt(m*(1+(m-1)*(rstar+0.2*sqrt(2/(m+1))*(1-rstar))))
        #auxz <- sum(logp)/sqrt(m)
        c(x,'(combined)'=2*pnorm(-abs(auxz)))
        } else {
        c(x,'(combined)'=NaN)
        }
      }
    #
    res1 <- lapply(res1, function(z){names(z)<-gr; z})
    names(res1) <- names(iadf)
    res1$p.value <- pvalComb(res1$p.value)
    res2 <- lapply(res2, function(z){names(z)<-gr; z})
    names(res2) <- names(ikpss)
    res2$p.value <- pvalComb(res2$p.value)
    res <- res1
    for(i in 1:3) {
      res[[i]] <- cbind(adf=res1[[i]],kpss=res2[[i]])
      }
    }
  res
  }

# function for adf test (auxiliary)  
adfFun <- function(x, max.lag) {
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
      res <- lm(yt~xt1+tt+yt1)
      } else {
      res <- lm(yt~xt1+tt)
      }
    suppressWarnings(
      res.sum <- summary.lm(res)$coefficients
      )
    if(nrow(res.sum)>=2) {
      STAT <- res.sum[2,1]/res.sum[2,2]
      table <- -1*cbind(c(4.38, 4.15, 4.04, 3.99, 3.98, 3.96),
                        c(3.95, 3.8, 3.73, 3.69, 3.68, 3.66),
                        c(3.6, 3.5, 3.45, 3.43, 3.42, 3.41),
                        c(3.24, 3.18, 3.15, 3.13, 3.13, 3.12),
                        c(1.14, 1.19, 1.22, 1.23, 1.24, 1.25),
                        c(0.8, 0.87, 0.9, 0.92, 0.93, 0.94),
                        c(0.5, 0.58, 0.62, 0.64, 0.65, 0.66),
                        c(0.15, 0.24, 0.28, 0.31, 0.32, 0.33))
      tablen <- dim(table)[2]
      tableT <- c(25, 50, 100, 250, 500, 1e+05)
      tablep <- c(0.01, 0.025, 0.05, 0.1, 0.9, 0.95, 0.975, 0.99)
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
kpssFun <- function(x, max.lag) {
  #
  doKPSS <- function(lag) {
    n <- length(x)
    #if(trend==T) {
    t <- 1:n
    e <- residuals.lm(lm(x ~ t))
    table <- c(0.216, 0.176, 0.146, 0.119)
    #  } else {
    #  e <- residuals.lm(lm(x ~ 1))
    #  table <- c(0.739, 0.574, 0.463, 0.347)
    #  }
    tablep <- c(0.01, 0.025, 0.05, 0.1)
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
    STAT <- eta/(s2+2*k/n)
    PVAL <- approx(table,tablep,STAT,rule=2)$y
    c(statistic=STAT, p.value=PVAL)
    }
  #
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
unirootTest <- function(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0, max.lag=NULL) {
  dataD <- preProcess(var.names=var.names, unit=unit, time=time, data=data, box.cox=box.cox, ndiff=ndiff)
    #imputation=FALSE)
  if(is.null(unit)) gr <- NULL else gr <- dataD[,unit]
  max.lag <- max.lag[1]
  if(!is.numeric(max.lag)) max.lag <- NULL else max.lag <- max(0,ceiling(max.lag))
  testList <- list()
  for(i in 1:length(var.names)) {
    iadf <- oneTest(x=dataD[,var.names[i]], unit=gr, max.lag=max.lag)
    iadf$box.cox <- unname(attr(dataD,"box.cox")[var.names[i]])
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
    if(is.matrix(pval)) pval["(combined)",] else pval
    })
  print(t(tab))
  }

# function for differencing (auxiliary)
diffFun <- function(var.names, time, data, ndiff) {
  newdat <- data
  if(!is.null(time)) newdat <- newdat[order(newdat[,time]),]
  n <- nrow(data)
  for(i in 1:length(var.names)) {
    if(ndiff[var.names[i]]>0) {
      idat <- linInterp(data[,var.names[i]])
      idat_lag <- c(rep(NA,ndiff[i]),idat)[1:length(idat)]
      newdat[,var.names[i]] <- idat-idat_lag
      }
    }
  attr(newdat,"ndiff") <- ndiff
  newdat
  }

# pre-processing
preProcess <- function(var.names, unit=NULL, time=NULL, data, box.cox=1, ndiff=0) {
  #imputation=TRUE, em.control=list(nlags=NULL,tol=1e-4,maxit=1000,quiet=FALSE)) {
  #
  if(missing(data)) stop("Argument 'data' is missing")
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  if(missing(var.names)) stop("Argument 'var.names' is missing")
  if(!is.character(var.names) || length(var.names)<1) stop("Argument 'var.names' must be a character vector of length 1 or greater")
  auxchk <- setdiff(var.names,colnames(data))  
  if(length(auxchk)>0) stop("Unknown variable '",auxchk[1],"' in argument 'var.names'")
  #
  unit <- unit[1]
  if(!is.null(unit)&&is.na(unit)) unit <- NULL
  if(!is.null(unit)) {
    if(length(setdiff(unit,colnames(data)))>0) stop("Unknown variable '",unit,"' provided to argument 'unit'")
    if(length(intersect(unit,var.names))>0) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'")
    if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values")
    }
  #
  time <- time[1]
  if(!is.null(time)&&is.na(time)) time <- NULL
  if(!is.null(time)) {
    if(length(setdiff(time,colnames(data)))>0) stop("Unknown variable '",time,"' provided to argument 'time'")
    if(length(intersect(time,var.names))>0) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'")
    if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    if(sum(is.na(data[,time]))>0) stop("Variable '",time,"' provided to argument 'time' contains missing values")
    }
  #
  if(!is.numeric(box.cox)) stop("Argument 'box.cox' must be a numeric value or vector")
  if(length(box.cox)==1) {
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
  if(sum(box.cox<0)>0) stop("Argument 'box.cox' must contain non-negative real values")
  #
  if(!is.numeric(ndiff)) stop("Argument 'ndiff' must be a numeric value or vector")
  if(length(ndiff)==1) {
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
  if(sum(ndiff<0|round(ndiff)!=ndiff)>0) stop("Argument 'ndiff' must contain non-negative integer values")
  #
  #imputation <- imputation[1]
  #if(is.na(imputation)||(!is.logical(imputation)|is.null(imputation))) imputation <- TRUE
  #
  dataL <- data
  for(i in 1:length(var.names)) {
    ilam <- box.cox[var.names[i]]
    if(ilam!=1 & sum(data[,var.names[i]]<0,na.rm=T)>0) {
      box.cox[var.names[i]] <- ilam <- 1
      warning("Box-Cox transformation not applied to variable '",var.names[i],"'",call.=F)
      }
    dataL[,var.names[i]] <- makeBoxCox(data[,var.names[i]],ilam)
    }
  #
  if(is.null(unit)) {
    n <- nrow(data)
    for(i in 1:length(var.names)) {
      if(ndiff[var.names[i]]>n-5) {
        ndiff[var.names[i]] <- 0
        warning("Differencing not applied to variable '",var.names[i],"'",call.=F)
        }
      }
    dataD <- diffFun(var.names=var.names, time=time, data=dataL, ndiff=ndiff)
    attr(dataD,"box.cox") <- box.cox
    attr(dataD,"ndiff") <- ndiff
    if(max(ndiff)>0) {
      res <- dataD[setdiff(1:nrow(data),1:max(ndiff)),,drop=F]
      } else {
      res <- dataD
      }
    } else {
    data[,unit] <- factor(data[,unit])
    dataD <- dataL
    isNA <- c()
    gr <- levels(data[,unit])
    val0 <- matrix(nrow=length(gr),ncol=length(var.names))
    rownames(val0) <- gr
    colnames(val0) <- var.names
    n_gr <- c()
    for(w in 1:length(gr)) {
      ind <- which(data[,unit]==gr[w])
      n_gr[w] <- length(ind)
      if(max(ndiff)>0) isNA <- c(isNA, ind[1]:ind[max(ndiff)])
      dataD[ind,] <- diffFun(var.names=var.names, time=time, data=dataL[ind,], ndiff=ndiff)
      val0[w,] <- as.numeric(data[ind[1],var.names])
      }
    n <- max(n_gr)
    for(i in 1:length(var.names)) {
      if(ndiff[var.names[i]]>n-5) {
        ndiff[var.names[i]] <- 0
        warning("Differencing not applied to variable '",var.names[i],"'",call.=F)
        }
      }
    attr(dataD,"box.cox") <- box.cox
    attr(dataD,"ndiff") <- ndiff
    res <- dataD[setdiff(1:nrow(data),isNA),,drop=F]
    }
  res
  #
  #if(imputation & sum(is.na(res[,var.names]))>0) {
  #  nlags <- em.control$nlags[1]
  #  if(!is.numeric(nlags)) {
  #    nlags <- NULL
  #    } else {
  #    nlags <- round(abs(nlags))
  #    }
  #  tol <- em.control$tol[1]
  #  if(!is.numeric(tol)|is.null(tol)) tol <- 1e-4
  #  if(tol<=0) tol <- 1e-4
  #  maxit <- em.control$maxit[1]
  #  if(!is.numeric(maxit)|is.null(maxit)) maxit <- 1000
  #  if(maxit<=0) maxit <- 1000 else maxit <- ceiling(maxit)
  #  quiet <- em.control$quiet[1]
  #  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  #  resI <- EMimput(var.names=var.names, unit=unit, time=time, data=res,
  #    nlags=nlags,tol=tol,maxit=maxit,quiet=quiet)
  #  resI$data.imputed
  #  } else {
  #  res
  #  }
  }


###  VAR MODELING  ###

# generate matrix of lags
LAG <- function(x, p, unit=NULL, cut=0, ...) {
  if(p>0) {
    #
    lfun <- function(v) {
      n <- length(v)
      res <- matrix(nrow=n,ncol=p)
      for(i in 1:p) {
        res[,i] <- c(rep(NA,i),v)[1:length(v)]
        }
      colnames(res) <- 1:p
      if(cut>0) {
        res[1:cut,] <- NA
        }
      res
      }
    #
    if(is.null(unit)) {
      lfun(x)
      } else {
      res <- matrix(nrow=length(x),ncol=p)
      gr <- levels(factor(unit))
      for(i in 1:length(gr)) {
        ind <- which(unit==gr[i])
        res[ind,] <- lfun(x[ind])
        }
      colnames(res) <- 1:p
      res
      }
    } else {
    x  
    }
  }

# fit a single regression (auxiliary)
regFit <- function(y.name, var.names, nlags, exogenous, data, unit, add.intercept, cut) {
  nomi <- c(y.name,var.names)
  names(nlags) <- nomi
  if(!is.null(exogenous)) {
    estr <- paste0(paste(exogenous,collapse="+"),"+")
    } else {
    estr <- ""
    }
  if(is.null(unit)) {
    xstr <- c()
    for(i in 1:length(nomi)) {
      if(nlags[i]>=1) {
        if(cut>0) {
          xstr[i] <- paste("LAG(",nomi[i],",",nlags[i],",cut=",cut,")", sep="")
          } else {
          xstr[i] <- paste("LAG(",nomi[i],",",nlags[i],")", sep="")
          }
        }
      }
    xstr <- na.omit(xstr)
    if(length(xstr)>0) {
      xstrOK <- paste(xstr, collapse="+")
      } else {
      xstrOK <- "1"  
      }
    if(add.intercept) {
      form <- formula(paste(y.name,"~",estr,xstrOK, sep=""))
      } else {
      form <- formula(paste(y.name,"~-1+",estr,xstrOK, sep=""))  
      }
    } else {
    xstr <- c()
    for(i in 1:length(nomi)) {
      if(nlags[i]>=1) {
        if(cut>0) {
          xstr[i] <- paste("LAG(",nomi[i],",",nlags[i],",",unit,",cut=",cut,")", sep="")
          } else {
          xstr[i] <- paste("LAG(",nomi[i],",",nlags[i],",",unit,")", sep="")
          }
        }
      }
    xstr <- na.omit(xstr)
    if(length(xstr)>0) {
      xstrOK <- paste(xstr, collapse="+")
      } else {
      xstrOK <- "1"  
      }
    if(add.intercept) {
      form <- formula(paste(y.name,"~-1+",unit,"+",estr,xstrOK, sep=""))
      } else {
      form <- formula(paste(y.name,"~-1+",estr,xstrOK, sep=""))        
      }
    }
  mod <- lm(form, data=data)
  mod$call$formula <- form
  mod$nlags <- nlags
  mod
  }

# variable selection in a single regression (auxiliary) 
lagSelect <- function(y.name, var.names, exogenous, data, max.nlags, unit, add.intercept, penaltyFun) {
  ic.current <- Inf
  lags.current <- rep(max.nlags, length(var.names)+1)
  fine <- 0
  while(fine==0) {
    icval <- c()
    for(i in 1:length(lags.current)) {
      lags <- lags.current
      lags[i] <- lags[i]-1
      mod <- regFit(y.name=y.name, var.names=var.names, nlags=lags, exogenous=exogenous, data=data, unit=unit, add.intercept=add.intercept, cut=max.nlags)
      ll <- logLik(mod)
      icval[i] <- -2*ll[1]+penaltyFun(attr(ll,"df"), nobs(mod))
      }
    ind <- which.min(icval)
    ic.test <- icval[ind]
    if(ic.test<ic.current) {
      lags.current[ind] <- lags.current[ind]-1
      ic.current <- ic.test
      if(sum(lags.current>=1)==0) fine <- 1
      } else {
      fine <- 1  
      }
    }
  regFit(y.name=y.name, var.names=var.names, nlags=lags.current, exogenous=exogenous, data=data, unit=unit, add.intercept=add.intercept, cut=0)
  }

# fit a VAR model (auxiliary)
fitVAR <- function(var.names, unit, exogenous, data, max.nlags, nlags, add.intercept, penaltyFun, cut, quiet=quiet) {
  mod <- list()
  if(is.null(nlags)) {
    for(i in 1:length(var.names)) {
      iy <- var.names[i]
      ix <- setdiff(var.names, iy)
      if(quiet==F) {
        cat("\r","Selecting lags for equation ",i,"/",length(var.names),sep="")
        flush.console()
        }
      mod[[i]] <- lagSelect(y.name=iy, var.names=ix, exogenous=exogenous, data=data, max.nlags=max.nlags, unit=unit, add.intercept=add.intercept, penaltyFun=penaltyFun)
      }
    if(quiet==F) cat("\n")
    } else {
    for(i in 1:length(var.names)) {
      iy <- var.names[i]
      ix <- setdiff(var.names, iy)
      mod[[i]] <- regFit(y.name=iy, var.names=ix, nlags=rep(nlags,length(var.names)), exogenous=exogenous, data=data, unit=unit, add.intercept=add.intercept, cut=cut)
      }
    }
  names(mod) <- var.names
  mod
  }

# compute maximum lag length (auxiliary)
lagCalc <- function(var.names, unit, data, add.intercept) {
  #
  dfCalc <- function(n) {
    xstr <- c()
    for(i in 1:length(var.names)) {
      xstr[i] <- paste("LAG(",var.names[i],",",n,",",unit,")", sep="")
      }
    if(add.intercept) {
      if(is.null(unit)) {
        form <- formula(paste("~",paste(xstr, collapse="+"), sep=""))
        } else {
        form <- formula(paste("~",unit,"+",paste(xstr, collapse="+"), sep=""))
        }
      } else {
      form <- formula(paste("~-1+",paste(xstr, collapse="+"), sep=""))  
      }
    M <- model.matrix(form, data=data)
    nrow(M)-ncol(M)
    }
  #
  if(is.null(unit)) {
    nlags <- nrow(data)-3
    } else {
    nlags <- min(sapply(split(data, data[,unit]),nrow))-3
    }
  df <- dfCalc(nlags)
  if(df>=3|nlags<=0) fine <- 1 else fine <- 0
  while(fine==0) {
    nlags <- nlags-1
    df <- dfCalc(nlags)
    if(df>=3|nlags<=0) fine <- 1
    }
  nlags
  }

# MASTER FUNCTION
feVAR <- function(var.names, unit=NULL, time=NULL, exogenous=NULL, data, max.nlags=NULL, nlags=NULL, add.intercept=TRUE, local.adapt=FALSE, ic="bic", quiet=FALSE) {
  #
  if(missing(data)) stop("Argument 'data' is missing")
  if(!identical(class(data),"data.frame")) stop("Argument 'data' must be a data.frame")
  if(missing(var.names)) stop("Argument 'var.names' is missing")
  if(!is.character(var.names) || length(var.names)<1) stop("Argument 'var.names' must be a character vector of length 1 or greater")
  auxchk <- setdiff(var.names,colnames(data))  
  if(length(auxchk)>0) stop("Unknown variable '",auxchk[1],"' in argument 'var.names'")
  #
  unit <- unit[1]
  if(!is.null(unit)&&is.na(unit)) unit <- NULL
  if(!is.null(unit)) {
    if(length(setdiff(unit,colnames(data)))>0) stop("Unknown variable '",unit,"' provided to argument 'unit'")
    if(length(intersect(unit,var.names))>0) stop("Variable '",unit,"' appears in both arguments 'var.names' and 'unit'")
    #if(sum(is.na(data[,unit]))>0) stop("Variable '",unit,"' provided to argument 'unit' contains missing values")
    }
  #
  time <- time[1]
  if(!is.null(time)&&is.na(time)) time <- NULL
  if(!is.null(time)) {
    if(length(setdiff(time,colnames(data)))>0) stop("Unknown variable '",time,"' provided to argument 'time'")
    if(length(intersect(time,var.names))>0) stop("Variable '",time,"' appears in both arguments 'var.names' and 'time'")
    if(length(intersect(time,unit))>0) stop("Variable '",time,"' appears in both arguments 'unit' and 'time'")
    if(!is.numeric(data[,time])&!identical(class(data[,time]),"Date")) stop("Variable '",time,"' must be numeric or of class 'Date'")
    }
  if(!is.null(exogenous)) {
    auxch1 <- setdiff(exogenous,colnames(data))
    if(length(auxch1)>0) stop("Unknown variable '",auxch1[1],"' provided to argument 'exogenous'")
    auxch2 <- intersect(exogenous,var.names)
    if(length(auxch2)>0) stop("Variable '",auxch2[1],"' appears in both arguments 'var.names' and 'exogenous'")
    auxch3 <- intersect(exogenous,unit)
    if(length(auxch3)>0) stop("Variable '",auxch3[1],"' appears in both arguments 'unit' and 'exogenous'")
    if(!is.null(time)) {
      auxch4 <- intersect(exogenous,time)
      if(length(auxch4)>0) stop("Variable '",auxch4[1],"' appears in both arguments 'time' and 'exogenous'")
      }
    }
  add.intercept <- add.intercept[1]
  if(is.na(add.intercept)||(!is.logical(add.intercept)|is.null(add.intercept))) add.intercept <- TRUE
  local.adapt <- local.adapt[1]
  if(is.na(local.adapt)||(!is.logical(local.adapt)|is.null(local.adapt))) local.adapt <- FALSE 
  quiet <- quiet[1]
  if(is.na(quiet)||(!is.logical(quiet)|is.null(quiet))) quiet <- FALSE
  #
  data <- data[complete.cases(data[,c(unit,time,var.names,exogenous)]),]
  #
  if(!is.null(unit)) data[,unit] <- factor(data[,unit])
  if(!is.null(time)) {
    if(is.null(unit)) {
      data <- data[order(data[,time]),]
      } else {
      cou <- levels(data[,unit])
      for(i in 1:length(cou)) {
        ind <- which(data[,unit]==cou[i])
        idat <- data[ind,]
        idat <- idat[order(idat[,time]),]
        data[ind,] <- idat
        }
      }
    }
  laglim <- lagCalc(var.names=var.names, unit=unit, data=data, add.intercept=add.intercept)
  if(laglim<=0) stop("Insufficient number of observations")
  if(is.null(nlags)) {
    if(is.null(max.nlags)) {
      max.nlags <- laglim
      } else {
      if(!is.numeric(max.nlags)) {
        warning("Argument 'max.nlags' has been set to the maximum possible: ",laglim,sep="",call.=F)
        max.nlags <- laglim
        } else {
        max.nlags <- round(max.nlags[1])
        }
      if(max.nlags>laglim|max.nlags<=0) warning("Argument 'max.nlags' has been set to the maximum possible: ",laglim,sep="",call.=F)
      max.nlags <- max(1,min(max.nlags,laglim))
      }
    } else {
    if(!is.numeric(nlags)) {
      warning("Argument 'nlags' has been set to the maximum possible: ",laglim,sep="",call.=F)
      nlags <- laglim
      } else {
      nlags <- round(nlags[1])
      }
    if(nlags>laglim|nlags<=0) warning("Argument 'nlags' has been set to the maximum possible: ",laglim,sep="",call.=F)
    nlags <- max(1,min(nlags,laglim))
    }
  p <- length(var.names)
  ic <- intersect(tolower(ic),c("bic","hqic","aic"))
  if(length(ic)==0) ic <- "bic" else ic <- ic[1]
  if(ic=="aic") {
    penaltyFun <- function(npar,n){2*npar}
    } else if(ic=="hqic") {
    penaltyFun <- function(npar,n){2*npar*log(log(n))}
    } else {
    penaltyFun <- function(npar,n){npar*log(n)}
    }
  if(is.null(nlags) & !is.null(max.nlags) & local.adapt==F) {
    icval <- c()
    for(i in 1:max.nlags) {
      if(quiet==F) cat("VAR with ",i," lags.",sep="")
      auxmod <- fitVAR(var.names=var.names, unit=unit, exogenous=exogenous, data=data, nlags=i, add.intercept=add.intercept, cut=max.nlags, quiet=T)
      ll <- sum(sapply(auxmod, function(x){logLik(x)[1]}))
      npar <- sum(sapply(auxmod, function(x){attr(logLik(x),"df")}))
      icval[i] <- -2*ll[1]+penaltyFun(npar, sum(sapply(auxmod,nobs)))
      if(quiet==F) cat(" ",toupper(ic),": ",icval[i],sep="","\n")
      }
    lagOK <- which.min(icval)
    if(quiet==F) cat("Number of lags selected: ",lagOK,sep="","\n")
    mod <- fitVAR(var.names=var.names, unit=unit, exogenous=exogenous, data=data, nlags=lagOK, add.intercept=add.intercept, cut=0, quiet=T)
    } else {
    mod <- fitVAR(var.names=var.names, unit=unit, exogenous=exogenous, data=data, max.nlags=max.nlags, nlags=nlags, add.intercept=add.intercept, penaltyFun=penaltyFun, cut=0, quiet=quiet)
    }
  modOK <- list(models=mod)
  modOK$call <- list(var.names=var.names, unit=unit, time=time, exogenous=exogenous, add.intercept=add.intercept)
  modOK$nlags <- lapply(mod, function(y){y$nlags})
  #
  if(!is.null(unit) & add.intercept) {
    modOK$intercepts <- do.call(cbind,lapply(modOK$models, function(x){
      b <- x$coefficients
      b[paste0(unit,levels(data[,unit]))]
      }))
    } else {
    modOK$intercepts <- sapply(modOK$models, function(x){
      unname(x$coefficients["(Intercept)"])
      })
    }
  #
  lagM <- max(unlist(modOK$nlags))
  if(lagM>0) {
    betaList <- vector("list",length=lagM)
    names(betaList) <- 1:lagM
    for(i in 1:lagM) {
      imat <- matrix(0,nrow=length(var.names),ncol=length(var.names))
      rownames(imat) <- colnames(imat) <- var.names
      betaList[[i]] <- imat
      }
    for(i in 1:length(var.names)) {
      imod <- modOK$models[[var.names[i]]]
      icoef <- lapply(paste0("LAG\\(",var.names,","), function(x){
        imod$coefficients[grep(x, names(imod$coefficients))]
        })
      for(j in 1:length(icoef)) {
        if(length(icoef[[j]])>0) {
          for(k in 1:length(icoef[[j]])) {
            betaList[[k]][j,i] <- icoef[[j]][k]
            }
          }
        }
      }
    modOK$Beta <- betaList
    }
  #
  resMat <- do.call(cbind, lapply(mod, function(x){
    x$residuals[rownames(data)]
    }))
  modOK$Sigma <- cov(resMat,use="pairwise.complete.obs")
  #
  ll <- sum(sapply(mod, function(x){logLik(x)[1]}))
  npar <- sum(sapply(mod, function(x){attr(logLik(x),"df")}))
  n <- sum(sapply(mod,nobs))
  ic <- c(aic=-2*ll+2*npar, bic=-2*ll+npar*log(n), hqic=-2*ll+2*npar*log(log(n)))
  modOK$ic <- ic
  #
  modOK$data <- data
  class(modOK) <- "feVAR"
  modOK
  }

# print method for class 'feVAR'
print.feVAR <- function(x, ...) {
  cat("Fixed effects VAR on ",length(x$call$var.names)," variables",sep="","\n")
  if(is.null(x$call$unit)) {
    ng <- 1
    nt <- nrow(x$data)
    } else {
    ng <- nlevels(x$data[,x$call$unit])
    nt <- max(sapply(split(x$data,x$data[,x$call$unit]),nrow))
    }
  cat("  number of units: ",ng,sep="","\n")
  cat("  number of time points: ",nt,sep="","\n")
  }

# coef method for class 'feVAR'
coef.feVAR <- function(object, ...) {
  lapply(object$models, function(x){x$coefficients})
  }

# summary method for class 'feVAR'
summary.feVAR <- function(object, ...) {
  lapply(object$models, summary, ...)
  }

# fitted method for class 'feVAR'
fitted.feVAR <- function(object, ...) {
  fit <- lapply(object$models, function(x){x$fitted.values})
  nomiObs <- rownames(object$data)
  for(i in 1:length(fit)) {
    ifit <- fit[[i]][nomiObs]
    names(ifit) <- nomiObs
    fit[[i]] <- ifit
    }
  tab <- do.call(cbind, fit)
  nomi <- colnames(tab)
  tabOK <- data.frame(object$data[,object$call$unit],
                      object$data[,object$call$time],tab)
  colnames(tabOK) <- c(object$call$unit,object$call$time,nomi)
  tabOK
  }

# residuals method for class 'feVAR'
residuals.feVAR <- function(object, ...) {
  res <- lapply(object$models, function(x){x$residuals})
  nomiObs <- rownames(object$data)
  for(i in 1:length(res)) {
    ires <- res[[i]][nomiObs]
    names(ires) <- nomiObs
    res[[i]] <- ires
    }
  tab <- do.call(cbind, res)
  nomi <- colnames(tab)
  tabOK <- data.frame(object$data[,object$call$unit],
                      object$data[,object$call$time],tab)
  colnames(tabOK) <- c(object$call$unit,object$call$time,nomi)
  tabOK
  }

# predict method for class 'feVAR'
predict.feVAR <- function(object, newdata=NULL, n.ahead=0, unit=NULL, level=0.95, ...) {
  if(!is.numeric(level)) level <- 0.95 else level <- min(c(1,max(c(level,0),na.rm=T)),na.rm=T)
  if(!is.numeric(n.ahead)) n.ahead <- 0 else n.ahead <- round(max(n.ahead,na.rm=T))
  if(is.na(n.ahead)|n.ahead<0|abs(n.ahead)=="Inf") n.ahead <- 0
  #
  ## CHECK ARGUMENTS: unit
  #
  if(n.ahead==0) {
    pred <- lapply(object$models, predict, newdata=newdata, interval="prediction", level=level)
    dat0 <- object$data
    nomiObs <- rownames(dat0)
    for(i in 1:length(pred)) {
      ipred <- data.frame(pred[[i]])
      ipredOK <- cbind(dat0[,object$call$unit],dat0[,object$call$time],ipred[nomiObs,])
      colnames(ipredOK) <- c(object$call$unit,object$call$time,"mean","lower","upper")
      rownames(ipredOK) <- nomiObs
      pred[[i]] <- ipredOK
      }
    if(!is.null(object$call$unit)&!is.null(unit)) {
      predList <- lapply(pred,function(z){z[which(z[,object$call$unit]%in%unit),]})
      } else {
      predList <- pred
      }
    } else {
    nomi <- object$call$var.names
    if(is.null(object$call$unit)|!is.null(unit)) {
      if(!is.null(object$call$unit)&!is.null(unit)) {
        dat0 <- object$data[which(object$data[,object$call$unit]%in%unit),]
        } else {
        dat0 <- object$data
        }
      obs0 <- dat0[nrow(dat0),]
      datOK <- rbind(dat0, obs0)
      fit <- sx <- dx <- data.frame(matrix(nrow=n.ahead,ncol=length(nomi)))
      colnames(fit) <- colnames(sx) <- colnames(dx) <- nomi
      for(j in 1:n.ahead) {
        pred <- lapply(object$models, predict, newdata=datOK, interval="prediction",level=level)
        fit[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),1]})
        sx[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),2]})
        dx[j,nomi] <- sapply(pred[nomi],function(x){x[nrow(x),3]})
        newobs <- obs0
        newobs[,nomi] <- fit[j,nomi]
        datOK <- rbind(dat0, newobs, obs0)
        }
      predList <- list()
      for(i in 1:length(nomi)) {
        itab <- data.frame(1:n.ahead,fit[,nomi[i]],sx[,nomi[i]],dx[,nomi[i]])
        colnames(itab) <- c("n.ahead","mean","lower","upper")
        predList[[i]] <- itab
        }
      names(predList) <- nomi
      } else {
      dataList <- split(object$data,object$data[,object$call$unit])
      glev <- levels(object$data[,object$call$unit])
      fitList <- sxList <- dxList <- list()
      for(i in 1:length(dataList)) {
        dat0 <- dataList[[i]]
        obs0 <- dat0[nrow(dat0),]
        datOK <- rbind(dat0, obs0)
        fit <- sx <- dx <- data.frame(matrix(nrow=n.ahead,ncol=length(nomi)))
        colnames(fit) <- colnames(sx) <- colnames(dx) <- nomi
        for(j in 1:n.ahead) {
          predList <- lapply(object$models, predict, newdata=datOK, interval="prediction", level=level)
          fit[j,nomi] <- sapply(predList[nomi],function(x){x[nrow(x),1]})
          sx[j,nomi] <- sapply(predList[nomi],function(x){x[nrow(x),2]})
          dx[j,nomi] <- sapply(predList[nomi],function(x){x[nrow(x),3]})
          newobs <- obs0
          newobs[,nomi] <- fit[j,nomi]
          datOK <- rbind(dat0, newobs, obs0)
          } 
        fitList[[i]] <- cbind(names(dataList)[i],1:n.ahead,fit)
        sxList[[i]] <- cbind(names(dataList)[i],1:n.ahead,sx)
        dxList[[i]] <- cbind(names(dataList)[i],1:n.ahead,dx)
        }
      fitOK <- do.call(rbind,fitList)
      sxOK <- do.call(rbind,sxList)
      dxOK <- do.call(rbind,dxList)
      fitOK[,1] <- factor(fitOK[,1],levels=glev)
      sxOK[,1] <- factor(sxOK[,1],levels=glev)
      dxOK[,1] <- factor(dxOK[,1],levels=glev)
      colnames(fitOK) <- colnames(sxOK) <- colnames(dxOK) <- c(object$call$unit,"n.ahead",nomi)
      predList <- list()
      for(i in 1:length(nomi)) {
        itab <- cbind(fitOK[,1:2],fitOK[,nomi[i]],sxOK[,nomi[i]],dxOK[,nomi[i]])
        colnames(itab) <- c(object$call$unit,"n.ahead","mean","lower","upper")
        predList[[i]] <- itab
        }
      names(predList) <- nomi
      }
    }
  attr(predList,"level") <- level
  predList
  }

# logLik method for class 'feVAR'
logLik.feVAR <- function(object, ...) {
  ll <- sum(sapply(object$models,function(x){logLik(x,...)[1]}))
  npar <- sum(sapply(object$models,function(x){attr(logLik(x,...),"df")}))
  attr(ll,"df") <- npar
  class(ll) <- "logLik"
  ll
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
  n <- sum(sapply(object$models,nobs))
  -2*ll[1]+attr(ll,"df")*log(n)
  }

# residual ACF (auxiliary)
residualACF <- function(x, max.lag=NULL, signif=0.05, print=TRUE, plot=TRUE, ylim=NULL, xlab="lag", ylab="ACF", titles=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), las=0, ...) {
  res <- residuals(x)
  nomi <- x$call$var.names
  if(is.null(x$call$unit)) {
    if(is.null(max.lag)) max.lag <- nrow(res)-3
    mat <- matrix(nrow=max.lag+1,ncol=length(nomi))
    mat[1,] <- 1
    for(i in 1:length(nomi)) {
      for(j in 1:max.lag) {
        ires <- res[,nomi[i]]
        ires_lag <- c(rep(NA,j),ires)[1:length(ires)]
        mat[j+1,i] <- cor(ires,ires_lag,use="pairwise") 
        }
      }  
    } else {
    resList <- split(res[,nomi], x$data[,x$call$unit])
    if(is.null(max.lag)) max.lag <- max(sapply(resList,nrow)-3)
    mat <- matrix(nrow=max.lag+1,ncol=length(nomi))
    mat[1,] <- 1
    for(i in 1:length(nomi)) {
      iresList <- lapply(resList,function(z){z[,nomi[i]]})
      ires <- do.call(c,iresList)
      for(j in 1:max.lag) {
        ires_lag <- do.call(c,lapply(iresList,function(z){
          c(rep(NA,j),z)[1:length(z)]
          }))
        mat[j+1,i] <- cor(ires,ires_lag,use="pairwise") 
        }
      }
    }
  colnames(mat) <- nomi
  rownames(mat) <- 0:max.lag
  matOK <- na.omit(mat)
  attributes(matOK)$na.action <- NULL
  lagOK <- nrow(matOK)-1
  if(plot) {
    if(is.null(titles)) titles <- nomi
    if(is.null(ylim)) ylim <- range(matOK)
    zval <- qnorm(1-signif/2)
    nt <- sapply(resList, nrow)
    opar <- par(no.readonly=T)
    on.exit(par(opar))
    par(mfrow=n2mfrow(length(nomi)), mar=mar, mgp=mgp, las=las)
    for(i in 1:length(nomi)) {
      plot(0:lagOK, matOK[,nomi[i]], type="h", xlab=xlab, ylab=ylab, ylim=ylim, main=titles[i], ...)
      abline(h=0)
      xseq <- seq(0,lagOK,length=100)
      bsx <- sapply(xseq, function(x){-zval/sqrt(sum(nt-x))})
      bdx <- sapply(xseq, function(x){zval/sqrt(sum(nt-x))})
      lines(xseq, bsx, lty=2, col=2)
      lines(xseq, bdx, lty=2, col=2)
      }
    }
  if(print) matOK
  }

# residual qqnorm (auxiliary)
residualQQ <- function(x, xlab="expected", ylab="observed", cex=0.6, titles=NULL, las=0, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  res <- residuals(x)
  nomi <- x$call$var.names
  if(is.null(titles)) titles <- nomi
  opar <- par(no.readonly=T)
  on.exit(par(opar))
  par(mfrow=n2mfrow(length(nomi)), mar=mar, mgp=mgp, las=las)
  for(i in 1:length(nomi)) {
    qqnorm(res[,nomi[i]], xlab=xlab, ylab=ylab, cex=cex, main=titles[i], ...)
    qqline(res[,nomi[i]])
    }
  }

# fitted vs residuals (auxiliary)
fittedVSresiduals <- function(x, xlab="fitted values", ylab="residuals", cex=0.8, titles=NULL, add.grid=TRUE, las=0, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  res <- residuals(x)
  fit <- fitted.values(x)
  nomi <- x$call$var.names
  if(is.null(titles)) titles <- nomi
  opar <- par(no.readonly=T)
  on.exit(par(opar))
  par(mfrow=n2mfrow(length(nomi)), mar=mar, mgp=mgp, las=las)
  for(i in 1:length(nomi)) {
    plot(fit[,nomi[i]], res[,nomi[i]], type="n", xlab=xlab, ylab=ylab, cex=cex, main=titles[i], ...)
    if(add.grid) grid()
    points(fit[,nomi[i]], res[,nomi[i]])
    abline(h=0)
    box()
    }
  }

# spaghetti plot of residuals (auxiliary)
spagResid <- function(x, xlab="time point", ylab="residuals", ylim=NULL, titles=NULL, add.grid=TRUE, las=0, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  res <- residuals(x)
  nomi <- x$call$var.names
  nt <- length(unique(x$data[,x$call$time]))
  if(is.null(titles)) titles <- nomi
  if(is.null(x$call$unit)) {
    opar <- par(no.readonly=T)
    on.exit(par(opar))
    par(mfrow=n2mfrow(length(nomi)), mar=mar, mgp=mgp, las=las)
    for(i in 1:length(nomi)) {
      if(is.null(ylim)) iylim <- range(res[,nomi[i]],na.rm=T)
      plot(res[,nomi[i]], xlim=c(1,nt), ylim=iylim, type="n", xlab=xlab, ylab=ylab, main=titles[i], ...)
      if(add.grid) grid()
      lines(res[,nomi[i]], col="grey70")
      abline(h=0)
      box()
      }
    } else {
    resList <- split(res,x$data[,x$call$unit])
    opar <- par(no.readonly=T)
    on.exit(par(opar))
    par(mfrow=n2mfrow(length(nomi)), mar=mar, mgp=mgp, las=las)
    for(i in 1:length(nomi)) {
      if(is.null(ylim)) iylim <- range(res[,nomi[i]],na.rm=T)
      plot(resList[[1]][,nomi[i]], xlim=c(1,nt), ylim=iylim, type="n", xlab=xlab, ylab=ylab, main=titles[i], ...)
      if(add.grid) grid()
      for(j in 1:length(resList)) {
        lines(resList[[j]][,nomi[i]], col="grey70")
        }
      abline(h=0)
      box()
      }
    }
  }

# graphical diagnostics of residuals
residualPlot <- function(x, type="ts", max.lag=NULL, signif=0.05, acf.print=TRUE, acf.plot=TRUE,
  cex=0.6, xlab=NULL, ylab=NULL, ylim=NULL, titles=NULL, add.grid=TRUE, las=0, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  #
  ## CHECK ARGUMENTS: max.lag=NULL
  #
  if(!identical(class(x),"feVAR")) stop("Argument 'x' must be an object of class 'feVAR'")
  add.grid <- add.grid[1]
  if(is.na(add.grid)||(!is.logical(add.grid)|is.null(add.grid))) add.grid <- TRUE 
  if(is.null(type)||is.na(type)) type <- 1
  if(type[1]=="ts"|type[1]==1) {
    if(is.null(xlab)) xlab <- "time point"
    if(is.null(ylab)) ylab <- "residuals"
    spagResid(x, xlab=xlab, ylab=ylab, ylim=ylim, titles=titles, add.grid=add.grid, las=las, mar=mar, mgp=mgp, ...)
    } else if(type[1]=="acf"|type[1]==2) {
    if(!is.numeric(signif)) signif <- 0.05 else signif <- min(c(1,max(c(signif,0),na.rm=T)),na.rm=T)
    acf.print <- acf.print[1]
    if(is.na(acf.print)||(!is.logical(acf.print)|is.null(acf.print))) acf.print <- TRUE 
    acf.plot <- acf.plot[1]
    if(is.na(acf.plot)||(!is.logical(acf.plot)|is.null(acf.plot))) acf.plot <- TRUE 
    if(is.null(xlab)) xlab <- "lag"
    if(is.null(ylab)) ylab <- "ACF"
    residualACF(x, max.lag=max.lag, signif=signif, ylim=ylim, xlab=xlab, ylab=ylab, titles=titles, las=las, mar=mar, mgp=mgp, print=acf.print, plot=acf.plot, ...)
    } else if(type[1]=="qq"|type[1]==3) {
    if(is.null(xlab)) xlab <- "expected"
    if(is.null(ylab)) ylab <- "observed"
    residualQQ(x, xlab=xlab, ylab=ylab, cex=cex, titles=titles, las=las, mar=mar, mgp=mgp, ...)
    } else if(type[1]=="fitVSres"|type[1]==4) {
    if(is.null(xlab)) xlab <- "fitted values"
    if(is.null(ylab)) ylab <- "residuals"
    fittedVSresiduals(x, xlab=xlab, ylab=ylab, cex=cex, titles=titles, add.grid=add.grid, las=las, mar=mar, mgp=mgp, ...)
    } else {
    stop("Argument 'type' must be one among 'ts' (1), 'acf' (2), 'qq' (3) and 'fitVSres' (4)")
    }
  }

# plot of prediction
predictPlot <- function(x, unit, n.ahead=0, newdata=NULL, cex=0.6, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, titles=NULL, add.grid=TRUE, in.col="grey40", out.col="red", in.lty=1, out.lty=2, bands.col="grey80", cex.axis=NULL, las=NULL, mar=c(3.5,3.5,2,2), mgp=c(2.3,0.8,0), ...) {
  if(!identical(class(x),"feVAR")) stop("Argument 'x' must be an object of class 'feVAR'")
  if(!is.numeric(n.ahead)) n.ahead <- 0 else n.ahead <- round(max(n.ahead,na.rm=T))
  if(is.na(n.ahead)|n.ahead<0|abs(n.ahead)=="Inf") n.ahead <- 0
  if(missing(unit)) {
    stop("Argument 'unit' is missing")
    } else {
    unit <- unit[1]
    if(is.null(unit)||is.na(unit)) stop("Argument 'unit' is missing")
    }
  glev <- levels(x$data[,x$call$unit])
  if((unit%in%glev)==F) stop("Unknown unit '",unit,"'")
  add.grid <- add.grid[1]
  if(is.na(add.grid)||(!is.logical(add.grid)|is.null(add.grid))) add.grid <- TRUE
  #
  ## CHECK ARGUMENTS: newdata=NULL
  #
  nomi <- x$call$var.names
  opar <- par(no.readonly=T)  
  if(is.null(x$call$unit)) {
    insam <- predict.feVAR(x)
    if(n.ahead>0) {
      outsam <- predict.feVAR(x, n.ahead=n.ahead)
      outnam <- paste0("(",1:n.ahead,")")
      } else {
      outsam <- outnam <- c()
      }
    dat <- x$data
    datnew <- newdata
    } else {
    insam <- predict.feVAR(x, unit=unit)
    if(n.ahead>0) {
      outsam <- predict.feVAR(x, n.ahead=n.ahead, unit=unit)
      outnam <- paste0("(",1:n.ahead,")")
      } else {
      outsam <- outnam <- c()
      }
    dat <- x$data[which(x$data[,x$call$unit] %in% unit),]
    datnew <- newdata[which(newdata[,x$call$unit] %in% unit),]
    }
  n <- nrow(dat)
  if(!is.null(datnew)&&nrow(datnew)==0) datnew <- NULL
  if(is.null(x$call$time)) {
    tnam <- 1:n
    } else {
    tnam <- sort(unique(dat[,x$call$time]))
    }
  if(is.null(titles)) titles <- nomi
  if(is.null(xlab)) xlab <- ""
  if(is.null(ylab)) ylab <- "" 
  if(is.null(cex.axis)|!is.numeric(cex.axis)) cex.axis <- c(1,1)
  if(length(cex.axis)<2) cex.axis <- rep(cex.axis,2)
  if(is.null(las)) {
    if(is.null(x$call$time)) las <- 1 else las <- 2
    }
  on.exit(opar)
  par(mfrow=n2mfrow(length(nomi)), mar=mar, mgp=mgp, las=las)
  for(i in 1:length(nomi)) {
    iobs <- dat[,nomi[i]]
    idat <- rbind(insam[[nomi[i]]][,c("mean","lower","upper")],
                  outsam[[nomi[i]]][,c("mean","lower","upper")])
    iN <- nrow(idat)
    if(is.null(ylim)) {
      ilimy <- range(c(iobs,idat,datnew[,nomi[i]]),na.rm=T)
      } else {
      ilimy <- ylim
      }
    if(is.null(xlim)) {
      ilimx <- c(1,iN)
      } else {
      ilimx <- xlim
      }
    plot(idat[,1], type="n", xlab=xlab, ylab=ylab, las=las, ylim=ilimy, xlim=ilimx, main=titles[i],
         cex.axis=cex.axis[2], xaxt="n", ...)
    if(add.grid) grid()
    polygon(c(1:iN,iN:1), c(idat[,2],rev(idat[,3])),
            border=NA, col=adjustcolor(bands.col,alpha.f=0.5))
    if(n.ahead>0) lines(n:iN, c(iobs[n],idat[(n+1):iN,1]), col=out.col, lty=out.lty)
    points(1:n, iobs, cex=cex, col=in.col)
    lines(1:n, iobs, col=in.col, lty=in.lty)
    if(!is.null(datnew)) {
      points((n+1):(n+nrow(datnew)), datnew[,nomi[i]], cex=cex, col=in.col)
      lines(n:(n+nrow(datnew)), c(iobs[n],datnew[,nomi[i]]), col=in.col, lty=in.lty)
      }
    axis(1, at=1:iN, labels=c(tnam,outnam), las=las, cex.axis=cex.axis[1])
    box()
    }
  }
