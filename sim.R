rm(list=ls())
libraries_check <- c("data.table","MatchIt","optmatch","geepack","MASS","lqa",
                     "Hmisc")
# install.packages("lqa_1.0-3.tar.gz", repos = NULL, type="source")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# simulation settings
simsets <- expand.grid(n=80,p=c(25,60),z=c(2,4),q=c(0,1),y=c(0,1))

# initialize for parallel MC jobs
args <- nrow(simsets)
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
  simsets <- simsets[rep(1:nrow(simsets),each=100),]
  nrow(simsets)
}
(seed <- as.integer(args[1]))

(n <- simsets[seed,"n"]) # sample size
(p <- simsets[seed,"p"]) # number of covariates
(z <- simsets[seed,"z"]) # number of instruments
(q <- simsets[seed,"q"]) # whether there is collider bias
(y <- simsets[seed,"y"]) # continuous or binary y

# indices for (C)onfounders, (I)nstruments and outcome (P)redictors
lC <- 1:2
lP <- 3:4
lI <- 5+(0:(z-1))

One_ObsData <- function() {
  l <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
  u <- matrix(0,nrow=n,ncol=2)
  if (q>0) {
    # unmeasured confounders to possibly induce collider bias
    u <- matrix(rnorm(n=n*2,sd=1/4),nrow=n,ncol=2)
    # covariates that induce collider bias and should not be adjusted for
    for (k in lI) {
      l[,k] <- rnorm(n=n,mean=2*u[,1]+2*u[,2],sd=sqrt(1/2))
    }
  }
  # round(apply(l,2,var),1); round(colMeans(l),1)
  
  # conditional probability of treatment: true propensity score
  g_raw <- rep(0,p)
  g_raw[lC] <- 1 # true confounders
  g_raw[lI] <- 1.6*(1-q) # instruments when q==0
  nu <- 2*q # to induce collider bias when q==1
  linear.ps <- (l %*% g_raw)[,1] + nu*u[,1]
  pA1 <- exp(linear.ps)/(1+exp(linear.ps))
  rm(linear.ps,g_raw)
  
  # uniformity trial outcomes
  b_raw <- rep(0,p)
  b_raw[lC] <- 0.8 # true confounders
  b_raw[lP] <- 0.8 # predictors of outcome only
  linear.out <- (l %*% b_raw)[,1] + nu*u[,2]
  if (y==0) {
    y0_true <- linear.out + rnorm(n=n,sd=4)
  } else {
    pY1 <- exp(linear.out)/(1+exp(linear.out))
    y0_true <- rbinom(n,1,pY1)
  }
  
  # Bernoulli sampling
  a <- rbinom(n,1,pA1)
  # table(a);boxplot(pA1~a)
  mydata <- data.frame("i"=1:n,"L"=l,"treat"=a,"Y"=y0_true)
  return(mydata)
}

## full PS model with all covariates
ps.full <- as.formula(paste0("treat~",
                             paste(paste0("L.",1:p),collapse="+")))
## target PS model with all predictors of outcome and excluding instruments
ps.targ <- as.formula(paste0("treat~",
                             paste(paste0("L.",c(lC,lP)),collapse="+")))

# load helper functions
source("fullmatch-funs.R")
source("oal.R")
source("other_selection_methods.R")

One_sim <- function(M=1e2) {
  mydata <- One_ObsData()
  
  res <- list()
  
  # randomization inference ignoring all covariates
  data.noL <- cbind(mydata[,c("i","treat","Y")],"stratum"=1)
  pv.noL <- One_pv_conditional(mydata=data.noL,n_resample=M)
  names(pv.noL) <- paste0(names(pv.noL),".noL")
  res <- c(res,pv.noL); rm(pv.noL,data.noL)
  # full matching with target PS model
  pv.cond <- One_pv_strata(mydata,ps_fit=ps.targ,n_resample=M)
  names(pv.cond) <- paste0(names(pv.cond),".target")
  res <- c(res,pv.cond); rm(pv.cond)
  
  # OAL =======================================================================
  var.list <- paste0("L.",1:p)
  Data <- mydata
  colnames(Data)[colnames(Data)=="treat"] <- "A"
  res.oal <- OneData_OAL(var.list,Data)
  ps.oal <- res.oal$Data[,paste("f.pA",names(res.oal$tt),sep="")]
  # Wald test with PS from OAL (GEE p-value only)
  pv.oal <- One_IPW_GEE(mydata,ps_fit=ps.oal)
  names(pv.oal) <- paste0(names(pv.oal),".OAL")
  res <- c(res,pv.oal); rm(pv.oal,ps.oal,res.oal,Data)
  
  # list for storing all ordered covariates
  res.FwdSel <- list()
  # order covariates based on priority to be confounders ======================
  res.FwdSel[["ds_stable"]] <- ForwardSelect_DS(
    Y=mydata$Y,X=mydata[,var.list],A=mydata$treat)
  
  # diagnostics over sequence of submodels (excluding empty submodel) =========
  diags.FwdSel <- lapply(res.FwdSel, function(res.fs) {
    DiagnosticsCovariatesOrdered(L.ordered=res.fs[["ordered"]],mydata,k=5)
  })
  
  # summaries of selected PS model using Forward Selection
  L.selected <- list()
  for (meth in 1:length(res.FwdSel)) {
    # how many confounders selected?
    pssize.meth <- as.integer(diags.FwdSel[[meth]][["selected_orbit"]][1])
    L.selected[[meth]] <- res.FwdSel[[meth]][["ordered"]][1:pssize.meth]
    rm(pssize.meth)
  }
  names(L.selected) <- names(res.FwdSel)
  
  # other covariate selection methods =========================================
  L.selected[["Boruta"]] <- OneData_Boruta(mydata,var.list)
  L.selected[["CovSel"]] <- OneData_CovSel(mydata,var.list)
  # functions that require data.frame in global environment
  L.selected[["ABE"]] <- OneData_ABE(mydata,var.list)
  rm(data.ABE,fitY.ABE)
  L.selected[["SignifReg"]]  <- OneData_SignifReg(mydata,var.list)
  rm(data.SignifReg,fitY.SignifReg)
  
  res.selected <- lapply(L.selected, function(l.select) {
    if (any(!grepl("L",l.select) | is.na(l.select))) {
      # randomization inference ignoring all covariates
      data.noL <- cbind(mydata[,c("i","treat","Y")],"stratum"=1)
      pv.select <- One_pv_conditional(mydata=data.noL,n_resample=M)
      rm(data.noL)
    } else {
      ## full matching using selected covariates
      ps.select <- as.formula(paste0("treat~",paste(l.select,collapse="+")))
      pv.select <- One_pv_strata(mydata,ps_fit=ps.select,n_resample=M)
      rm(ps.select)
    }
    return(unlist(list(
      pv.select,
      "pssize"=sum(grepl("L",l.select)),
      "lcselect"=sum(paste0("L.",lC) %in% l.select),
      OneDRstd_Est(Ls=l.select,mydata=mydata,varest=TRUE)
    )))
  })
  
  pv.selected <- unlist(lapply(res.selected, "[", "conditional"))
  names(pv.selected) <- paste0("conditional.",names(res.selected))
  res <- c(res,pv.selected)
  rm(pv.selected)
  
  # how many confounders selected?
  PSsizes <- unlist(lapply(res.selected, "[", "pssize"))
  names(PSsizes) <- unlist(strsplit(names(PSsizes),split=".pssize"))
  names(PSsizes) <- paste0("PSsizes.",names(PSsizes))
  # how many true confounders selected?
  LC.selected <- unlist(lapply(res.selected, "[", "lcselect"))
  names(LC.selected) <- unlist(strsplit(names(LC.selected),split=".lcselect"))
  names(LC.selected) <- paste0("selected.",names(LC.selected))
  
  # ttmt effect estimator given selected covariates
  ATE.selected <- lapply(res.selected, "[", c("ate","se"))
  
  rm(res.selected)
  
  ## posterior p-values
  if (p<50) {
    out_BCEE <- OneData_BCEE(mydata,var.list)
    pv.bcee <- out_BCEE[[1]]  
    rm(out_BCEE)
    
    out_bacr <- OneData_bacr(mydata,var.list)
    pv.bacr <- out_bacr[[1]]
    rm(out_bacr)
  } else {
    pv.bcee <- NA
    pv.bacr <- NA
  }
  names(pv.bcee) <- "posterior.BCEE"
  res <- c(res,pv.bcee); rm(pv.bcee)
  names(pv.bacr) <- "posterior.bacr"
  res <- c(res,pv.bacr); rm(pv.bacr)
  
  ## p-values from HDM
  pv.hdm <- OneData_hdm(mydata,var.list)
  pv.hdm.list <- lapply(1:nrow(pv.hdm), function(x) pv.hdm[x,c("est","se")])
  names(pv.hdm.list) <- rownames(pv.hdm)
  ATE.selected[["hdm"]] <- unlist(pv.hdm.list)
  pv.hdm <- pv.hdm[,"pv"]
  names(pv.hdm) <- paste0("hdm.",names(pv.hdm))
  res <- c(res,pv.hdm); rm(pv.hdm,pv.hdm.list)
  
  if (p<50) {
    ## p-values from CTMLE (without and with LASSO)
    pv.ctmle.try <- tryCatch(
      pv.ctmle <- OneData_ctmle(mydata,var.list),
      error=function(cond) return(NA))
  } else {
    pv.ctmle.try <- NA
  }
  names.ctmle <- c("nolasso","lasso")
  if (any(is.na(pv.ctmle.try))) {
    pv.ctmle.list  <- lapply(1:length(names.ctmle), function(x) 
      c("est"=NA,"se"=NA))
    pv.ctmle <- rep(NA,length(names.ctmle))
  } else {
    pv.ctmle.list <- lapply(1:nrow(pv.ctmle), function(x) 
      pv.ctmle[x,c("est","se")])
    pv.ctmle <- pv.ctmle[,"pv"]
  }
  names(pv.ctmle.list) <- names.ctmle
  ATE.selected[["ctmle"]] <- unlist(pv.ctmle.list)
  names(pv.ctmle) <- paste0("ctmle.",names.ctmle)
  res <- c(res,pv.ctmle); rm(pv.ctmle)
  
  ATE.selected <- unlist(ATE.selected)
  names(ATE.selected) <- paste0("treat.",names(ATE.selected))
  
  return( c(unlist(simsets[seed,]),res,PSsizes,LC.selected,ATE.selected) )
}

n_sims <- 2e1
ptm=proc.time()[3]
sim_res <- replicate(n=n_sims,expr=One_sim(M=1e3),simplify=FALSE)
proc.time()[3]-ptm
save(sim_res,file=paste0("stability-sim-",seed,".Rdata"))
q()
