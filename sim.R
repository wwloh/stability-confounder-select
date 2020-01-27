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

# initialize for parallel MC jobs
args <- 1
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'  
}
(seed <- as.integer(args[1]))

simsets <- expand.grid(n=80,p=25,q=c(1,1.8))
if (!grepl("apple",sessionInfo()[[1]]$platform)) {
  simsets <- simsets[rep(1:nrow(simsets),each=400),]
}

(n <- simsets[seed,"n"]) # sample size
(p <- simsets[seed,"p"]) # number of covariates
(q <- simsets[seed,"q"]) # coefficient of confounder in PS model
rm(simsets)

# indices for (C)onfounders, (I)nstruments and outcome (P)redictors
lC <- 1:2
lI <- 5:6
lP <- 3:4

# conditional probability of treatment: true propensity score
ProbTreat <- function(l) {
  gamma0 <- 0
  g_raw <- rep(0,p)
  g_raw[lC] <- q # true confounders
  g_raw[lI] <- 1.8 # instruments
  gamma <- g_raw
  linear.ps <- gamma0 + (l %*% gamma)[,1]
  pA1 <- exp(linear.ps)/(1+exp(linear.ps))
  return(pA1)
}

# uniformity trial outcomes
Y0out <- function(l) {
  beta0 <- 0
  b_raw <- rep(0,p)
  b_raw[lC] <- 0.6 # true confounders
  b_raw[lP] <- 0.6
  beta <- b_raw
  linear.out <- beta0 + (l %*% beta)[,1]
  y0_true <- linear.out + rnorm(n=n,sd=4)
  return(y0_true)
}

One_ObsData <- function() {
  l <- matrix(rnorm(n=n*p),nrow=n,ncol=p)
  # Normlize covariates to have mean 0 and standard deviation 1
  l <- scale(l,center=TRUE,scale=TRUE)
  pA1 <- ProbTreat(l=l)
  # Bernoulli sampling
  a <- rbinom(n,1,pA1)
  # table(a);boxplot(pA1~a)
  y <- Y0out(l=l)
  # y <- y + q*a # ATE
  mydata <- data.frame("i"=1:n,"L"=l,"pA1.true"=pA1,"treat"=a,"Y"=y)
  return(mydata)
}

## full PS model with all covariates
ps.full <- as.formula(paste0("treat~",
                             paste(paste0("L.",1:p),collapse="+")))
## true PS model with only confounders
ps.true <- as.formula(paste0("treat~",
                             paste(paste0("L.",lC),collapse="+")))
## target PS model with all covariates excluding instruments
ps.targ <- as.formula(paste0("treat~",
                             paste(paste0("L.",c(lC,lP)),collapse="+")))

# load helper functions
source("fullmatch-funs.R")
source("oal.R")
source("other_selection_methods.R")

One_sim <- function(M=1e2) {
  check_positivity <- FALSE
  # generate a dataset where the full PS model can be fitted with no issues
  while(!check_positivity) {
    mydata <- One_ObsData()
    check_glm.try <- tryCatch(
      check_glm <- glm(ps.full,family=binomial(link = "logit"),data=mydata),
      warning=function(cond) return(NA))
    if (any(is.na(check_glm.try))) {
      check_positivity <- FALSE
    } else {
      fit_glm <- predict.glm(check_glm,type="response")
      m.out <- matchit(ps.full, data = mydata, method = "full")
      check_positivity <- 
        # fitted model converged
        check_glm$converged & 
        # no perfect separation
        all(pmin(fit_glm,(1-fit_glm))>.Machine$double.eps*1e1) &
        sum(m.out$discarded)==0  
    }
  }
  res <- list()
  # parametric bootstrap with full PS model
  pv <- One_pv_paraboot(mydata,ps_fit=ps.full,n_resample=0)
  pv <- pv["gee"] # keep GEE p-value only
  res <- c(res,pv); rm(pv)
  # parametric bootstrap with target PS model
  pv <- One_pv_paraboot(mydata,ps_fit=ps.targ,n_resample=0)
  pv <- pv["gee"] # keep GEE p-value only
  names(pv) <- paste0(names(pv),".target")
  res <- c(res,pv); rm(pv)
  
  # full matching with full PS model
  pv.cond <- One_pv_strata(mydata,ps_fit=ps.full,n_resample=M)
  res <- c(res,pv.cond); rm(pv.cond)
  # full matching with target PS model
  pv.cond <- One_pv_strata(mydata,ps_fit=ps.targ,n_resample=M)
  names(pv.cond) <- "conditional.target"
  res <- c(res,pv.cond); rm(pv.cond)
  # full matching with true PS
  pv.cond <- One_pv_strata(mydata,ps_fit=mydata$pA1.true,n_resample=M)
  names(pv.cond) <- "conditional.oraclePS"
  res <- c(res,pv.cond); rm(pv.cond)
  # full matching with empty PS model
  data.noL <- cbind(mydata[,c("i","treat","Y")],"stratum"=1)
  pv.noL <- One_pv_conditional(mydata=data.noL,n_resample=M)
  names(pv.noL) <- "conditional.noL"
  res <- c(res,pv.noL); rm(pv.noL)
  
  # OAL =======================================================================
  var.list <- paste0("L.",1:p)
  Data <- mydata
  colnames(Data)[colnames(Data)=="treat"] <- "A"
  res.oal <- OneData_OAL(var.list,Data)
  ps.oal <- res.oal$Data[,paste("f.pA",names(res.oal$tt),sep="")]
  # parametric bootstrap with PS from OAL
  pv.oal <- One_pv_paraboot_OAL(mydata,ps_fit=ps.oal,n_resample=0)
  pv.oal <- pv.oal["gee"] # keep GEE p-value only
  names(pv.oal) <- paste0(names(pv.oal),".OAL")
  res <- c(res,pv.oal); rm(pv.oal)
  # full matching with PS from OAL
  pv.oal <- One_pv_strata(mydata,ps_fit=ps.oal,n_resample=M)
  names(pv.oal) <- "conditional.OAL"
  res <- c(res,pv.oal); rm(pv.oal)
  
  # list for storing all ordered covariates
  res.FwdSel <- list() 
  # order covariates based on 'priority' to be confounders ====================
  res.FwdSel[["minimin"]] <- ForwardSelect_minimax(
    Y=mydata$Y,X=mydata[,grepl("L",colnames(mydata))],Z=mydata$treat,
    criterion="minimin")
  
  # other methods =============================================================
  ## posterior p-values
  out_BCEE <- OneData_BCEE(mydata,var.list)
  pv.bcee <- out_BCEE[[1]]
  names(pv.bcee) <- "posterior.BCEE"
  res <- c(res,pv.bcee); rm(pv.bcee)
  
  out_bacr <- OneData_bacr(mydata,var.list)
  pv.bacr <- out_bacr[[1]]
  names(pv.bacr) <- "posterior.bacr"
  res <- c(res,pv.bacr); rm(pv.bacr)
  
  ## full matching using selected covariates
  L.selected <- list()
  L.selected[["Boruta"]] <- OneData_Boruta(mydata,var.list)
  L.selected[["CovSel"]] <- OneData_CovSel(mydata,var.list)
  # functions that require data.frame in global environment
  L.selected[["ABE"]] <- OneData_ABE(mydata,var.list)
  rm(data.ABE,fitY.ABE)
  L.selected[["SignifReg"]]  <- OneData_SignifReg(mydata,var.list)
  rm(data.SignifReg,fitY.SignifReg)
  
  # how many confounders selected?
  PSsizes <- list()
  # how many true confounders selected?
  LC.selected <- list()
  # point estimates of ATE
  ATE.selected <- list()
  ATE.selected[["OAL"]] <- res.oal$ATE
  
  for (ll in 1:length(L.selected)) {
    l.select <- L.selected[[ll]]
    if (any(is.na(l.select))) {
      pv.select <- NA
      PSsizes[[ll]] <- NA
      LC.selected[[ll]] <- NA
      ATE.selected[[ll]] <- NA
    } else {
      ps.select <- as.formula(paste0("treat~",paste(l.select,collapse="+")))
      pv.select <- One_pv_strata(mydata,ps_fit=ps.select,n_resample=M)
      PSsizes[[ll]] <- sum(grepl("L",l.select))
      LC.selected[[ll]] <- sum(paste0("L.",lC) %in% l.select)
      ## outcome model conditional on treatment and selected covariates
      ATE.selected[[ll]] <- ATEhat_lm(mydata,Ls=l.select,point.only=FALSE)
    }
    names(pv.select) <- paste0(names(pv.select),".",names(L.selected)[ll])
    res <- c(res,pv.select)
    rm(l.select,ps.select,pv.select)
  }
  names(PSsizes) <- names(L.selected)
  names(LC.selected) <- names(L.selected)
  names(ATE.selected) <- names(L.selected)
  
  ## p-values from CTMLE (without and with LASSO)
  pv.ctmle.try <- tryCatch(
    pv.ctmle <- OneData_ctmle(mydata,var.list),
    error=function(cond) return(NA))
  names.ctmle <- c("nolasso","lasso1","lasso2")
  if (any(is.na(pv.ctmle.try))) {
    ctmle_NAs <- rep(NA,3)
    names(ctmle_NAs) <- names.ctmle
    pv.ctmle <- ctmle_NAs
    ATE.selected[["ctmle"]] <- rep(NA,6)
  } else {
    ATE.selected[["ctmle"]] <- c(pv.ctmle[,"est"],pv.ctmle[,"se"])
    pv.ctmle <- pv.ctmle[,"pv"]
  }
  names(ATE.selected[["ctmle"]]) <- c(paste0(names.ctmle,".est"),
                                      paste0(names.ctmle,".se"))
  names(pv.ctmle) <- paste0("ctmle.",names(pv.ctmle))
  res <- c(res,pv.ctmle); rm(pv.ctmle)
  
  ## p-values from HDM
  pv.hdm <- OneData_hdm(mydata,var.list)
  pv.hdm.list <- lapply(1:nrow(pv.hdm), function(x) pv.hdm[x,c("est","se")])
  names(pv.hdm.list) <- rownames(pv.hdm)
  ATE.selected[["hdm"]] <- unlist(pv.hdm.list)
  pv.hdm <- pv.hdm[,"pv"]
  names(pv.hdm) <- paste0("hdm.",names(pv.hdm))
  res <- c(res,pv.hdm); rm(pv.hdm)
  
  ATE.selected <- unlist(ATE.selected)
  
  # diagnostics over sequence of submodels (excluding empty submodel) =========
  diags.FwdSel <- lapply(res.FwdSel, function(res.fs) {
    DiagnosticsCovariatesOrdered(L.ordered=res.fs[[1]],mydata,M)
  })
  
  # p-values
  pv.FwdSel <- lapply(diags.FwdSel, function(diags.fs) {
    pv.fs <- diags.fs[[1]][diags.fs[[2]],1]
    names(pv.fs) <- paste0("pv.",names(diags.fs[[2]]))
    return(pv.fs)
  })
  
  # summaries of selected PS model using Forward Selection
  PSsizes.FwdSel <- list()
  LC.selected.FwdSel <- list()
  ATE.selected.FwdSel <- list()
  for (meth in 1:length(res.FwdSel)) {
    # how many confounders selected?
    PSsizes.FwdSel[[meth]] <- diags.FwdSel[[meth]][[2]]["score.vs.allL.Q.3"]
    # orders of true confounders in ordered priorities
    LC.order <- OrderCovariateIndices(res.FwdSel[[meth]][[1]])[lC,2]
    # how many true confounders selected?
    LC.selected.FwdSel[[meth]] <- sum(sapply(LC.order, function(lc)
      lc <= PSsizes.FwdSel[[meth]]))
    # point estimate of ATE
    ATE.selected.FwdSel[[meth]] <- diags.FwdSel[[meth]][[1]][,c(
      "treat","treat.se.lm","treat.se.score")][PSsizes.FwdSel[[meth]],]
  }
  names(PSsizes.FwdSel) <- names(res.FwdSel)
  names(LC.selected.FwdSel) <- names(res.FwdSel)
  names(ATE.selected.FwdSel) <- names(res.FwdSel)

  PSsizes <- c(PSsizes,unlist(PSsizes.FwdSel))
  names(PSsizes) <- paste0("PSsizes.",names(PSsizes))
  LC.selected <- c(LC.selected,unlist(LC.selected.FwdSel))
  names(LC.selected) <- paste0("selected.",names(LC.selected))
  ATE.selected <- c(ATE.selected,unlist(ATE.selected.FwdSel))
  names(ATE.selected) <- paste0("treat.",names(ATE.selected))
  
  return(c(n=n,p=p,q=q,res,unlist(pv.FwdSel),
           PSsizes,LC.selected,ATE.selected))
}

n_sims <- 2e1
ptm=proc.time()[3]
sim_res <- replicate(n=n_sims,expr=One_sim(M=1e3),simplify=FALSE)
proc.time()[3]-ptm
save(sim_res,file=paste0("sim-",seed,".Rdata"))
q()
