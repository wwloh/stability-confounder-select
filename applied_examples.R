rm(list=ls())
libraries_check <- c("data.table","MatchIt","optmatch","geepack","MASS","lqa",
                     "Hmisc","CovSel","speff2trial","locfit")
# install.packages("lqa_1.0-3.tar.gz", repos = NULL, type="source")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# load helper functions
source("fullmatch-funs.R")
source("oal.R")
source("other_selection_methods.R")

dataset.names <- c("lalonde","speff2","rhc")
for (data.name in dataset.names) {
  source(paste0("data-prep-",data.name,".R"))
  
  # mydata: data.frame with 
  ## "i"=observation number
  ## "L"=covariates labelled L.1, ... L.p
  ## "treat"=observed treatment
  ## "Y"=observed outcome

  p <- sum(grepl("L",colnames(mydata)))
  n <- nrow(mydata)
  M <- 2e3

  ## full PS model with all covariates
  ps.full <- as.formula(paste0("treat~",paste(paste0("L.",1:p),collapse="+")))
  
  res <- list()
  # parametric bootstrap with full PS model
  pv <- One_pv_paraboot(mydata,ps_fit=ps.full,n_resample=0)
  pv <- pv["gee"] # keep GEE p-value only
  res <- c(res,pv); rm(pv)
  
  # full matching with full PS model
  pv.cond <- One_pv_strata(mydata,ps_fit=ps.full,n_resample=M)
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
  ## full matching using selected covariates
  L.selected <- list()
  L.selected[["Boruta"]] <- OneData_Boruta(mydata,var.list)
  if(data.name!="rhc") {
    L.selected[["CovSel"]] <- OneData_CovSel(mydata,var.list,
                                             cov.sel.type="np",cov.sel.alg=2)
  } else {
    L.selected[["CovSel"]] <- NA
  }
  # functions that require data.frame in global environment
  L.selected[["ABE"]] <- OneData_ABE(mydata,var.list)
  rm(data.ABE,fitY.ABE)
  L.selected[["SignifReg"]]  <- OneData_SignifReg(mydata,var.list)
  rm(data.SignifReg,fitY.SignifReg)
  
  # how many confounders selected?
  PSsizes <- list()
  # point and standard error estimates of ATE
  ATE.selected <- list()
  
  for (ll in 1:length(L.selected)) {
    l.select <- L.selected[[ll]]
    if (any(is.na(l.select))) {
      pv.select <- NA
      PSsizes[[ll]] <- NA
      ATE.selected[[ll]] <- NA
    } else {
      ps.select <- as.formula(paste0("treat~",paste(l.select,collapse="+")))
      pv.select <- One_pv_strata(mydata,ps_fit=ps.select,n_resample=M)
      PSsizes[[ll]] <- sum(grepl("L",l.select))
      ATE.selected[[ll]] <- rbind(
        ## outcome model conditional on treatment and selected covariates
        "lm"=ATEhat_lm(mydata,Ls=l.select,point.only=FALSE)[1:2],
        ## IPW estimator with sandwich standard error estimate
        "gee"=ATEhat_gee(mydata,ps_fit=ps.select))
      colnames(ATE.selected[[ll]]) <- c("est","se")
    }
    names(pv.select) <- paste0(names(pv.select),".",names(L.selected)[ll])
    res <- c(res,pv.select)
    rm(l.select,ps.select,pv.select)
  }
  names(PSsizes) <- names(L.selected)
  names(ATE.selected) <- names(L.selected)
  
  ## CTMLE (without and with LASSO)
  pv.ctmle.try <- tryCatch(
    pv.ctmle <- OneData_ctmle(mydata,var.list,pvalues=FALSE),
    error=function(cond) return(NA))
  if (any(is.na(pv.ctmle.try))) {
    ctmle_NAs <- matrix(NA,nrow=3,ncol=2)
    rownames(ctmle_NAs) <- c("nolasso","lasso1","lasso2")
    ATE.selected[["ctmle"]] <- ctmle_NAs
  } else {
    ATE.selected[["ctmle"]] <- pv.ctmle[,c("est","se")]
  }
  colnames(ATE.selected[["ctmle"]]) <- names(ATE.selected[[1]])[1:2]
  
  ## HDM
  pv.hdm <- OneData_hdm(mydata,var.list,pvalues=TRUE)
  ATE.selected[["hdm"]] <- pv.hdm[,c("est","se")]
  pv.hdm <- pv.hdm[,"pv"]
  names(pv.hdm) <- paste0("hdm.",names(pv.hdm))
  res <- c(res,pv.hdm); rm(pv.hdm)
  
  # IPW estimator with sandwich standard error estimate
  ## empty PS model
  ATE.selected[["empty"]] <- ATEhat_gee(mydata,ps_fit=as.formula("treat~1"))
  colnames(ATE.selected[["empty"]]) <- colnames(ATE.selected[[1]])
  ## full PS model
  ATE.selected[["full"]] <- ATEhat_gee(mydata,ps_fit=ps.full)
  colnames(ATE.selected[["full"]]) <- colnames(ATE.selected[[1]])
  
  ATE.selected <- do.call(rbind,ATE.selected)
  ctmle.idx <- grep("ctmle",rownames(ATE.selected))
  hdm.idx <- grep("hdm",rownames(ATE.selected))
  rownames(ATE.selected)[ctmle.idx] <- paste0(
    "CTMLE \n (",rownames(ATE.selected)[ctmle.idx], ")")
  rownames(ATE.selected)[hdm.idx] <- paste0(
    "HDM \n (",rownames(ATE.selected)[hdm.idx], ")")
  ATE.selected <- ATE.selected[order(rownames(ATE.selected)),]
  
  # specific set of covariates to adjust for
  L.selected <- sapply(L.selected, function(l.select) {
    if(any(is.na(l.select))) {
      return(l.select)
    } else {
      return(Lnames[sort(na.omit(as.integer(unlist(strsplit(
        l.select,"L."))))),])  
    }
    
  }, simplify=FALSE)
  
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
  L.selected.FwdSel <- list()
  ATE.selected.FwdSel <- list()
  for (meth in 1:length(res.FwdSel)) {
    # how many confounders selected?
    PSsizes.FwdSel[[meth]] <- diags.FwdSel[[meth]][[2]]["score.vs.allL.Q.3"]
    
    # specific set of covariates to adjust for
    L.selected.FwdSel[[meth]] <- sapply(PSsizes.FwdSel[[meth]], function(pssize) {
      l.selected <- res.FwdSel[[meth]][[1]][1:pssize]
      return(Lnames[sort(na.omit(as.integer(unlist(strsplit(
        l.selected,"L."))))),])
    }, simplify=FALSE)
    
    # ATE from outcome model conditional on treatment and selected covariates
    ATE.selected.FwdSel[[meth]] <- diags.FwdSel[[meth]][[1]][
      ,c("treat","treat.se.lm")]
    
    # nested PS models (based on ordered covariates)
    res.nested <- sapply(1:length(res.FwdSel[[meth]][[1]]), function(pssize) {
      l.nested <- res.FwdSel[[meth]][[1]][1:pssize]
      ps.nested <- as.formula(paste0("treat~",paste(l.nested,collapse="+")))
      ## IPW estimator with sandwich standard error estimate
      ATEhat_gee(mydata,ps_fit=ps.nested)
    })
    ATE.selected.FwdSel[[meth]] <- cbind(ATE.selected.FwdSel[[meth]],
                                         t(res.nested))
    
  }
  names(PSsizes.FwdSel) <- names(res.FwdSel)
  names(L.selected.FwdSel) <- names(res.FwdSel)
  names(ATE.selected.FwdSel) <- names(res.FwdSel)
  
  # save results ==============================================================s
  PSsizes.FwdSel <- c(unlist(PSsizes.FwdSel),PSsizes)
  names(PSsizes.FwdSel) <- paste0("PSsizes.",names(PSsizes.FwdSel))
  
  L.selected.FwdSel <- c(L.selected.FwdSel,L.selected)
  ATE.selected.FwdSel <- c(ATE.selected.FwdSel,"others"=list(ATE.selected))
  
  res <- c(res,unlist(pv.FwdSel))
  res_filename <- paste0("data-",data.name,".Rdata")
  save(res,res.FwdSel,diags.FwdSel,
       PSsizes.FwdSel,
       L.selected.FwdSel,
       ATE.selected.FwdSel,
       file=res_filename)
}
q()
