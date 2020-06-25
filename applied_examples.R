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

dataset.names <- c("rhc5","rhc7","rhc9","speff2","lalonde")
for (data.name in dataset.names) {
  if (grepl("rhc",data.name)) {
    windowsize <- as.integer(strsplit(data.name,"rhc")[[1]][2])
    if (is.na(windowsize)) windowsize <- 5
    data.name <- "rhc"
  } else {
    windowsize <- 3
  }
  
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
  
  # randomization inference ignoring all covariates
  data.noL <- cbind(mydata[,c("i","treat","Y")],"stratum"=1)
  pv.noL <- One_pv_conditional(mydata=data.noL,n_resample=M)
  names(pv.noL) <- paste0(names(pv.noL),".noL")
  res <- c(res,pv.noL); rm(pv.noL,data.noL)
  # full matching with full PS model
  pv.cond <- One_pv_strata(mydata,ps_fit=ps.full,n_resample=M)
  names(pv.cond) <- paste0(names(pv.cond),".full")
  res <- c(res,pv.cond); rm(pv.cond)
  
  var.list <- paste0("L.",1:p)
  
  # list for storing all ordered covariates
  res.FwdSel <- list()
  # order covariates based on priority to be confounders ======================
  res.FwdSel[["ds_stable"]] <- ForwardSelect_DS(
    Y=mydata$Y,X=mydata[,var.list],A=mydata$treat)
  
  # diagnostics over sequence of submodels (excluding empty submodel) =========
  diags.FwdSel <- lapply(res.FwdSel, function(res.fs) {
    DiagnosticsCovariatesOrdered(L.ordered=res.fs[["ordered"]],mydata,
                                 k=windowsize)
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
  if(data.name!="rhc") {
    L.selected[["CovSel"]] <- OneData_CovSel(mydata,var.list,
                                             cov.sel.type="np",cov.sel.alg=2)
  } else {
    L.selected[["CovSel"]] <- "1"
  }
  # functions that require data.frame in global environment
  L.selected[["ABE"]] <- OneData_ABE(mydata,var.list)
  rm(data.ABE,fitY.ABE)
  L.selected[["SignifReg"]]  <- OneData_SignifReg(mydata,var.list)
  rm(data.SignifReg,fitY.SignifReg)
  
  # empty and full models
  L.selected[["Empty"]] <- "1"
  L.selected[["Full"]] <- var.list
  
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
  
  # ttmt effect estimator given selected covariates
  ATE.selected <- lapply(res.selected, "[", c("ate","se"))
  
  rm(res.selected)
  
  # save results ==============================================================
  if (data.name=="rhc") {
    data.name <- paste0(data.name,windowsize)
  }
  res_filename <- paste0("stability-applied_examples-",data.name,".Rdata")
  save(res,res.FwdSel,diags.FwdSel,PSsizes,L.selected,ATE.selected,
       file=res_filename)
}
q()
