library("Boruta")
library("CovSel")
library("ctmle")
# install.packages("BayesPen_1.0.tar.gz", repos = NULL, type="source")
library("BayesPen")
library("BCEE")
library("bacr")
library("hdm")
library("abe")
library("SignifReg")


OneData_Boruta <- function(mydata,var.list) {
  data.boruta <- mydata[,c("treat",var.list)]
  # fitA.boruta <- Boruta(treat~.,data=data.boruta)
  data.boruta <- cbind("Y"=mydata$Y, data.boruta)
  fitY.boruta <- Boruta(Y~.,data=data.boruta)
  
  # all covariates that are confirmed in outcome model
  select.boruta <- fitY.boruta$finalDecision[
    names(fitY.boruta$finalDecision) %in% var.list]=="Confirmed"
  
  if (sum(select.boruta)>0) {
    select.boruta <- names(fitY.boruta$finalDecision)[select.boruta]
    select.boruta <- grep(pattern="L[.]",select.boruta,value=TRUE)
    if (length(select.boruta)==0) {
      select.boruta <- "1" # intercept only
    }
  } else {
    select.boruta <- "1" # intercept only
  }
  return(select.boruta)
}

OneData_CovSel <- function(mydata,var.list,
                           cov.sel.type="dr",cov.sel.alg=1) {
  cov.sel.out <- tryCatch(
    out <- cov.sel(T = mydata$treat, Y = mydata$Y, X = mydata[,var.list], 
                   type = cov.sel.type, alg = cov.sel.alg, trace = 0),
    error=function(cond) return(NA))
  if (any(is.na(cov.sel.out))) {
    select.CovSel <- "1" # intercept only
  } else {
    if (cov.sel.alg==1) {
      select.CovSel <- var.list[(var.list %in% out$Q.0) |
                                  (var.list %in% out$Q.1)]
    } else if (cov.sel.alg==2) {
      select.CovSel <- var.list[(var.list %in% out$Z.0) |
                                  (var.list %in% out$Z.1)]
    }
    if (length(select.CovSel)==0) {
      select.CovSel <- "1" # intercept only
    }
  } 
  return(select.CovSel)
}

OneData_ABE <- function(mydata,var.list) {
  data.ABE <<- mydata[,c("Y","treat",var.list)]
  Y.binary <- all(0<=data.ABE$Y & data.ABE$Y<=1)
  if (Y.binary) {
    fitY.ABE <<- glm(Y~.,family=binomial(link = "logit"),
                     x=TRUE,y=TRUE,data=data.ABE)
  } else {
    fitY.ABE <<- lm(Y~.,x=TRUE,y=TRUE,data=data.ABE)
  }
  abe.fit <- abe(fit=fitY.ABE,data=data.ABE,include="treat",
                 exp.beta=Y.binary,active=var.list,verbose=FALSE)
  select.ABE <- names(abe.fit$coefficients)
  select.ABE <- grep(pattern="L[.]",select.ABE,value=TRUE)
  
  if (length(select.ABE)>0) {
    select.ABE <- paste0("L.",sort(as.integer(gsub("L[.]","",select.ABE))))
  } else {
    select.ABE <- "1" # intercept only
  }
  return(select.ABE)
}

OneData_SignifReg <- function(mydata,var.list) {
  data.SignifReg <<- mydata[,c("Y","treat",var.list)]
  
  # forward selection can only start with intercept-only model
  fitY.SignifReg <<- lm(Y~1,data=data.SignifReg)
  lm.SignifReg <- SignifReg(fitY.SignifReg, scope=as.formula(
    paste("~.+",paste(c("treat",var.list),collapse="+"))))
  select.SignifReg <- names(lm.SignifReg$coefficients)
  select.SignifReg <- grep(pattern="L[.]",select.SignifReg,value=TRUE)
  
  if (length(select.SignifReg)>0) {
    select.SignifReg <- paste0("L.",sort(
      as.integer(gsub("L[.]","",select.SignifReg))))
  } else {
    select.SignifReg <- "1" # intercept only
  }
  return(select.SignifReg)
}

OneData_BCEE <- function(mydata,var.list) {
  Y.binary <- all(0<=mydata$Y & mydata$Y<=1)
  # Using ABCEE to estimate the causal exposure effect
  if (Y.binary) {
    results <- ABCEE(X=mydata$treat, Y=mydata$Y, U=as.matrix(mydata[,var.list]),
                     family.X=binomial(link = "logit"),
                     omega = sqrt(nrow(mydata))*500)
  } else {
    results <- ABCEE(X=mydata$treat, Y=mydata$Y, U=as.matrix(mydata[,var.list]),
                     family.X=gaussian(link = "identity"),
                     omega = sqrt(nrow(mydata))*500)
  }
  # Tail prob. of 0 in the posterior dist. of the exposure effect:
  pvalue.bcee <- mean(results$betas>=0)
  #The posterior inclusion probability of each covariate:
  select.out.bcee <- colMeans(results$models.Y)
  # ordered covariates by decreasing posterior prob. of inclusion
  select.bcee <- var.list[order(select.out.bcee,decreasing=TRUE)]
  return(list(pvalue.bcee,select.bcee))
}

OneData_bacr <- function(mydata,var.list) {
  Y.binary <- all(0<=mydata$Y & mydata$Y<=1)
  bac.fam <- ifelse(Y.binary,"binomial","gaussian")
  
  ##### run BAC  #################
  result = bac(data=mydata, exposure="treat", outcome="Y", 
               confounders=var.list,
               interactors=NULL, familyX="binomial", familyY=bac.fam,
               num_its=5000,burnM=500,burnB=500,thin=10)
  # Tail prob. of 0 in the posterior dist. of the exposure effect:
  pvalue.bacr <- mean(result$ACE>=0)
  #The posterior inclusion probability of each covariate:
  select.out.bacr <- apply(result$models[[1]], 2, weighted.mean, 
                           w=result$models[[2]])
  # remove treatment
  select.out.bacr <- select.out.bacr[-length(select.out.bacr)]
  # ordered covariates by decreasing posterior prob. of inclusion
  select.bacr <- var.list[order(select.out.bacr,decreasing=TRUE)]
  return(list(pvalue.bacr,select.bacr))
}


OneData_ctmle <- function(mydata,var.list,pvalues=TRUE) {
  # C-TMLE LASSO for model selection of LASSO
  # from https://cran.r-project.org/web/packages/ctmle/vignettes/vignette.html
  
  # With initial estimate of Q
  Q <- cbind(rep(mean(mydata$Y[mydata$treat==0]), nrow(mydata)), 
             rep(mean(mydata$Y[mydata$treat==1]), nrow(mydata)))
  
  Y <- mydata$Y
  A <- mydata$treat
  W <- mydata[,var.list]
  
  ## discrete fit =============================================================
  ctmle_discrete_fit1 <- ctmleDiscrete(Y = Y, A = A, W = W, Q = Q,
                                       preOrder = FALSE, detailed = TRUE)
  
  glmnet_fit <- cv.glmnet(y = A, x = as.matrix(W),
                          family = 'binomial', nlambda = 20)
  
  ## algorithm 1 ==============================================================
  lambdas <- glmnet_fit$lambda[
    (which(glmnet_fit$lambda==glmnet_fit$lambda.min)):length(glmnet_fit$lambda)]
  
  Y.binary <- all(0<=Y & Y<=1)
  ctmleGlmnet.fam <- ifelse(Y.binary,"binomial","gaussian")
  
  ctmle_fit1 <- ctmleGlmnet(Y = Y, A = A, W = W, Q = Q,
                            lambdas = lambdas, ctmletype=1,
                            family=ctmleGlmnet.fam, gbound=0.025, V=5)

  ctmle_res <- rbind(
    c(ctmle_discrete_fit1$est,sqrt(ctmle_discrete_fit1$var.psi)),
    c(ctmle_fit1$est,sqrt(ctmle_fit1$var.psi))
    )
  rownames(ctmle_res) <- c("nolasso","lasso")
  colnames(ctmle_res) <- c("est","se")
  if (pvalues==TRUE) {
    ctmle_pv <- c(ctmle_discrete_fit1$pvalue,
                  ctmle_fit1$pvalue
                  )
    ctmle_res <- cbind("pv"=ctmle_pv,ctmle_res)
  }
  return(ctmle_res)
}

OneData_hdm <- function(mydata,var.list,pvalues=TRUE) {
  # Following example in \S4.4 of 
  # https://cran.r-project.org/web/packages/hdm/vignettes/hdm.pdf
  y = mydata$Y
  d = mydata$treat
  X = as.matrix(mydata[,var.list])
  
  lasso.effect = rlassoEffect(x = X, y = y, d = d, method = "partialling out")
  doublesel.effect = rlassoEffect(x = X, y = y, d = d, method = "double selection")
  rlasso.ate = rlassoATE(X,d,y)
  
  hdm_res <- rbind(
    c(lasso.effect$alpha,lasso.effect$se),
    as.numeric(c(doublesel.effect$alpha,doublesel.effect$se)),
    c(rlasso.ate$te,rlasso.ate$se))
  rownames(hdm_res) <- c("lasso","ds","rlasso")
  colnames(hdm_res) <- c("est","se")
  if (pvalues==TRUE) {
    hdm_pv <- c(lasso.effect$pval,
                as.numeric(doublesel.effect$pval),
                pnorm(q=(rlasso.ate$te/rlasso.ate$se)[1,1],lower.tail=FALSE))
    hdm_res <- cbind("pv"=hdm_pv,hdm_res)
  }
  return(hdm_res)
}
