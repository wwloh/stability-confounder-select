# function to calculate the IP weighted estimator (test statistic)
One_IPWest <- function(y_ipw,a_ipw,pl_ipw) {
  # weights using fitted propensity score
  ipw <- a_ipw/pl_ipw + (1-a_ipw)/(1-pl_ipw)
  return( sum(ipw*a_ipw*y_ipw)/sum(ipw*a_ipw)-
            sum(ipw*(1-a_ipw)*y_ipw)/sum(ipw*(1-a_ipw)) )
}

# function to calculate a p-value using a parametric bootstrap
One_pv_paraboot <- function(mydata,ps_fit,ps_true=NA,n_resample) {
  # ps_fit = propensity score model formula
  # ps_true = true PS
  if (any(is.na(ps_true))) {
    propscore.obs <- glm(ps_fit,family=binomial(link = "logit"),data=mydata)
    pA1_hat.obs <- predict.glm(propscore.obs,type="response")
  } else {
    pA1_hat.obs=ps_true
  }
  ts.obs <- One_IPWest(y_ipw=mydata$Y,a_ipw=mydata$treat,pl_ipw=pA1_hat.obs)
  ts <- sapply(1:n_resample, function(x) {
    xa <- rbinom(n,1,pA1_hat.obs)
    xa.dat <- mydata
    xa.dat$treat <- xa
    propscore.xa <- glm(ps_fit,family=binomial(link = "logit"),data=xa.dat)
    pA1_hat.xa <- predict.glm(propscore.xa,type="response")  
    if(!propscore.xa$converged) return(NA)
    ts.xa <- One_IPWest(y_ipw=mydata$Y,a_ipw=xa,pl_ipw=pA1_hat.xa)
    return(ts.xa)
  })
  pv <- mean(abs(ts) >= abs(ts.obs)-.Machine$double.eps*1e2, na.rm=TRUE)
  
  ## return min and max weights
  ipw.obs <- mydata$treat/pA1_hat.obs + (1-mydata$treat)/(1-pA1_hat.obs)
  ## return asymptotic p-value
  gee.obj <- geeglm(Y~treat, data=mydata, weights=ipw.obs, id=i,
                    corstr="independence")
  pv.gee <- summary(gee.obj)$coef["treat",4]
  
  res <- c(pv,pv.gee,range(ipw.obs))
  names(res) <- c("paraboot","gee","ipw.para.min","ipw.para.max")
  return(res)
}

# function to calculate a p-value given strata
One_pv_conditional <- function(mydata,n_resample) {
  strata <- sort(unique(mydata$stratum))
  strata_size <- as.integer(by(data=mydata$treat,INDICES=mydata$stratum,length))
  m_strata <- as.integer(by(data=mydata$treat,INDICES=mydata$stratum,sum))
  pA1_strata <- m_strata/strata_size
  names(pA1_strata) <- strata
  pA1_hat.strata <- pA1_strata[mydata$stratum]
  ts.obs <- One_IPWest(y_ipw=mydata$Y,a_ipw=mydata$treat,pl_ipw=pA1_hat.strata)
  ts <- sapply(1:n_resample, function(x) {
    xa <- rep(0L,n)
    xa1 <- sapply(strata, function(s) {
      s.idx <- which(mydata$stratum==s)
      sample(x=s.idx,size=m_strata[s],replace=FALSE)
    })
    xa[unlist(xa1)] <- 1L # table(xa,subdata$subclass)
    ts.xa <- One_IPWest(y_ipw=mydata$Y,a_ipw=xa,pl_ipw=pA1_hat.strata)
    return(ts.xa)
  })
  pv <- mean(abs(ts) >= abs(ts.obs)-.Machine$double.eps*1e2, na.rm=TRUE)
  
  res <- pv
  names(res) <- "conditional"
  return(res)
}

# function to calculate a p-value using full matching
One_pv_strata <- function(mydata,ps_fit,...) {
  if (class(ps_fit)=="formula") {
    # form strata based on full matching --------------------
    m.out <- matchit(ps_fit, data = mydata, method = "full")
  } else {
    # form strata based on full matching of provided PS
    mydata[,"pA1_hat"] <- ps_fit
    m.out <- matchit(treat~pA1_hat, data = mydata, method = "full")
  }
  # sum(m.out$discarded)
  # table(m.out$subclass,m.out$treat)
  # plot(predict.glm(glm(ps.fit,data=mydata),type="response")~m.out$subclass)
  subdata <- cbind(mydata,"subclass"=m.out$subclass)
  
  # conditional inference within each subclass ------------
  subdata_cond <- subdata
  names(subdata_cond)[which(names(subdata_cond)=="subclass")] <- "stratum"
  pv_cond <- One_pv_conditional(mydata=subdata_cond,...)
  return(pv_cond)
}

ForwardSelect_minimax <- function(Y,X,Z,criterion="minimin") {
  X.names <- colnames(X)
  p <- length(X.names)
  n <- length(Y)
  X.curr <- NULL
  X.crit <- NULL
  while(length(X.curr) < min(n,p)) {
    # find covariate among those not already in the model with minimax p-value
    X.cands <- X.names[!(X.names %in% X.curr)]
    Xj.crit <- sapply(X.cands, function(X.cand) {
      # create dataset with current and candidate covariates
      data.X.cand <- cbind("i"=1:n,"Y"=Y,"treat"=Z,
                           X[,c(X.curr,X.cand),drop=FALSE])
      # propensity score and/or outcomes models with the included covariates
      if (criterion!="Outcome_only") {
        ps.X.cand <- as.formula(paste0("treat~",paste(
          c(X.curr,X.cand),collapse="+")))
        ps.X.cand.lm <- lm(ps.X.cand,data=data.X.cand)
        ps.X.cand.pv <- coef(summary(ps.X.cand.lm))[X.cand,"Pr(>|t|)"]
      }
      if (criterion!="PS_only") {
        out.X.cand <- as.formula(paste0("Y~",paste(
          c("treat",X.curr,X.cand),collapse="+")))
        out.X.cand.lm <- lm(out.X.cand,data=data.X.cand)
        out.X.cand.pv <- coef(summary(out.X.cand.lm))[X.cand,"Pr(>|t|)"]
      }
      # evaluate criterion
      if (criterion=="minimax") {
        max(ps.X.cand.pv,out.X.cand.pv)
      } else if (criterion=="minimin") {
        min(ps.X.cand.pv,out.X.cand.pv)
      } else if (criterion=="Outcome_only") {
        out.X.cand.pv
      } else if (criterion=="PS_only") {
        ps.X.cand.pv
      }
    })
    X.add <- X.cands[which.min(Xj.crit)]
    X.curr <- c(X.curr,X.add)
    X.crit <- c(X.crit,min(Xj.crit))
  }
  return(list(X.curr,X.crit))
}

OrderCovariateIndices <- function(L.ordered) {
  # L.ordered: vector of ordered covariates of form L.j, where j = integer
  Lj.ordered <- as.integer(unlist(
    lapply(strsplit(x=L.ordered,split="L."),"[[",2)))
  data.frame("L"=Lj.ordered[order(Lj.ordered)],
             "order"=order(Lj.ordered))
}

TreatmentHeterogeneity <- function(x,s,k=NA) {
  if (is.na(k)) k <- length(x)
  # indices for subsets of size k
  windows <- sapply(0:(length(x)-k), function(i) (1:k)+i, simplify="matrix")
  # calculate Cochran's Q for each subset
  window.summ <- apply(windows, 2, function(i) {
    xi <- x[i]
    wi <- s[i]^(-2) # weights
    wi[is.na(wi)] <- 0 # zero weights if zero std err
    xbarw <- sum(wi*xi)/sum(wi)
    CochranQ <- sum(wi*((xi-xbarw)^2))
    I2 <- max(0,1-(length(i)-1)/CochranQ)
    return(c("Q"=CochranQ,"I2"=I2))
  })
  # index for center of each window
  colnames(window.summ) <- windows[ceiling(k/2),]
  return(window.summ)
}

DiagnosticsCovariatesOrdered <- function(L.ordered,mydata,M) {
  # L.ordered: sequence of ordered covariate indices
  diags.list <- lapply(1:length(L.ordered), function(x){
    Ls <- L.ordered[1:x]
    # PS model (linear)
    ps.Ls <- as.formula(paste0("treat~",paste(Ls,collapse="+")))
    ps.lm <- lm(ps.Ls,data=mydata)
    # Outcome model (linear)
    out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
    out.lm <- lm(out.Ls,data=mydata)
    # individual contributions to score equations
    score.Ls <- (mydata$treat-predict.lm(ps.lm))*(mydata$Y-predict.lm(out.lm))/
      sum((mydata$treat-predict.lm(ps.lm))^2)
    list(
      c(One_pv_strata(mydata,ps_fit=ps.Ls,n_resample=M),
        coef(out.lm)["treat"],
        "treat.se.lm"=summary(out.lm)$coef["treat","Std. Error"],
        "treat.se.score"=sqrt(sum(score.Ls^2))),
      score.Ls)
  })
  
  diags <- do.call(rbind,lapply(diags.list, "[[", 1))
  score.list <- do.call(rbind,lapply(diags.list, "[[", 2))
  
  # full PS model
  ATE.hat.allL <- diags[nrow(diags),"treat"]
  score.allL <- score.list[nrow(score.list),]
  
  diags <- data.frame(cbind(
    diags,
    "score.SE.vs.allL"=apply(score.list, 1, function(x) {
      sqrt( sum((x^2)+(score.allL^2)-(2*x*score.allL)) )})
  ))
  diags <- cbind(
    diags,
    # standardized diff in ttmt effect ests
    "score.vs.allL"=(diags$treat-ATE.hat.allL)/diags$score.SE.vs.allL
  )
  diags[nrow(diags),c("score.vs.allL","score.SE.vs.allL")] <- 
    c(0,NA) # set NaN score to zero and zero SE to NA
  
  # Q statistics for window width 3
  QI2.allL <- TreatmentHeterogeneity(
    x=diags$treat-ATE.hat.allL,s=diags$score.SE.vs.allL,k=3)
  crits <- c("score.vs.allL.Q.3"=as.integer(names(which.min(QI2.allL["Q",]))))
  return(list(diags,crits))
}

ATEhat_lm <- function(mydata, Ls, point.only=TRUE) {
  # Outcome model
  out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
  out.lm <- lm(out.Ls,data=mydata)
  if (point.only==FALSE) {
    # PS model (linear)
    ps.Ls <- as.formula(paste0("treat~",paste(Ls,collapse="+")))
    ps.lm <- lm(ps.Ls,data=mydata)
    # individual contributions to score equations
    score.Ls <- (mydata$treat-predict.lm(ps.lm))*(mydata$Y-predict.lm(out.lm))/
      sum((mydata$treat-predict.lm(ps.lm))^2)
    return(c(coef(out.lm)["treat"],
             "treat.se.lm"=summary(out.lm)$coef["treat","Std. Error"],
             "treat.se.score"=sqrt(sum(score.Ls^2))))
  } else {
    return(as.numeric(coef(out.lm)["treat"]))
  }
}


ATEhat_gee <- function(mydata,ps_fit,ps_true=NA) {
  # ps_fit = propensity score model formula
  # ps_true = true PS
  if (any(is.na(ps_true))) {
    propscore.obs <- glm(ps_fit,family=binomial(link = "logit"),data=mydata)
    pA1_hat.obs <- predict.glm(propscore.obs,type="response")
  } else {
    pA1_hat.obs=ps_true
  }
  
  ## IP weights
  ipw.obs <- mydata$treat/pA1_hat.obs + (1-mydata$treat)/(1-pA1_hat.obs)
  ## return point and sandwich variance estimates
  gee.obj <- geeglm(Y~treat, data=mydata, weights=ipw.obs, id=i,
                    corstr="independence")
  res <- summary(gee.obj)$coef["treat",1:2]
  names(res) <- paste0("gee.",c("est","se"))
  return(res)
}
