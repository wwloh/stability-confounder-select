# function to calculate the test statistic for randomization inference
One_TS <- function(y_ts,a_ts,x_ts) {
  return( sum(by(data=y_ts*a_ts, INDICES=x_ts, FUN=function(ya) 
    sum(ya)*length(ya)))/length(y_ts) )
}

# function to calculate a p-value for a Wald test using GEE variance
One_IPW_GEE <- function(mydata,ps_fit) {
  ipw.obs <- mydata$treat/ps_fit + (1-mydata$treat)/(1-ps_fit)
  ## return asymptotic p-value
  gee.obj <- geeglm(Y~treat, data=mydata, weights=ipw.obs, id=i,
                    corstr="independence")
  pv.gee <- summary(gee.obj)$coef["treat","Pr(>|W|)"]
  
  res <- pv.gee
  names(res) <- "gee"
  return(res)
}

# function to calculate a p-value given strata
One_pv_conditional <- function(mydata,n_resample) {
  ts.obs <- One_TS(y_ts=mydata$Y,a_ts=mydata$treat,x_ts=mydata$stratum)
  strata <- sort(unique(mydata$stratum))
  m_strata <- as.integer(by(data=mydata$treat,INDICES=mydata$stratum,sum))
  ts <- sapply(1:n_resample, function(x) {
    xa <- rep(0L,n)
    xa1 <- sapply(strata, function(s) {
      s.idx <- which(mydata$stratum==s)
      sample(x=s.idx,size=m_strata[s],replace=FALSE)
    })
    xa[unlist(xa1)] <- 1L
    # check same number assigned to treatment within each stratum
    ts.xa <- NA
    if (all(m_strata==as.integer(by(data=xa,INDICES=mydata$stratum,sum)))) {
      ts.xa <- One_TS(y_ts=mydata$Y,a_ts=xa,x_ts=mydata$stratum)  
    }
    return(ts.xa)
  })
  pv <- mean(abs(ts) >= (abs(ts.obs)-.Machine$double.eps*1e2), na.rm=FALSE)
  
  res <- pv
  names(res) <- "conditional"
  return(res)
}

# function to calculate a p-value using full matching
One_pv_strata <- function(mydata,ps_fit,...) {
  # form strata based on full matching --------------------
  m.out <- matchit(ps_fit, data = mydata, method = "full")
  # table(m.out$subclass,m.out$treat); length(unique(m.out$subclass))
  # table(rowSums(table(m.out$subclass,m.out$treat)))
  subdata <- cbind(mydata,"stratum"=m.out$subclass)
  
  # conditional inference within each subclass ------------
  pv_cond <- One_pv_conditional(mydata=subdata,...)
  return(pv_cond)
}

ForwardSelect_DS <- function(Y,X,A) {
  Y.binary <- all(0<=Y & Y<=1)
  X.names <- colnames(X)
  p <- length(X.names)
  n <- length(Y)
  X.curr <- NULL
  X.crit <- NULL
  while(length(X.curr) < p) {
    # find covariate among those not already in the model with minimin p-value
    X.cands <- X.names[!(X.names %in% X.curr)]
    Xj.crit <- sapply(X.cands, function(X.cand) {
      # create dataset with current and candidate covariates
      data.X.cand <- cbind("i"=1:n,"Y"=Y,"treat"=A,
                           X[,c(X.curr,X.cand),drop=FALSE])
      
      # treatment and outcomes models with the included covariates
      ps.X.cand <- as.formula(paste0("treat~",paste(
        c(X.curr,X.cand),collapse="+")))
      ps.X.cand.glm <- glm(ps.X.cand,family=binomial(link = "logit"),
                           data=data.X.cand)
      ps.X.cand.pv <- coef(summary(ps.X.cand.glm))[X.cand,"Pr(>|z|)"]
      
      out.X.cand <- as.formula(paste0("Y~",paste(
        c("treat",X.curr,X.cand),collapse="+")))
      
      if (Y.binary) {
        out.X.cand.lm <- glm(out.X.cand,family=binomial(link = "logit"),
                             data=data.X.cand)
        out.X.cand.pv <- coef(summary(out.X.cand.lm))[X.cand,"Pr(>|z|)"]
      } else {
        out.X.cand.lm <- lm(out.X.cand,data=data.X.cand)
        out.X.cand.pv <- coef(summary(out.X.cand.lm))[X.cand,"Pr(>|t|)"]
      }
      
      # evaluate criterion
      ps.OK <- all( !is.na(coef(summary(ps.X.cand.glm))[-1,"Pr(>|z|)"]) &
                      (coef(summary(ps.X.cand.glm))[-1,"Pr(>|z|)"]>0) )
      if (!ps.OK) {
        ps.X.cand.pv <- Inf
      }
      ## check for perfect separation in outcome model (if Y binary)
      out.OK <- TRUE
      if (Y.binary) {
        out.OK <- all( !is.na(coef(summary(out.X.cand.lm))[-1,"Pr(>|z|)"]) &
                         (coef(summary(out.X.cand.lm))[-1,"Pr(>|z|)"]>0) )
      }
      if (!out.OK) {
        out.X.cand.pv <- Inf
      }
      return( min(ps.X.cand.pv,out.X.cand.pv,na.rm=TRUE) )
    })
    if (min(Xj.crit)<=1) {
      X.add <- X.cands[which.min(Xj.crit)] # valid p-values
    } else {
      X.add <- sample(x=X.cands,size=1) # randomly select a covariate
    }
    X.curr <- c(X.curr,X.add)
    X.crit <- c(X.crit,min(Xj.crit))
  }
  return(list("ordered"=X.curr,"pvalues"=X.crit))
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
    wi[is.infinite(wi)] <- 0 # zero weights if zero std err
    xbarw <- weighted.mean(x=xi,w=wi)
    CochranQ <- sum(wi*((xi-xbarw)^2))
    return( CochranQ )
  })
  # index for center of each window
  names(window.summ) <- windows[median(1:k),]
  return(window.summ)
}

OneDRstd_Est <- function(Ls,mydata,varest=FALSE) {
  # PS model
  ps.Ls <- as.formula(paste0("treat~",paste(Ls,collapse="+")))
  ps.glm <- glm(ps.Ls,family=binomial(link = "logit"),data=mydata)
  ps.glm.hat <- predict.glm(ps.glm,type="response")
  ipw.ps.glm <- mydata$treat/ps.glm.hat + (1-mydata$treat)/(1-ps.glm.hat)
  rm(ps.Ls,ps.glm,ps.glm.hat)
  
  # trim observations with extreme weights from sample
  ipw.ps.glm[ipw.ps.glm > 1/(.Machine$double.eps*1e1)] <- 0
  
  # Outcome model fitted via weighted regression
  out.Ls <- as.formula(paste0("Y~",paste(c("treat",Ls),collapse="+")))
  if (all(0<=mydata$Y & mydata$Y<=1)) {
    # logistic regression model
    out.lm <- glm(out.Ls,family=quasibinomial(link = "logit"),data=mydata,
                  weights=ipw.ps.glm)
  } else {
    # linear regression model
    out.lm <- lm(out.Ls,data=mydata,weights=ipw.ps.glm)
  }
  # for predicting counterfactuals
  mydata.A1 <- mydata.A0 <- mydata
  mydata.A1[,"treat"] <- 1
  mydata.A0[,"treat"] <- 0
  infn.est <- (2*mydata$treat-1)*ipw.ps.glm*
    (mydata$Y - predict(out.lm,type="response")) + # residuals
    predict(out.lm,newdata=mydata.A1,type="response") -
    predict(out.lm,newdata=mydata.A0,type="response")
  # one-step plug-in estimator
  ate.hat <- mean(infn.est)
  # individual influence functions
  infn.est <- infn.est - ate.hat
  
  if (varest==TRUE) {
    return( c("ate"=ate.hat,"se"=sqrt(var(infn.est)/length(infn.est))) ) 
  } else {
    return( list("ate"=ate.hat,"infn"=infn.est) )  
  }
}

DiagnosticsCovariatesOrdered <- function(L.ordered,mydata,k) {
  # L.ordered: sequence of ordered covariate indices
  diags.list <- lapply(1:length(L.ordered), function(x){
    OneDRstd_Est(Ls=L.ordered[1:x],mydata=mydata)
  })
  
  diags <- do.call(rbind,lapply(diags.list, "[[", "ate"))
  score.list <- do.call(rbind,lapply(diags.list, "[[", "infn"))
  # all(abs(rowMeans(score.list)) < .Machine$double.eps) # should be almost zero
  
  diags <- cbind(diags, 
                 # variance est
                 "var"=apply(score.list, 1, var)/ncol(score.list),
                 # diff in ttmt effect ests
                 "dif"=diags[,1]-diags[nrow(diags),1])
  # var est of diff
  var.diff <- sapply(1:(nrow(score.list)-1), function(j) {
    var(score.list[j,]-score.list[nrow(score.list),])/ncol(score.list)
  })
  diags <- data.frame(cbind(diags,"var.dif"=c(var.diff,0)))
  colnames(diags)[1] <- "ate"
  
  # Q statistics for given window width
  Qk <- TreatmentHeterogeneity(x=diags$dif,s=sqrt(diags$var.dif),k=k)
  crits <- as.integer(names(which.min(Qk)))
  names(crits) <- paste0("score.vs.allL.Q.",k)
  return(list("est"=diags,"selected_orbit"=crits,Qk))
}
