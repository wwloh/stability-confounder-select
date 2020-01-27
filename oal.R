OneData_OAL <- function(var.list,Data) {
  
# set vector of possible lambda's to try
lambda_vec = c( -10, -5, -2, -1, -0.75, -0.5, -0.25, 0.25, 0.49)
names(lambda_vec) = as.character(lambda_vec)
# lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
gamma_convergence_factor = 2
# get the gamma value for each value in the lambda vector that corresponds to convergence factor
gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
names(gamma_vals) = names(lambda_vec)

### define some functions for generating data, ATE estimates, and the wAMD,
expit = function(x){ 
 pr = ( exp(x) / (1+exp(x)) ) 
 return(pr)
}
ATE_est = function(fY,fw,fA){
 t_ATE = fY*fw
 tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
 return(tt_ATE) 
}
create_weights = function(fp,fA,fw){
 fw = (fp)^(-1)
 fw[fA==0] = (1 - fp[fA==0])^(-1)
 return(fw)
}
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
 trt = untrt = diff_vec = rep(NA,length(beta)) 
 names(trt) = names(untrt) = names(diff_vec) = varlist
 for(jj in 1:length(varlist)){ 
  this.var = paste("w",varlist[jj],sep="") 
  DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
  trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
  untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
  diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
 } 
 wdiff_vec = diff_vec * abs(beta) 
 wAMD = c( sum(wdiff_vec))
 ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
 return(ret) 
}

# estimate outcome model
y.form = formula(paste("Y~A+",paste(var.list,collapse="+")))
lm.Y = lm(y.form,data=Data)
betaXY = coef(lm.Y)[var.list] 

## Want to save ATE, wAMD and propensity score coefficients for each lambda value
ATE = wAMD_vec = rep(NA, length(lambda_vec))
names(ATE) = names(wAMD_vec) = names(lambda_vec)
coeff_XA = as.data.frame(matrix(NA,nrow=length(var.list),ncol=length(lambda_vec)))
names(coeff_XA) = names(lambda_vec)
rownames(coeff_XA) = var.list

######################################################################################
#####  Run outcome adaptive lasso for each lambda value 
######################################################################################
# weight model with all possible covariates included, this is passed into lasso function
w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
for( lil in names(lambda_vec) ){
  il = lambda_vec[lil]
  ig = gamma_vals[lil]

  ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
  oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
  ### run outcome-adaptive lasso model with appropriate penalty
  logit_oal = lqa( w.full.form, data=Data, penalty=oal_pen, family=binomial(logit) )
  # generate propensity score
  Data[,paste("f.pA",lil,sep="")] = predict(logit_oal)$mu.new
  # save propensity score coefficients
  coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
  # create inverse probability of treatment weights
  Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
  # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
  wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
			wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
  # save ATE estimate for this lambda value
  ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
} # close loop through lambda values

# print out wAMD for all the lambda values tried
wAMD_vec
# find the lambda value that creates the smallest wAMD
tt = which.min(wAMD_vec)
# print out ATE corresponding to smallest wAMD value
ATE[tt]
# print out the coefficients for the propensity score that corresponds with smalles wAMD value 
coeff_XA[,tt]

return(list("Data"=Data,"coeff_XA"=coeff_XA,"tt"=tt, "ATE"=ATE[tt]))
}

# function to calculate a p-value using a parametric bootstrap with OAL
One_pv_paraboot_OAL <- function(mydata,ps_fit,n_resample) {
  # ps_fit = estimated propensity scores from OAL
  pA1_hat.obs <- ps_fit
  ts.obs <- One_IPWest(y_ipw=mydata$Y,a_ipw=mydata$treat,pl_ipw=pA1_hat.obs)
  ts <- sapply(1:n_resample, function(x) {
    xa <- rbinom(n,1,pA1_hat.obs)
    xa.dat <- mydata
    xa.dat$A <- xa
    xa.var.list <- names(mydata)[grep("L[.]",names(mydata))]
    xa.res_oal <- OneData_OAL(xa.var.list,xa.dat)
    pA1_hat.xa <- xa.res_oal$Data[,paste("f.pA",names(xa.res_oal$tt),sep="")]
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