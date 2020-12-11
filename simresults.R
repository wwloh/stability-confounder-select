rm(list=ls())
library("data.table")

# results ---------------------------------------------------------------------
subfolder <- "sims-v7/"
myfiles <- list.files(subfolder)
myfiles <- grep(pattern=".Rdata",myfiles,value=TRUE)
res_list <- list()
for (ll in myfiles) {
  # results
  load(paste0(subfolder,ll))
  sim_res.dt <- data.table(do.call(rbind,lapply(sim_res, unlist)))
  setkey(sim_res.dt)
  res_list <- c(res_list,list(sim_res.dt))
}
sim_res_all <- rbindlist(res_list)
setkey(sim_res_all)
# number of sims for each setting
(simset <- sim_res_all[,.N,by=list(n,p,z,q,y)])
simset[, N := NULL]

# helper function for adding an ECDF step function to a plot
my_ecdf <- function(x, ...) {
  Fn <- ecdf(x)
  xx <- sort(unique(x))
  lines(xx, Fn(xx), col.01line=NULL, ...)
}

for (ll in c(1)) {
  sim_res <- NULL
  sim_res <- sim_res_all[sim_res_all$n==simset[ll,n] & 
                           sim_res_all$p==simset[ll,p] &
                           sim_res_all$z==simset[ll,z] &
                           sim_res_all$q==simset[ll,q] & 
                           sim_res_all$y==simset[ll,y],]
  filename <- paste(names(simset),simset[ll],sep="_",collapse="-")
  filename <- gsub(pattern="[.]",replacement="_",x=filename)
  
  # selected covariates using different methods
  ## prob. of selecting 0,>=1,2 confounders vs. 
  ## av./sd of number of selected covariates
  sel_true <- sim_res[,grep("selected.",names(sim_res_all),value=TRUE),with=FALSE]
  sel_size <- sim_res[,grep("PSsizes",names(sim_res_all),value=TRUE),with=FALSE]
  meth_names <- unlist(lapply(strsplit(names(sel_true),"selected."),"[",2))
  sel.summ <- NULL
  for (mm in 1:length(meth_names)) {
    sel.summ[[mm]] <- c(
      "select_geq1"=mean(unlist(sel_true[,..mm])>0,na.rm=TRUE),
      "select_both"=mean(unlist(sel_true[,..mm])==2,na.rm=TRUE),
      "size_av"=mean(unlist(sel_size[,..mm]),na.rm=TRUE),
      "size_lq"=quantile(unlist(sel_size[,..mm]),probs=0.25,na.rm=TRUE),
      "size_uq"=quantile(unlist(sel_size[,..mm]),probs=0.75,na.rm=TRUE))
  }
  sel.summ <- data.frame(do.call(rbind,sel.summ))
  meth_names <- gsub("ds_stable","Stability",meth_names)
  sel.summ <- cbind("meth"=meth_names,sel.summ)
  
  for (ss in 1:2) {
    pdf(paste0("plots-sims-selectVsPSsize-",ss,"-",filename,".pdf"),
        height=5,width=5)
    plot(0,0,ylim=c(0,1),xlim=c(0,simset[ll,p]),type="n",
         xlab="Average number of covariates selected",
         ylab=ifelse(ss==1, 
                     "Prob. of selecting >= 1 confounder",
                     "Prob. of selecting all confounders"))
    for (mm in 1:nrow(sel.summ)) {
      points(x=sel.summ[mm,"size_av"],
             y=sel.summ[mm,ss+1],
             pch=4,cex=.7)
      lines(x=sel.summ[mm,c("size_lq.25.","size_uq.75.")],
            y=rep(sel.summ[mm,ss+1],2),
            type="o",pch=1,cex=.5,
            col="grey50")
      text.pos.x <- sel.summ[mm,"size_av"]
      if (sel.summ[mm,"meth"]=="SignifReg" | sel.summ[mm,"meth"]=="CovSel") {
        text.pos <- 3
      } else if (sel.summ[mm,"meth"]=="Boruta" | sel.summ[mm,"meth"]=="Stability") {
        text.pos <- 1
      } else {
        text.pos <- 4
        text.pos.x <- sel.summ[mm,"size_uq.75."]
      }
      text(x=text.pos.x,
           y=sel.summ[mm,ss+1],
           pos=text.pos,
           label=sel.summ[mm,"meth"],cex=.7)
    }
    abline(a=0,b=1/simset[ll,p],lty=2)
    dev.off()
  }

  # ECDFs of p-values =========================================================
  pdf(paste0("plots-sims-typeI-",filename,".pdf"),
      height=5,width=10)
  par(mfrow=c(1,2))
  plot.ecdf(sim_res$conditional.target, 
            do.points=FALSE,col.01line=NULL,
            xlim=c(0,1),ylim=c(0,1), 
            xlab="p-value",ylab="Pr(X<=x)",
            main="Known PS model")
  abline(a=0,b=1,col="grey50")
  my_ecdf(sim_res$conditional.noL, lty=2)
  legend("bottomright",lty=c(1:3),bty="n",
         legend=c("Target PS model",
                  "Empty PS model"
                  ))
  plot.ecdf(sim_res$conditional.ds_stable, 
            do.points=FALSE,col.01line=NULL,
            xlim=c(0,1),ylim=c(0,1), 
            xlab="p-value",ylab="Pr(X<=x)",
            main="Selected PS model")
  abline(a=0,b=1,col="grey50")
  my_ecdf(sim_res$conditional.ABE, lty=2)
  my_ecdf(sim_res$conditional.CovSel, lty=3)
  legend("bottomright",lty=c(1:5),bty="n",
         legend=c("Stability",
                  "ABE",
                  "CovSel"))
  dev.off()
  
  # other covariate selection methods
  pdf(paste0("plots-sims-typeI-others1-",filename,".pdf"),
      height=5,width=10)
  par(mfrow=c(1,2))
  plot.ecdf(sim_res$conditional.Boruta, 
            do.points=FALSE,col.01line=NULL,
            xlim=c(0,1),ylim=c(0,1), 
            xlab="p-value",ylab="Pr(X<=x)",
            main="Full Matching")
  abline(a=0,b=1,col="grey50")
  my_ecdf(sim_res$conditional.SignifReg, lty=2)
  legend("bottomright",lty=c(1:2),bty="n",
         legend=c("Boruta",
                  "SignifReg"))
  
  plot.ecdf(sim_res$gee.OAL, 
            do.points=FALSE,col.01line=NULL,
            xlim=c(0,1),ylim=c(0,1), 
            xlab="p-value",ylab="Pr(X<=x)",
            main="Wald tests")
  abline(a=0,b=1,col="grey50")
  my_ecdf(sim_res$hdm.ds, lty=2)
  if (!all(is.na(sim_res$ctmle.nolasso)) & 
      !all(is.na(sim_res$ctmle.lasso))) {
    my_ecdf(sim_res$ctmle.nolasso, lty=3)
    my_ecdf(sim_res$ctmle.lasso, lty=4, lwd=1.2)
    legend("bottomright",lty=c(1:4),bty="n",
           legend=c("OAL",
                    "HDM",
                    "CTMLE: discrete",
                    "CTMLE: LASSO"))
  } else {
    legend("bottomright",lty=c(1,4),bty="n",
           legend=c("HDM","OAL"))
  }
  dev.off()
  
  if (!all(is.na(sim_res$posterior.bacr)) & 
      !all(is.na(sim_res$posterior.BCEE))) {
    pdf(paste0("plots-sims-typeI-others2-",filename,".pdf"),
        height=5,width=5)
    par(mfrow=c(1,1))
    plot.ecdf(sim_res$posterior.bacr, 
              do.points=FALSE,col.01line=NULL,
              xlim=c(0,1),ylim=c(0,1), 
              xlab="p-value",ylab="Pr(X<=x)",
              main="Bayesian methods")
    abline(a=0,b=1,col="grey50")
    my_ecdf(sim_res$posterior.BCEE, lty=2)
    legend("bottomright",lty=1:2,bty="n",
           legend=c("BACR","BCEE"))
    dev.off()
  }
}

mymerge = function(x,y) merge(x,y,all=TRUE)

# treatment effect estimates post-selection
sel_ate <- sim_res_all[
  ,.SD,
  .SDcols=grep("treat.",names(sim_res_all),value=TRUE),
  by=list(n,q,p,z,y)]
meth_names <- unlist(lapply(strsplit(names(sel_ate),"treat."),"[",2))
meth_names <- unlist(strsplit(grep(pattern="[.]se",meth_names,value=TRUE),"[.]se"))
# meth_names <- meth_names[grepl("hdm.ds",meth_names) | grepl("ctmle",meth_names)]
meth_names
sel.summ <- NULL
for (mm in 1:length(meth_names)) {
  sim_est_meth <- sel_ate[
    ,.SD,
    .SDcols=grep(meth_names[mm],names(sel_ate),value=TRUE),
    by=list(n,q,p,z,y)]
  setnames(sim_est_meth,old=ncol(sim_est_meth)-(1:0),new=c("est","se"))
  setkey(sim_est_meth)
  
  sim_est_meth <- sim_est_meth[,list("mean.pt"=mean(est,na.rm=TRUE),
                                     "ese"=sd(est,na.rm=TRUE),
                                     "mean.se"=mean(se,na.rm=TRUE),
                                     "ase"=sd(se,na.rm=TRUE),
                                     "mse.pt"=mean(sqrt(est^2+se^2),na.rm=TRUE))
                               ,by=list(n,q,p,z,y)]
  
  # NAs
  sim_est_meth.na <- sim_res_all[,lapply(.SD,function(x) mean(is.na(x))),
                                 .SDcols=grep(meth_names[mm],names(sel_ate),value=TRUE)[1],
                                 by=list(n,q,p,z,y)]
  setnames(sim_est_meth.na,ncol(sim_est_meth.na),"nas")
  setkey(sim_est_meth.na)
  
  # p-values
  sim_est_meth.pv <- sim_res_all[,lapply(.SD,function(x) mean(x<=0.05,na.rm=TRUE)),
                                 .SDcols=grep(meth_names[mm],names(sim_res_all),value=TRUE)[1],
                                 by=list(n,q,p,z,y)]
  setnames(sim_est_meth.pv,ncol(sim_est_meth.pv),"pv")
  setkey(sim_est_meth.pv)
  
  sel.summ[[mm]] <- cbind("meth"=meth_names[mm],
                          Reduce(mymerge,list(
                            sim_est_meth.na,sim_est_meth.pv,sim_est_meth)))
}
sel.summ <- do.call(rbind,sel.summ)
setcolorder(sel.summ,c(2:6,1,7:ncol(sel.summ)))
setkey(sel.summ)

round(dcast.data.table(sel.summ, formula=q+p+z+y ~ meth, value.var = "nas"),2)

library("xtable")
xtable(sel.summ[q==0 & p==25 & z==2 & y==0])

print(xtable(dcast.data.table(
  sel.summ, 
  formula=q+p+z+y ~ meth, 
  value.var = c("pv","mean.pt","ese","mean.se","ase"))
  ),include.rownames=FALSE)

# covariate selection methods
meth_names <-  c("noL","target","ds_stable","Boruta","CovSel","ABE","SignifReg")
sel.summ <- NULL
for (mm in 1:length(meth_names)) {
  sim_est_meth <- sim_res_all[
    ,.SD,
    .SDcols=grep(meth_names[mm],names(sim_res_all),value=TRUE),
    by=list(n,q,p,z,y)]
  if (ncol(sim_est_meth)==6) {
    sim_est_meth <- cbind(sim_est_meth,NA,NA,NA,NA)
  }
  setnames(sim_est_meth,old=6:ncol(sim_est_meth),
           new=c("pv","pssize","lsel","est","se"))  
  
  sim_est_meth <- sim_est_meth[,list(
    "l_selboth"=mean(lsel==2,na.rm=FALSE),
    "l_selsize.mean"=mean(pssize,na.rm=FALSE),
    "typeI_5"=mean(pv<=0.05,na.rm=TRUE))
    ,by=list(n,q,p,z,y)]
  sel.summ[[mm]] <- cbind("meth"=meth_names[mm],sim_est_meth)
}
sel.summ <- do.call(rbind,sel.summ)
setcolorder(sel.summ,c(2:6,1,7:ncol(sel.summ)))
setkey(sel.summ)
sel.summ[meth=="noL", meth := "AA_noL"]
sel.summ[meth=="target", meth := "z_target"]
setkey(sel.summ)

sel.summ.wide <- dcast.data.table(
  sel.summ, formula=q+p+z+y ~ meth, 
  value.var = c("l_selboth","l_selsize.mean","typeI_5"))
sel.summ.wide[,c("l_selboth_AA_noL","l_selboth_z_target",
                 "l_selsize.mean_AA_noL","l_selsize.mean_z_target") := NULL]

print(xtable(sel.summ.wide),include.rownames=FALSE)
