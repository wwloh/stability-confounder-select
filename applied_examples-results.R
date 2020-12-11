rm(list=ls())
libraries_check <- c("data.table","MatchIt","optmatch","geepack","MASS","lqa",
                     "Hmisc","CovSel","speff2trial","locfit","xtable")
# install.packages("lqa_1.0-3.tar.gz", repos = NULL, type="source")
for (libs in libraries_check) {
  if(!libs %in% rownames(installed.packages())) {
    install.packages(libs,repos="http://lib.ugent.be/CRAN/")
  }
  library(libs, character.only=TRUE)
}
rm(libraries_check, libs)

# results =====================================================================
dataset.names <- c("rhc5","rhc7","rhc9","speff2","lalonde")
for (data.name in dataset.names) {
  if (grepl("rhc",data.name)) {
    source("data-prep-rhc.R")
  } else {
    source(paste0("data-prep-",data.name,".R"))
  }
  res_filename <- paste0("stability-applied_examples-",data.name,".Rdata")
  load(file=res_filename)
  
  # ordering based on minimin
  Lnames.ordered <- Lnames[na.exclude(as.integer(unlist(
    strsplit(res.FwdSel$ds_stable$ordered,split="L.")))),]
  L.selected <- lapply(L.selected, function(covselected) {
    if (any(covselected=="1")) {
      l.selected <- NULL
    } else {
      l.selected <- na.exclude(as.integer(unlist(strsplit(covselected,"L."))))
    }
    cbind("order"=1:nrow(Lnames.ordered),
          Lnames.ordered,
          "selected"=Lnames.ordered[,"L.idx"] %in% l.selected)
  })
  L.selected.dt <- data.table(L.selected$ds_stable)
  for (ll in 2:length(L.selected)) {
    if (!any(is.na(L.selected[[ll]]))) {
      ll.dt <- as.data.table(L.selected[[ll]])
      setnames(ll.dt,"selected",names(L.selected)[ll])
      setkey(ll.dt)
      L.selected.dt <- merge(L.selected.dt,ll.dt,
                             by=c("order","L.idx","L.names"), all=TRUE)
      setkey(L.selected.dt)
      rm(ll.dt)
    }
  }
  col.order <- c(1:4,c(7,5:6),8:10)
  setcolorder(L.selected.dt,col.order)
  if (data.name=="lalonde") {
    L.selected.dt[, Empty := NULL]
  } else if (grepl("rhc",data.name)) {
    L.selected.dt[, c("CovSel","Empty") := NULL]
    if (file.exists("rhc-varnames.csv")) {
      rhc_varnames <- read.csv2("rhc-varnames.csv",as.is=TRUE)
      rhc_varnames.full <- rep(NA,nrow(L.selected.dt))
      for (vn in 1:nrow(L.selected.dt)) {
        rhc_varnames.full[vn] <- 
        rhc_varnames[rhc_varnames$labelname==
                       L.selected.dt[,as.character(L.names)][vn],
                     "tablename"]
      }
      L.selected.dt[ ,L.names := rhc_varnames.full]
    } 
  }
  print(xtable(L.selected.dt[,-c(1:2)]),include.rownames=FALSE)
  print(colSums(L.selected.dt[,-(1:3),with=FALSE]))
  print(unlist(PSsizes)[(col.order-3)[-(1:3)]])
  print(unlist(res))
  print(ATE.selected)

  # plots ===================================================================
  diags.fs <- data.frame(diags.FwdSel$ds_stable[[1]])
  diags.fs <- cbind("j"=1:nrow(diags.fs),diags.fs)
  diags.fs <- cbind(diags.fs, 
                    "std.dif"=diags.fs[,"dif"]/sqrt(diags.fs[,"var.dif"]))
  diags.fs[nrow(diags.fs),"std.dif"] <- 0
  plot_filename <- paste0("plots-",data.name,"-diff_traj.pdf")
  pdf(plot_filename,width=8,height=5)
  plot(diags.fs$std.dif,
       main="Std. diff. between treatment effects",
       xlab="Orbit j",ylab="Std. diff.",
       pch=4)
  abline(h=0,lty=3)
  if (max(abs(range(diags.fs$std.dif))) > 2) {
    abline(h=c(-1,1)*qnorm(.975),lty=2)  
  }
  if (grepl("rhc",data.name)) {
    lines(locfit(std.dif~lp(j,nn=0.3),data=diags.fs))
  } else {
    lines(locfit(std.dif~lp(j),data=diags.fs))  
  }
  
  points(PSsizes["PSsizes.ds_stable"],
         diags.fs$std.dif[PSsizes["PSsizes.ds_stable"]],
         pch=17, cex=1.5)
  dev.off()

  # Cochran's Q
  plot_filename <- paste0("plots-",data.name,"-Q.pdf")
  pdf(plot_filename,width=8,height=5)
  if (grepl("rhc",data.name)) {
    windowsize <- as.integer(strsplit(data.name,"rhc")[[1]][2])
  } else {
    windowsize <- 3
  }
  Qvals <- c(rep(NA,(windowsize-1)/2), diags.FwdSel$ds_stable[[3]],
             rep(NA,(windowsize-1)/2))  
  plot(log10(Qvals), 
       main="Cochran's Q",
       xlab="Orbit j",ylab=expression(paste(log[10],"(Q)")),
       pch=4)
  points(PSsizes["PSsizes.ds_stable"],
         log10(Qvals[PSsizes["PSsizes.ds_stable"]]),
         pch=17, cex=1.5)
  dev.off()
}
