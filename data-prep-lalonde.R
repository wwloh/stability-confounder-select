# observational study
data(lalonde)
dim(lalonde) # number of obs
summary(lalonde)
sapply(lalonde, class)

lalonde.dt <- data.table(sapply(lalonde, as.numeric))
setkey(lalonde.dt)

# covariates
l <- lalonde.dt[,!(colnames(lalonde.dt) %in% c("treat","re78")),with=FALSE]
Lnames <- data.frame(cbind("L.idx"=1:ncol(l),"L.names"=names(l)))
l <- as.data.frame(l)
colnames(l) <- NULL
mydata <- data.frame("i"=1:nrow(l),
                     "L"=l,
                     "treat"=as.integer(lalonde.dt$treat),
                     "Y"=as.integer(lalonde.dt$re78))

# ITT effect
itt <- sum(lalonde.dt$re78*lalonde.dt$treat, na.rm=TRUE)/sum(lalonde.dt$treat) - 
  sum(lalonde.dt$re78*(1-lalonde.dt$treat),na.rm=TRUE)/sum(1-lalonde.dt$treat)
