# randomized trial
data(ACTG175)
table(ACTG175$treat,ACTG175$arms)
### recode treatment 
ACTG175$treat <- ifelse(ACTG175$arms==0 | ACTG175$arms==3, 0, 1)
table(ACTG175$treat,ACTG175$arms)

nrow(ACTG175)
### a dichotomous response is created with missing values maintained
ACTG175$cd496bin <- ifelse(ACTG175$cd496 > 250, 1, 0)
complete_cases <- which(!is.na(ACTG175$cd496bin))
length(complete_cases)

### candidate covariates 
l <- ACTG175[complete_cases,
             c("age","wtkg","hemo","homo","drugs","karnof","oprior","z30",
               "zprior","preanti","race","gender","str2","strat","symptom",
               "cd40","cd80")]
### any singular covariates
which(apply(l,2,sd)==0)
l <- l[, -which(apply(l,2,sd)==0)]
any(apply(l,2,sd)==0)

(Lnames <- data.frame(cbind("L.idx"=1:ncol(l),"L.names"=names(l))))
summary(l)
l <- as.data.frame(l)
colnames(l) <- NULL
mydata <- data.frame("i"=1:nrow(l),
                     "L"=l,
                     "treat"=as.integer(ACTG175$treat[complete_cases]),
                     "Y"=as.integer(ACTG175$cd496bin[complete_cases]))

# ITT effect
itt <- sum(ACTG175$cd496bin*ACTG175$treat, na.rm=TRUE)/sum(ACTG175$treat) - 
  sum(ACTG175$cd496bin*(1-ACTG175$treat),na.rm=TRUE)/sum(1-ACTG175$treat)
