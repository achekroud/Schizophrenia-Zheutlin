# Schizophrenia-linked genetic influences on cognition
# Author: AMC & ABZ, Nov 2015

#### Housekeeping
libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "pROC", "permute", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "readxl", "psych", "e1071", "dplyr", "nlme")
lapply(libs, require, character.only = TRUE)
registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/data/ML_genetics/Schizophrenia-Zheutlin")
workDir <- getwd()
dataDir <- paste0(workDir, "/dfs/") 

# Load custom functions from adam's library
source("/data/ML_genetics/functions/functions.R")


### Read in pre-processed data 
sweden_snp  <- read.table(paste0(dataDir,"GWASinput2-swe.txt"), header=TRUE, as.is=TRUE)
cnp         <- read.table(paste0(dataDir,"GWASinput2-cnp.txt"), header=TRUE, as.is=TRUE)

# exclude symbol / spatial span
sweden_snp <- sweden_snp[,c(1:11,13:92)]
cnp <- cnp[,c(1:8,10:89)]

# Names of snps and outcomes
snps       <- names(cnp)[12:88] 
outcomes   <- names(cnp)[5:11]

# Rename swedish outcomes to match cnp
names(sweden_snp)[8:14] <- outcomes 
# Some proxy snps used, need to rename
names(sweden_snp)[71:73] <- c("rs2535627", "rs4388249", "rs2905426")

# Remove subjects with missing outcomes in CNP
# na_count <- sapply(cnp, function(y) sum(length(which(is.na(y)))))
# na_count <- data.frame(na_count)

cnp.df     <- cnp[(complete.cases(cnp$CRT_TIME1)),]
swe.df     <- sweden_snp 

# set predictors and targets
predictors <- cnp.df %>% dplyr::select(one_of(snps)) %>% as.matrix() 
targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()

## Scale variables manually (so that we can apply these params to replication set)
#  record m/sd
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)

#  apply scaling
targets       <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame

# get residuals
# targets <- cbind(targets, cnp.df[,2:3])
# resid <- lapply(outcomes, function(x) {
#   lm(eval(substitute(i ~ age + gender, list(i = as.name(x)))),data = cnp.df)
# })
# 
# resid_vals <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
# names(resid_vals) <- outcomes
# 
# targets <- resid_vals

### Machine Learning preparation
## Set up cross validation
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 3,
                           number = 10,
                           savePredictions = TRUE,
                           selectionFunction = "oneSE")

## Hyperparameter grid search
gbmGrid    <- expand.grid(.interaction.depth = c(1,2),
                          .n.trees = seq(10,4000,by=25),
                          .shrinkage = c(0.01, 0.05),
                          .n.minobsinnode = 5)
rfGrid     <- expand.grid(.mtry = c(3,4,10)) 



##### ##### ##### ##### ##### 
##### Model training
##### ##### ##### ##### ##### 

glmPipeline <- function(response, Xmat = predictors, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "glm",
                 metric = "Rsquared",
                 trControl = cvpar)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}

glm.loop     <- lapply(targets, function(i) glmPipeline(i))
glm.models   <- lapply(glm.loop, function(i) i[[1]])
glm.perfs    <- lapply(glm.loop, function(i) i[[2]])



gbmPipeline <- function(response, Xmat = predictors, grid = gbmGrid, cvpar = fitControl,...){
    set.seed(1)
    model <- train(x = Xmat, y = response,
                   method = "gbm",
                   metric = "Rsquared",
                   tuneGrid = grid,
                   trControl = cvpar)
    perf <- getTrainPerf(model)
    return(list(model,perf))
}

gbm.loop     <- lapply(targets, function(i) gbmPipeline(i))
gbm.models   <- lapply(gbm.loop, function(i) i[[1]])
gbm.perfs    <- lapply(gbm.loop, function(i) i[[2]])

rfPipeline <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "rf",
                 metric = "Rsquared",
                 tuneGrid = grid,
                 trControl = cvpar,
                 importance = TRUE)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}


rf.loop     <- lapply(targets, function(i) rfPipeline(i))
rf.models   <- lapply(rf.loop, function(i) i[[1]])
rf.perfs    <- lapply(rf.loop, function(i) i[[2]])


plot(varImp(rf.models$CRT_TIME1, scale=TRUE), top=15)


# dplyr::setdiff(c1$`rownames(crt1.vars)`, c2$`rownames(crt1.vars)`)






#### External Replication
raw_rep <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),] #rm NAs

rep.predictors <- raw_rep %>% dplyr::select(one_of(snps)) %>% as.matrix
rep.outcomes   <- raw_rep %>% dplyr::select(one_of(outcomes))

# multiply outcomes to fit CNP scales
# vocab: CNP = out of 50, Swedish = out of 40; 50/40 = 1.25
# digit symbol: CNP = out of 48, Swedish = out of 30; 48/30 = 1.6
# VR I and II: CNP = out of 43, Swedish = out of 108; 43/108 = .398
rep.outcomes$VOC_TOTALRAW   <- rep.outcomes$VOC_TOTALRAW * 1.25
rep.outcomes$DS_TOTALRAW    <- rep.outcomes$DS_TOTALRAW * 1.6
rep.outcomes$VR1IR_TOTALRAW <- rep.outcomes$VR1IR_TOTALRAW * .398
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW * .398



# Apply scaling transformations from CNP
rep.outcomes.s   <- (rep.outcomes - means)/sds
rep.outcomes.s   <- apply(rep.outcomes.s, 2, function(i) psych::winsor(i, trim=0.01))


gbm.pred.tr     <- lapply(gbm.models, function(i) predict(i, newdata=rep.predictors))
rf.pred.tr        <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr       <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))

out.cors <- list()
for (i in 1:7){
  glm.cor <- cor(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
  rf.cor  <- cor(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
  gbm.cor <- cor(gbm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
  out.cors[[i]] <- c(glm.cor, gbm.cor, rf.cor) 
}
names(out.cors) <- names(rep.outcomes)

out.ps <- list()
for (i in 1:7){
  glm.cor <- cor.test(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
  rf.cor  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
  gbm.cor <- cor.test(gbm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
  out.ps[[i]] <- c(glm.cor, gbm.cor, rf.cor)
}
names(out.ps) <- names(rep.outcomes)


##### ##### ##### ##### ##### 
##### Mixed effects models controlling for FID
##### ##### ##### ##### ##### 

mix <- cbind(rep.outcomes.s, rf.pred.tr, glm.pred.tr, raw_rep[,c("FID", "age", "sex")]) %>% as.data.frame

outcomes.rf <- NULL  # start empty list
outcomes.glm <- NULL
for (i in 1:7){outcomes.rf[i] <- paste0(outcomes[i],".rf")}  # rename cols to prevent clash later
for (i in 1:7){outcomes.glm[i] <- paste0(outcomes[i],".glm")}

names(mix)[1:7]   <- outcomes            # renaming using the list of names
names(mix)[8:14]  <- outcomes.rf
names(mix)[15:21] <- outcomes.glm
names(mix)[22:24] <- c("fid","age","sex")

outstats <- matrix(nrow=7,ncol=5)
for (i in 1:length(outcomes)){
  f1 <- formula(paste0(outcomes[i],"~",outcomes[i],".rf","+ age + sex"))
  f2 <- formula(paste0(outcomes[i],"~",outcomes[i],".glm","+ age + sex"))
  m1 <- lme(fixed=f1, random=~1|fid,data=mix,na.action=na.omit)
  m2 <- lme(fixed=f2, random=~1|fid,data=mix,na.action=na.omit)
  outstats[i,1] <- outcomes[i]
  outstats[i,2] <- summary(m1)$tTable[2,4]
  outstats[i,3] <- summary(m1)$tTable[2,5]
  outstats[i,4] <- summary(m2)$tTable[2,4]
  outstats[i,5] <- summary(m2)$tTable[2,5]
}

outstats <- as.data.frame(outstats)
names(outstats) <- c("var","tval.rf","pval.rf","tval.glm","pval.glm")

# CNP correlations and p-values
glm.r2     <- unlist(lapply(glm.perfs, function(x) x$TrainRsquared))
rf.r2      <- unlist(lapply(rf.perfs, function(x) x$TrainRsquared))

r2   <- rbind(glm.r2,unname(rf.r2))
rval <- apply(r2,2,sqrt)
tval <- apply(r2,2,function(x) sqrt((737*x)/(1-x)))
pval <- apply(tval,2,function(x) 2*(1-abs(pt(x,df=737))))

# to load later
save.image(file=paste0(workDir,"amanda.results4.Rdata"))

#################### figure
library(ggplot2); library(stringr)

# snp-gene table
code       <- read.csv("scz2.anneal.108.genes.csv")
proxy      <- read.table("index-proxy-table_allsnps.txt",header=T)
proxy$gene <- code[match(proxy$index,code$bestsnp),8]

# trails1
fig_trails1        <- caret::varImp(rf.models$CRT_TIME1, scale=TRUE)[[1]] %>% as.data.frame
fig_trails1        <- cbind(rownames(fig_trails1), fig_trails1) %>% arrange(-Overall)
names(fig_trails1) <- c("SNP", "Importance")
fig_trails1$Gene   <- proxy[ match(fig_trails1$SNP, proxy$proxy), 3]
fig_trails1$SNP    <- factor(fig_trails1$SNP,
                             levels = fig_trails1$SNP[order(fig_trails1$Importance)])

# visual reproduction - delayed recall
fig_vrdr           <- caret::varImp(rf.models$VR2DR_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vrdr           <- cbind(rownames(fig_vrdr), fig_vrdr) %>% arrange(-Overall)
names(fig_vrdr)    <- c("SNP", "Importance")
fig_vrdr$Gene      <- proxy[match(fig_vrdr$SNP, proxy$proxy), 3]
fig_vrdr$SNP       <- factor(fig_vrdr$SNP,
                             levels = fig_vrdr$SNP[order(fig_vrdr$Importance)])

# code by shared vs. task-specific genes
fig_trails1$shared         <- ifelse(fig_trails1$Gene %in% fig_vrdr[1:20,3],1,2) %>% as.factor()
fig_vrdr$shared            <- ifelse(fig_vrdr$Gene %in% fig_trails1[1:20,3],1,2) %>% as.factor()
levels(fig_trails1$shared) <- c("both tasks","task-specific")
levels(fig_vrdr$shared)    <- c("both tasks","task-specific")

ggplot(fig_trails1[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_trails1[1:20,3]), width=60)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=13, face="bold"),
        axis.text = element_text(size=9,color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(50,100))

ggplot(fig_vrdr[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.75) +
  scale_x_discrete(labels= function(x) str_wrap(rev(fig_vrdr[1:20,3]), width=60)) +
  scale_fill_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=13, face="bold"),
        axis.text = element_text(size=9,color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black")) +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(50,100))

