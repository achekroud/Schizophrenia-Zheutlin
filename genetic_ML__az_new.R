# Schizophrenia-linked genetic influences on cognition
# Author: AMC & ABZ, Nov 2015

#### Housekeeping
# Data manipulation and plotting
library("dplyr"); library("ggplot2"); library("RColorBrewer")

# Statistical packages
library("glmnet") # https://urldefense.proofpoint.com/v2/url?u=http-3A__web.stanford.edu_-7Ehastie_glmnet_glmnet-5Falpha.html&d=AwIGAg&c=-dg2m7zWuuDZ0MUcV7Sdqw&r=1eBIqVAhmnkZRWlCX5iU07Sqp55OtRSRubwrZkpK9R4&m=x0sYKnmPdCd6tz6RllAyYboYJJPTi-LQXgCBMesayGo&s=W3iRrgwrQgyU6pVDZrYbPv8yYV80ZJeVa7JpwoHDfqM&e= 
library("e1071"); library("caret"); library("pROC")
library("permute"); library("gbm"); library("klaR")
library("dplyr")

# Parallel Computing 
library("doMC"); registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
workDir <- "~/Google Drive/ABZ/Yale/ML_PGC-SZ/"
setwd(workDir)
#dataDir <- paste0(workDir, "/dfs/") 

# Load custom functions from adam's library
source('AMC_functions.R')


### Read in pre-processed data 
# sweden_eQTL <- read.table(paste0(dataDir,"eQTL-snps-genes.txt"), header=TRUE, as.is=TRUE)
# sweden_snp  <- read.table(paste0(dataDir,"GWAS-snps-swe-input.txt"), header=TRUE, as.is=TRUE)
# cnp         <- read.table(paste0(dataDir,"GWAS-snps-cnp-input.txt"), header=TRUE, as.is=TRUE)
sweden_snp  <- read.table(paste0(workDir,"GWASinput2-swe.txt"), header=TRUE, as.is=TRUE)
cnp         <- read.table(paste0(workDir,"GWASinput2-cnp.txt"), header=TRUE, as.is=TRUE)

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
targets$VR2DR_TOTALRAW <- targets$VR2DR_TOTALRAW +1 #someone scored 0 and you can't log 0 and it isnt normal so we added 1 to all scores (range 0-43)

# get residuals
# resid <- lapply(outcomes, function(x) {
#   lm(eval(substitute(i ~ age + gender, list(i = as.name(x)))),data = cnp.df)
# })
# 
# resid_vals <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
# names(resid_vals) <- outcomes


## Scale variables manually (so that we can apply these params to replication set)
#  record m/sd
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)

#  apply scaling
#  log before scaling
log.targets   <- apply(targets, 2, log) %>% as.data.frame # alt:log transform 
targets       <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame
# s.log.targets <- apply(log.targets, 2, scale) %>% as.data.frame

# rep.outcomes <- apply(rep.outcomes, 2, function(i) psych::winsor(i, trim=0.025))



### Machine Learning preparation
## Set up cross validation
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 3,
                           number = 10,
                           savePredictions = TRUE)

## Hyperparameter grid search
# gbmGrid    <- expand.grid(.interaction.depth = c(1,2),
#                           .n.trees = seq(10,1000,by=50),
#                           .shrinkage = c(0.01),
#                           .n.minobsinnode = 5)
rfGrid     <- expand.grid(.mtry = c(3,4,10)) 

# Train models

### just checking which models are now available for regression/dual use
t <- getModelInfo()
m <- list();
for (i in names(t)){
  if (t[[i]]$type != "Classification"){
    m <- c(m, t[i])
  }
}


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
glm.loop.log <- lapply(log.targets, function(i) glmPipeline(i))
# glm.loop.log.s <- lapply(s.log.targets, function(i) glmPipeline(i))
glm.models   <- lapply(glm.loop, function(i) i[[1]])
glm.perfs    <- lapply(glm.loop, function(i) i[[2]])
glm.models.log  <- lapply(glm.loop.log, function(i) i[[1]])
glm.perfs.log   <- lapply(glm.loop.log, function(i) i[[2]])
# glm.models.log.s  <- lapply(glm.loop.log.s, function(i) i[[1]])
# glm.perfs.log.s   <- lapply(glm.loop.log.s, function(i) i[[2]])


# gbmPipeline <- function(response, Xmat = predictors, grid = gbmGrid, cvpar = fitControl,...){
#     set.seed(1)
#     model <- train(x = Xmat, y = response,
#                    method = "gbm",
#                    metric = "Rsquared",
#                    tuneGrid = grid,
#                    trControl = cvpar)
#     perf <- getTrainPerf(model)
#     return(list(model,perf))
# }
# 
# gbm.loop     <- lapply(targets, function(i) gbmPipeline(i))
# gbm.loop.log <- lapply(log.targets, function(i) gbmPipeline(i))
# gbm.models   <- lapply(gbm.loop, function(i) i[[1]])
# gbm.perfs    <- lapply(gbm.loop, function(i) i[[2]])
# gbm.models.log  <- lapply(gbm.loop.log, function(i) i[[1]])
# gbm.perfs.log   <- lapply(gbm.loop.log, function(i) i[[2]])

rfPipeline <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "rf",
                 metric = "Rsquared",
                 tuneGrid = grid,
                 trControl = cvpar,
                 importance=TRUE)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}


rf.loop     <- lapply(targets, function(i) rfPipeline(i))
rf.loop.log <- lapply(log.targets, function(i) rfPipeline(i))
# rf.loop.log.s <- lapply(s.log.targets, function(i) rfPipeline(i))
rf.models   <- lapply(rf.loop, function(i) i[[1]])
rf.perfs    <- lapply(rf.loop, function(i) i[[2]])
rf.models.log  <- lapply(rf.loop.log, function(i) i[[1]])
rf.perfs.log   <- lapply(rf.loop.log, function(i) i[[2]])
# rf.models.log.s  <- lapply(rf.loop.log.s, function(i) i[[1]])
# rf.perfs.log.s   <- lapply(rf.loop.log.s, function(i) i[[2]])

plot(varImp(rf.models$CRT_TIME1, scale=TRUE), top=15)


# crt1.out  <- promising[[1]]
# crt1.vars <- varImp(crt1.out[[1]], scale = FALSE)[[1]]
# crt1.vars <- cbind(crt1.vars, rownames(crt1.vars)) %>% as.data.frame
# crt2.out  <- promising[[2]]
# crt2.vars <- varImp(crt2.out[[1]], scale = FALSE)[[1]]
# crt2.vars <- cbind(crt2.vars, rownames(crt2.vars)) %>% as.data.frame
# 
# c1 <- crt1.vars %>% dplyr::filter(Overall > 0)
# c2 <- crt2.vars %>% dplyr::filter(Overall > 0)

# write.csv(c1,file="crt1.csv")
# write.csv(c2,file="crt2.csv")
# 
# dplyr::setdiff(c1$`rownames(crt1.vars)`, c2$`rownames(crt1.vars)`)






#### External Replication
raw_rep <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),]

rep.predictors <- raw_rep %>% dplyr::select(one_of(snps)) %>% as.matrix
rep.outcomes   <- raw_rep %>% dplyr::select(one_of(outcomes))
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW +1

# multiply outcomes to fit CNP scales
# vocab: CNP = out of 50, Swedish = out of 40; 50/40 = 1.25
# digit symbol: CNP = out of 48, Swedish = out of 30; 48/30 = 1.6
# VR I and II: CNP = out of 43, Swedish = out of 108; 43/108 = .398
rep.outcomes$VOC_TOTALRAW <- rep.outcomes$VOC_TOTALRAW * 1.25
rep.outcomes$DS_TOTALRAW <- rep.outcomes$DS_TOTALRAW * 1.6
rep.outcomes$VR1IR_TOTALRAW <- rep.outcomes$VR1IR_TOTALRAW * .398
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW * .398

# twin1 <- raw_rep %>% dplyr::filter(tvab == 1) #twin specific
# twin2 <- raw_rep %>% dplyr::filter(tvab == 2)
# twin1.pred <- twin1 %>% dplyr::select(one_of(snps)) %>% as.matrix
# twin2.pred <- twin2 %>% dplyr::select(one_of(snps)) %>% as.matrix
# twin1.out.log  <- twin1 %>% dplyr::select(CRT_TIME1, CRT_TIME2) %>%
#   apply(2, log) %>% as.data.frame
# twin2.out.log  <- twin2 %>% dplyr::select(CRT_TIME1, CRT_TIME2) %>%
#   apply(2, log) %>% as.data.frame

# Apply scaling transformations from CNP
rep.outcomes.log <- apply(rep.outcomes, 2, log)
rep.outcomes.s   <- (rep.outcomes - means)/sds
rep.outcomes.s   <- apply(rep.outcomes.s, 2, function(i) psych::winsor(i, trim=0.01))
# rep.outcomes.log.s <- (rep.outcomes.log - means)/sds

# scale within sample
# means_swe    <- apply(rep.outcomes,2,mean,na.rm=T)
# sds_swe      <- apply(rep.outcomes,2,sd,na.rm=T)
# 
# rep.outcomes.s2   <- (rep.outcomes - means_swe)/sds_swe
# rep.outcomes.s2   <- apply(rep.outcomes.s2, 2, function(i) psych::winsor(i, trim=0.01))
# rep.outcomes.log.s2 <- (rep.outcomes.log - means_swe)/sds_swe

# rep.outcomes <- apply(rep.outcomes, 2, function(i) psych::winsor(i, trim=0.03))

# gbm.pred.tr     <- lapply(gbm.models, function(i) predict(i, newdata=rep.predictors))
# gbm.pred.log.tr <- lapply(gbm.models.log, function(i) predict(i, newdata=rep.predictors))
rf.pred.tr        <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
rf.pred.log.tr    <- lapply(rf.models.log, function(i) predict(i, newdata=rep.predictors))
# rf.pred.log.s.tr  <- lapply(rf.models.log.s, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr       <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))
glm.pred.log.tr   <- lapply(glm.models.log, function(i) predict(i, newdata=rep.predictors))
# glm.pred.log.s.tr <- lapply(glm.models.log.s, function(i) predict(i, newdata=rep.predictors))

out.cors <- list()
for (i in 1:7){
   glm.cor <- cor(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
   rf.cor  <- cor(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
   glm.log.cor <- cor(glm.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")
   rf.log.cor  <- cor(rf.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")
#   glm.log.s.cor <- cor(glm.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")
#   rf.log.s.cor  <- cor(rf.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")
#   glm.cor2 <- cor(glm.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")
#   rf.cor2  <- cor(rf.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")
#   glm.log.s.cor2 <- cor(glm.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")
#   rf.log.s.cor2  <- cor(rf.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")
#   out.cors[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor, glm.log.s.cor, 
#                      rf.log.s.cor, glm.cor2, rf.cor2, glm.log.s.cor2, rf.log.s.cor2)
out.cors[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor)
}
names(out.cors) <- names(rep.outcomes)

out.ps <- list()
for (i in 1:7){
   glm.cor <- cor.test(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
   rf.cor  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
   glm.log.cor <- cor.test(glm.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")$p.value
   rf.log.cor  <- cor.test(rf.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")$p.value
#   glm.log.s.cor <- cor.test(glm.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")$p.value
#   rf.log.s.cor  <- cor.test(rf.pred.log.s.tr[[i]], rep.outcomes.log.s[,i], use="pairwise")$p.value
#   glm.cor2 <- cor.test(glm.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")$p.value
#   rf.cor2  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")$p.value
#   glm.log.s.cor2 <- cor.test(glm.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")$p.value
#   rf.log.s.cor2  <- cor.test(rf.pred.log.s.tr[[i]], rep.outcomes.log.s2[,i], use="pairwise")$p.value
#   out.ps[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor, glm.log.s.cor, 
#                    rf.log.s.cor, glm.cor2, rf.cor2, glm.log.s.cor2, rf.log.s.cor2)
out.ps[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor)
}
names(out.ps) <- names(rep.outcomes)

#dsz <- ifelse(raw_rep$dx %in% c(1,3), "pat", "con")

# correlations controlling for FID for models that work
# library(lme4)
library(nlme)

mix <- cbind(rep.outcomes.log[,7], rep.outcomes.s[,3], rf.pred.log.tr[[7]], rf.pred.tr[[3]], raw_rep[,c(2,6,7)]) %>% as.data.frame
names(mix) <- c("out.vr2.log","out.tr1", "pred.vr.log", "pred.tr1", "fid","age","sex")

# mx.tr1 <- lmer(out.tr1.log ~ pred.tr1 + (1 | fid), data=mix)
# mx.voc <- lmer(outvoc ~ pred.voc + (1 | fid), data=mix)
# mx.vr <- lmer(outvr ~ pred.vr + (1 | fid), data=mix)

mx.tr1 <- nlme::lme(out.tr1 ~ pred.tr1 + age + sex, random = ~1 | fid, data = mix, na.action = na.omit)
# mx.voc <- nlme::lme(outvoc ~ pred.voc + age + sex, random = ~1 |fid, data = mix, na.action=na.omit)
mx.vr <- nlme::lme(out.vr2.log ~ pred.vr.log + age + sex, random = ~1 | fid, data = mix, na.action = na.omit)

# to load later
save.image(file=paste0(workDir,"amanda.results4.Rdata"))


# cor(gbm.pred.tr[[1]], rep.outcomes[,1])
# cor(rf.pred.tr[[1]], rep.outcomes[,1])
# cor(gbm.pred.log.tr[[1]], rep.outcomes.log[,1])
# cor(rf.pred.log.tr[[1]], rep.outcomes.log[,1])
# 
# cor(gbm.pred.tr[[2]], rep.outcomes[,2])
# cor(rf.pred.tr[[2]], rep.outcomes[,2])
# cor(gbm.pred.log.tr[[2]], rep.outcomes.log[,2])
# cor(rf.pred.log.tr[[2]], rep.outcomes.log[,2])


# Twin specific
# t1.gbm.pred.log.tr <- lapply(list(gbm.models.log[[5]],gbm.models.log[[6]]),
#                              function(i) predict(i, newdata=twin1.pred))
# t1.rf.pred.log.tr  <- lapply(list(rf.models.log[[5]],rf.models.log[[6]]),
#                              function(i) predict(i, newdata=twin1.pred))
# t2.gbm.pred.log.tr <- lapply(list(gbm.models.log[[5]],gbm.models.log[[6]]),
#                              function(i) predict(i, newdata=twin2.pred))
# t2.rf.pred.log.tr  <- lapply(list(rf.models.log[[5]],rf.models.log[[6]]),
#                              function(i) predict(i, newdata=twin2.pred))
# cor(t1.gbm.pred.log.tr[[1]], twin1.out.log[,1])
# cor(t1.rf.pred.log.tr[[1]],  twin1.out.log[,1])
# cor(t2.gbm.pred.log.tr[[1]], twin2.out.log[,1])
# cor(t2.rf.pred.log.tr[[1]],  twin2.out.log[,1])
# 
# 
# 
# rf.varimp.log.t1  <- lapply(rf.models.log, varImp)
# gbm.varimp.log.t1 <- lapply(gbm.models.log, varImp)
# plot(rf.varimp.log.t1[[5]])
# plot(gbm.varimp.log.t1[[5]])
# 
