# Schizophrenia-linked genetic influences on cognition
# Author: AMC & ABZ, Nov 2015

#### Housekeeping
# Data manipulation and plotting
library("dplyr"); library("ggplot2"); library("RColorBrewer")

# Statistical packages
library("glmnet"); library("nlme")
library("e1071"); library("caret"); library("pROC")
library("permute"); library("gbm"); library("klaR")
library("dplyr")

# Parallel Computing 
library("doMC"); registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
workDir <- "/Users/adamchekroud/Desktop/ABZ"
setwd(workDir)
dataDir <- paste0(workDir, "/dfs/") 

# Load custom functions from adam's library
source('/Users/adamchekroud/functions/functions.R')


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

# get residuals
resid <- lapply(outcomes, function(x) {
  lm(eval(substitute(i ~ age + gender, list(i = as.name(x)))), data = cnp.df)
})

resid_vals <- lapply(resid, function(x) residuals(x)) %>% as.data.frame()
names(resid_vals) <- outcomes

targets <- resid_vals

# targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()




## Scale variables manually (so that we can apply these params to replication set)
#  record m/sd
# means    <- apply(targets,2,mean)
# sds      <- apply(targets,2,sd)
# 
# #  apply scaling
# # log.targets   <- apply(targets, 2, log) %>% as.data.frame # alt:log transform 
# targets       <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame

# rep.outcomes <- apply(rep.outcomes, 2, function(i) psych::winsor(i, trim=0.025))



### Machine Learning preparation
## Set up cross validation
# parallel seeds
cv.rep  = 3
cv.fold = 10
set.seed(1)
seeds <- vector(mode = "list", length = cv.rep*cv.fold + 1)
for(i in 1:(length(seeds)-1)) seeds[[i]] <- sample.int(1000, cv.rep*cv.fold)
seeds[[length(seeds)]] <- sample.int(1000,1)

fitControl <- trainControl(method = "repeatedcv",
                           repeats = cv.rep,
                           number = cv.fold,
                           seeds = seeds,
                           savePredictions = TRUE)

## Hyperparameter grid search
gbmGrid    <- expand.grid(.interaction.depth = c(1,2),
                          .n.trees = seq(300,2000,by=100),
                          .shrinkage = c(0.005),
                          .n.minobsinnode = c(5,10))
# rfGrid     <- expand.grid(.mtry = c(3,4,10)) 

# Train models

### just checking which models are now available for regression/dual use
# t <- getModelInfo()
# m <- list();
# for (i in names(t)){
#   if (t[[i]]$type != "Classification"){
#     m <- c(m, t[i])
#   }
# }


glmPipeline <- function(response, Xmat = predictors, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "glm",
                 metric = "Rsquared",
                 trControl = cvpar)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}

glm.loop        <- lapply(targets, function(i) glmPipeline(i))
# glm.loop.log    <- lapply(log.targets, function(i) glmPipeline(i))
glm.models      <- lapply(glm.loop, function(i) i[[1]])
glm.perfs       <- lapply(glm.loop, function(i) i[[2]])
# glm.models.log  <- lapply(glm.loop.log, function(i) i[[1]])
# glm.perfs.log   <- lapply(glm.loop.log, function(i) i[[2]])





rfPipeline <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
  set.seed(1)
  model <- train(x = Xmat, y = response,
                 method = "gbm",
                 metric = "Rsquared",
                 tuneGrid = gbmGrid,
                 trControl = cvpar)
  perf <- getTrainPerf(model)
  return(list(model,perf))
}


rf.loop        <- lapply(targets, function(i) rfPipeline(i))
# rf.loop.log    <- lapply(log.targets, function(i) rfPipeline(i))
rf.models      <- lapply(rf.loop, function(i) i[[1]])
rf.perfs       <- lapply(rf.loop, function(i) i[[2]])
# rf.models.log  <- lapply(rf.loop.log, function(i) i[[1]])
# rf.perfs.log   <- lapply(rf.loop.log, function(i) i[[2]])


# plot(varImp(rf.models$CRT_TIME1, scale=TRUE), top=15)


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

# multiply outcomes to fit CNP scales
# vocab: CNP = out of 50, Swedish = out of 40; 50/40 = 1.25
# digit symbol: CNP = out of 48, Swedish = out of 30; 48/30 = 1.6
# VR I and II: CNP = out of 43, Swedish = out of 108; 43/108 = .398
rep.outcomes$VOC_TOTALRAW <- rep.outcomes$VOC_TOTALRAW * 1.25
rep.outcomes$DS_TOTALRAW <- rep.outcomes$DS_TOTALRAW * 1.6
rep.outcomes$VR1IR_TOTALRAW <- rep.outcomes$VR1IR_TOTALRAW * .398
rep.outcomes$VR2DR_TOTALRAW <- rep.outcomes$VR2DR_TOTALRAW * .398


# Apply scaling transformations from CNP
# rep.outcomes.log <- apply(rep.outcomes, 2, log)
rep.outcomes.s   <- rep.outcomes
# rep.outcomes.s   <- (rep.outcomes - means)/sds
rep.outcomes.s   <- apply(rep.outcomes.s, 2, function(i) psych::winsor(i, trim=0.01))


# gbm.pred.tr     <- lapply(gbm.models, function(i) predict(i, newdata=rep.predictors))
# gbm.pred.log.tr <- lapply(gbm.models.log, function(i) predict(i, newdata=rep.predictors))
rf.pred.tr        <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
# rf.pred.log.tr    <- lapply(rf.models.log, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr       <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))
# glm.pred.log.tr   <- lapply(glm.models.log, function(i) predict(i, newdata=rep.predictors))

out.cors <- list()
for (i in 1:7){
   glm.cor <- cor(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
   rf.cor  <- cor(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
   # glm.log.cor <- cor(glm.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")
   # rf.log.cor  <- cor(rf.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")
#   glm.cor2 <- cor(glm.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")
#   rf.cor2  <- cor(rf.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")
#   out.cors[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor, glm.cor2, rf.cor2)
out.cors[[i]] <- c(glm.cor, rf.cor) #, glm.log.cor, rf.log.cor)
}
names(out.cors) <- names(rep.outcomes)

out.ps <- list()
for (i in 1:7){
   glm.cor <- cor.test(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
   rf.cor  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
#    glm.log.cor <- cor.test(glm.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")$p.value
#    rf.log.cor  <- cor.test(rf.pred.log.tr[[i]], rep.outcomes.log[,i], use="pairwise")$p.value
#   glm.cor2 <- cor.test(glm.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")$p.value
#   rf.cor2  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s2[,i], use="pairwise")$p.value
#   out.ps[[i]] <- c(glm.cor, rf.cor, glm.log.cor, rf.log.cor, glm.cor2, rf.cor2)
out.ps[[i]] <- c(glm.cor, rf.cor) #, glm.log.cor, rf.log.cor)
}
names(out.ps) <- names(rep.outcomes)

#dsz <- ifelse(raw_rep$dx %in% c(1,3), "pat", "con")

# correlations controlling for FID for models that work


mix.df <- cbind(rep.outcomes.s, raw_rep[,c(2,6,7)]) %>% as.data.frame

mixed.mods <- lapply(1:7, function(x) {
    mix.df <- cbind(rep.outcomes.s[,x], rf.pred.tr[[x]], raw_rep[,c(2,6,7)]) %>% as.data.frame
    names(mix.df)[1:2] <- c("response", "prediction")
    nlme::lme(response ~ prediction + age + gender, random = ~1 | FID, na.action = na.omit, data = mix.df)
})
names(mixed.mods) <- names(rep.outcomes)






# to load later
save.image(file=paste0(workDir,"amanda.results4.Rdata"))





# # gbmPipeline <- function(response, Xmat = predictors, grid = gbmGrid, cvpar = fitControl,...){
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
