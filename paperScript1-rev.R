# Schizophrenia-linked genetic influences on cognition

#### Housekeeping
libs <- c("ggplot2", "RColorBrewer", "glmnet", "caret", "pROC", "permute", "gbm", "randomForest",
          "plyr", "foreach", "doMC", "psych", "e1071", "dplyr", "nlme", 
          "psych","stringr")
invisible(lapply(libs, require, character.only = TRUE))
registerDoMC(detectCores()-1)

### Directories and seeds
set.seed(1) # Reading from the same sheet
setwd("/data/ML_genetics/Schizophrenia-Zheutlin")
workDir <- getwd()
dataDir <- paste0(workDir, "/dfs/") 

# # Load custom functions from library
# source("/data/ML_genetics/functions/functions.R")

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
cnp.df     <- cnp[(complete.cases(cnp$CRT_TIME1)),]
swe.df     <- sweden_snp 

# set predictors and targets
predictors <- cnp.df %>% dplyr::select(one_of(snps)) %>% as.matrix() 
targets    <- cnp.df %>% dplyr::select(one_of(outcomes)) %>% apply(2,as.numeric) %>% as.data.frame()

#  apply scaling
means    <- apply(targets,2,mean)
sds      <- apply(targets,2,sd)
targets  <- apply(targets, 2, function(i) (i-mean(i))/sd(i)) %>% as.data.frame

# recode genotypes to dosage of minor allele (continuous)
# 1/1 = 1, 1/2 = 2, 2/2 = 3
# predictors <- apply(predictors,2,function(x) ifelse(x==11,1,ifelse(x==12,2,3)))

## Scale variables manually (so that we can apply these params to replication set)
#  scale based on controls
#  record m/sd
# targets.c <- cnp.df[cnp.df$Group=="control",] %>% dplyr::select(one_of(outcomes)) %>% 
#   apply(2,as.numeric) %>% as.data.frame()
# means     <- apply(targets.c,2,mean)
# sds       <- apply(targets.c,2,sd)
# 
# #  apply scaling
# targets   <- (targets - means)/sds

### Machine Learning preparation
## Set up cross validation
fitControl <- trainControl(method = "repeatedcv",
                           repeats = 3,
                           number = 10,
                           savePredictions = TRUE,
                           selectionFunction = "oneSE")

## Hyperparameter grid search
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

#### Permutation testing for the RF models

rfPipelinePerm <- function(response, Xmat = predictors, grid = rfGrid, cvpar = fitControl,...){
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

bootMe <- function(permutations, target, xmat = predictors,...){
  out_perm <- list()
  for (perm in 1:permutations){
    if(perm %% 5 == 0){print(paste((perm - 1), "perms done so far"))}
    set.seed(perm)
    permutedTarget <- target[shuffle(length(target))]
    mod <- rfPipelinePerm(permutedTarget, Xmat = xmat)[[1]]
    fit <- rfPipelinePerm(permutedTarget, Xmat = xmat)[[2]]
    out_perm[[perm]] <- list(mod,fit)
  }
  return(out_perm)
}



set.seed(1)
 test <- bootMe(2, targets$CVLT_TOTCOR, predictors)


set.seed(1)
perm.loop     <- lapply(select(targets, -CRT_TIME2), function(i) bootMe(10, i, predictors))
save.image(file = "permutations_for_amanda.RData")


model <- test[[1]][[1]]$finalModel
test.rf.pred <- predict(model, newdata = as.matrix(rep.predictors))







#### External Replication
raw_rep <- swe.df[!apply(swe.df[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),] #rm NAs

rep.predictors <- raw_rep %>% dplyr::select(one_of(snps)) %>% as.matrix
rep.predictors <- apply(rep.predictors,2,function(x) ifelse(x==11,1,ifelse(x==12,2,3)))
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


rf.pred.tr        <- lapply(rf.models, function(i) predict(i, newdata=rep.predictors))
glm.pred.tr       <- lapply(glm.models, function(i) predict(i, newdata=rep.predictors))

# out.cors <- list()
# for (i in 1:7){
#   glm.cor <- cor(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
#   rf.cor  <- cor(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")
#   out.cors[[i]] <- c(glm.cor, rf.cor) 
# }
# names(out.cors) <- names(rep.outcomes)
# 
# out.ps <- list()
# for (i in 1:7){
#   glm.cor <- cor.test(glm.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
#   rf.cor  <- cor.test(rf.pred.tr[[i]], rep.outcomes.s[,i], use="pairwise")$p.value
#   out.ps[[i]] <- c(glm.cor, rf.cor)
# }
# names(out.ps) <- names(rep.outcomes)
# 

##### ##### ##### ##### ##### 
##### Mixed effects models controlling for family ID
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

# general ability check 
# calculated from all tasks excepting model-specific task
# scaled, winsorized variables 

# all 6 tasks
g_all <- principal(rep.outcomes.s[,c(1:3,5:7)], nfactors=1, rotation="none",scores=T)
g_all # fit = .97, proportion var = .61
mix$g_all <- as.numeric(g_all$scores)

# leave out trails1
g_tr1 <- principal(rep.outcomes.s[,c(1:2,5:7)], nfactors=1, rotation="none",scores=T)
g_tr1 # fit = .98, proportion var = .72
mix$g_tr1 <- as.numeric(g_tr1$scores)

# leave out vr2
g_vr <- principal(rep.outcomes.s[,c(1:3,5:6)], nfactors=1, rotation="none",scores=T)
g_vr # fit = .96, proportion var = .59
mix$g_vr <- as.numeric(g_vr$scores)

# mixed effect models with g vars
mx.tr1.g <- lme(g_tr1 ~ CRT_TIME1.rf + age + sex, random = ~1 | fid, data = mix, na.action = na.omit)
mx.vr.g <- lme(g_vr ~ VR2DR_TOTALRAW.rf + age + sex, random = ~1 | fid, data = mix, na.action = na.omit)

# to load later
save.image(file=paste0(workDir,"rev_results.Rdata"))

#################### figure
load(paste0(workDir,"rev_results.Rdata"))
library(ggplot2); library(stringr); library(dplyr)
     
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
fig_trails1$Gene2 <- fig_trails1$Gene %>% as.character()
fig_trails1[1,5] <- "C1orf132, CD46*"
fig_trails1[3,5] <- "GLT8D1, GNL3*"
fig_trails1[11,5] <- "AC005609.1, CD14*"
fig_trails1[13,5] <- "C2orf82, EFHD1*"
fig_trails1[17,5] <- "ACD, C16orf86*"
fig_trails1[20,5] <- "CILP2, GATAD2A*"

# visual reproduction - delayed recall
fig_vrdr           <- caret::varImp(rf.models$VR2DR_TOTALRAW, scale=TRUE)[[1]] %>% as.data.frame
fig_vrdr           <- cbind(rownames(fig_vrdr), fig_vrdr) %>% arrange(-Overall)
names(fig_vrdr)    <- c("SNP", "Importance")
fig_vrdr$Gene      <- proxy[match(fig_vrdr$SNP, proxy$proxy), 3]
fig_vrdr$SNP       <- factor(fig_vrdr$SNP,
                             levels = fig_vrdr$SNP[order(fig_vrdr$Importance)])
fig_vrdr$Gene2 <- fig_vrdr$Gene %>% as.character()
fig_vrdr[1,5] <- "C4orf27, CLCN3*"
fig_vrdr[7,5] <- "AC005609.1, CD14*"
fig_vrdr[8,5] <- "GLT8D1, GNL3*"
fig_vrdr[14,5] <- "ADAMTSL3, GOLGA6L4*"
fig_vrdr[15,5] <- "SGSM2, SMG6*"
fig_vrdr[20,5] <- "C2orf82, EFHD1*"
  
# code by shared vs. task-specific genes
fig_trails1$shared         <- ifelse(fig_trails1$Gene %in% fig_vrdr[1:20,3],1,2) %>% as.factor()
fig_vrdr$shared            <- ifelse(fig_vrdr$Gene %in% fig_trails1[1:20,3],1,2) %>% as.factor()
levels(fig_trails1$shared) <- c("both tasks","task-specific")
levels(fig_vrdr$shared)    <- c("both tasks","task-specific")

# trails 1 figure
ggplot(fig_trails1[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.85) +
  geom_text(stat="identity",y=46, hjust=0, aes(label=Gene2), colour="white", size=4) +
  scale_fill_manual(values=c("#B21212","#1485CC")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.06,vjust=1.8,size=18,face="bold")) +
  ggtitle("A) Trails 1/A") +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(45,100))

# vis reprod figure
ggplot(fig_vrdr[1:20,], aes(x=SNP,y=Importance)) + 
  geom_bar(aes(fill=shared),stat="identity",width=.85) +
  geom_text(stat="identity",y=46, hjust=0, aes(label=Gene2), colour="white", size=4) +
  scale_fill_manual(values=c("#B21212","#1485CC")) +  
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=13,color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust=-.06,vjust=1.8,size=18,face="bold")) +
  ggtitle("B) Visual Reproduction II") +
  ylab("\nVariable Importance") +
  xlab("") +
  coord_flip(ylim = c(45,100))


### population stratification within CNP 
cnp_pcs <- read.table("cnp_MLsubs_pca.eigenvec",header=F)
cnp2    <- merge(cnp.df, cnp_pcs[,c(1,3:22)], by.x="ptid",by.y="V1")
cormat  <- cor(cnp2[,c(5:7,9:11,89:108)]) %>% as.data.frame
cormat_edit <- cormat[7:26,1:6]
summary(cormat_edit)