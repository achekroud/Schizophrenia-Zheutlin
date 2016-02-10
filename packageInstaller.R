# Adam's package installer

localLib <- "/home/azheutlin/Rlibs"

# 
# "dplyr", "ggplot2", "RColorBrewer", "ape", "dynamicTreeCut","cluster",
# "nFactors","clustsig", "MASS", "glmnet", "caret", "pROC", "permute", "gbm",
# "klaR", "lme4", "nlme","gee", "plyr", "dplyr", "foreach", "doMC")

basics   <- c("dplyr", "ggplot2", "RColorBrewer")
parallel <- c("foreach", "doMC")


wants <- c("dplyr", "ggplot2", "RColorBrewer", "ape", "dynamicTreeCut","cluster",
           "nFactors","clustsig", "glmnet", "caret", "pROC", "permute", "gbm",
           "klaR", "lme4", "nlme","gee", "plyr", "dplyr", "foreach", "doMC")
has   <- wants %in% rownames(installed.packages())

# if(any(!has)) install.packages(wants[!has], lib.loc = localLib)
if(any(!has)) install.packages(wants[!has])


basics   <- c("dplyr", "ggplot2", "RColorBrewer")
parallel <- c("foreach", "doMC")
mixed    <- c("lme4", "nlme", "gee")

loadBasics <- function(pac = basics){
  library(pac)
}

loadParallel <- function(pac = parallel){
  library(pac)
}

loadMixed <- function(pac = mixed){
  library(pac)
}

loadAll <- function(pac = wants){
  library(pac)
}


