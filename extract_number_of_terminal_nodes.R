library(caret); library(randomForest); library(dplyr)


cvlt_treenums <- list()
for (TREENUM in 1:500){
  gt <- getTree(rf.models$CVLT_TOTCOR$finalModel, k = TREENUM)
  cvlt_treenums[[TREENUM]] <- table(gt[,"status"])["-1"]
}
cvlt_treenums <- unlist(cvlt_treenums)
summary(cvlt_treenums)

voc_treenums <- list()
for (TREENUM in 1:500){
  gt <- getTree(rf.models$VOC_TOTALRAW$finalModel, k = TREENUM)
  voc_treenums[[TREENUM]] <- table(gt[,"status"])["-1"]
}
voc_treenums <- unlist(voc_treenums)
summary(voc_treenums)

crt1_treenums <- list()
for (TREENUM in 1:500){
  gt <- getTree(rf.models$CRT_TIME1$finalModel, k = TREENUM)
  crt1_treenums[[TREENUM]] <- table(gt[,"status"])["-1"]
}
crt1_treenums <- unlist(crt1_treenums)
summary(crt1_treenums)

ds_treenums <- list()
for (TREENUM in 1:500){
  gt <- getTree(rf.models$DS_TOTALRAW$finalModel, k = TREENUM)
  ds_treenums[[TREENUM]] <- table(gt[,"status"])["-1"]
}
ds_treenums <- unlist(ds_treenums)
summary(ds_treenums)

vr1_treenums <- list()
for (TREENUM in 1:500){
  gt <- getTree(rf.models$VR1IR_TOTALRAW$finalModel, k = TREENUM)
  vr1_treenums[[TREENUM]] <- table(gt[,"status"])["-1"]
}
vr1_treenums <- unlist(vr1_treenums)
summary(vr1_treenums)

vr2_treenums <- list()
for (TREENUM in 1:500){
  gt <- getTree(rf.models$VR2DR_TOTALRAW$finalModel, k = TREENUM)
  vr2_treenums[[TREENUM]] <- table(gt[,"status"])["-1"]
}
vr2_treenums <- unlist(vr2_treenums)
summary(vr2_treenums)



