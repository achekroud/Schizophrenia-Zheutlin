# let's try to make plink be R

workdir <- "/data/swe_gwas/ABZ/ML/"
genedir <- "/data/ML_genetics/Schizophrenia-Zheutlin/dfs/"

# select only subjects used in project
swedes <- read.table(paste0(genedir,"GWASinput2-swe.txt"),header=T)
cnp <- read.table(paste0(genedir,"GWASinput2-cnp.txt"),header=T)
swedes <- swedes[,c(1:11,13:92)]
cnp <- cnp[,c(1:8,10:89)]

cnp.df <- cnp[(complete.cases(cnp$CRT_TIME1)),]
outcomes <- names(cnp.df)[5:11]
names(swedes)[8:14] <- outcomes 
swe.df <- swedes[!apply(swedes[,c(outcomes)], 1, function(x) (mean(is.na(x))==1)),]

write.table(cnp.df[,c(1,1)],paste0(workdir,"cnp_subs.txt"),col.names=F,
            row.names=F,quote=F,sep="\t")
write.table(swe.df[,c(2,1)],paste0(workdir,"swe_subs.txt"),col.names=F,
            row.names=F,quote=F,sep="\t")

# map population stratification: cnp vs. swe
mds <- read.table(paste0(workdir,"all_MLsubs2.mds"),header=T)
cnpfam <- read.table(paste0(workdir,"cnp_callgenos_merged.fam"),header=F)
names(cnpfam)[1] <- "IID"
mds$sample <- ifelse(mds$IID %in% cnpfam$IID,"USA","Sweden") %>% as.factor()

ggplot(mds,aes(x=C1,y=C2)) +
  geom_point(aes(colour=sample)) + 
  scale_colour_manual(values=c("#0000FF","#CC0000")) +
  theme(axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=10,color="black"),
        panel.background=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(color = "black"))