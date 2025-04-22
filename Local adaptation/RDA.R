setwd("/Ro/GEA/pRDA/")
.libPaths("D:/R/R-4.3.1/library")
library(vcfR)
library(adegenet)
library(hierfstat)
library(stringr)
library(dplyr)
library(rdacca.hp)
library(vegan)
library("data.table")
library("dplyr")

clim.points <- read.table("./input_file/clim.points",sep="\t",header=T) #Includes grouping information for corresponding individual groups.  
pr.env<-clim.points

# PCA on intergenic region
pr.vcf <- read.vcfR("./input_file/PR.intergenic.vcf") 
df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)] 
all(colnames(pr.vcf@gt)[-1] == pr.env$sample)
pop(df) <- pr.env$pop #set populationp
df$other <-pr.env
df.genpop <- adegenet::genind2genpop(df)
df.genpop$other <- pr.env[,c("pop","longitude", "latitude")] |> distinct()
all(rownames(df.genpop$tab) == df.genpop$other$pop)
Xcount <- adegenet::tab(df.genpop, NA.method="zero")
freq_mean <- adegenet::makefreq(df.genpop, missing="mean") 
freq_Q <- freq_mean[,(!is.na(str_extract(colnames(freq_mean),".0$")))] 
freq_pr_intergenic <- freq_Q

# Running a PCA on neutral genetic markers
pca <- rda(freq_pr_intergenic, scale=T)

# Screeplot of the PCA eigenvalues
screeplot(pca, type = "barplot", npcs=10, main="PCA Eigenvalues")

# Neutral population structure table
PCs <- scores(pca, choices=c(1:3), display="sites", scaling=0)
PopStruct <- data.frame(Population = row.names(freq_pr), PCs)
colnames(PopStruct)

clim.pop.points <- read.table("./input_file/clim.pop.points",sep= "\t",header=T)
geo <- scale(clim.pop.points[,c(1,2)]) #Coordinates of the sampling locations
Variables <- data.frame(PopStruct, clim.pop.points[,c(8,12,15)],geo)

# pop freq for total SNP
pr.vcf <- read.vcfR("./input_file/PR.vcf")
df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)] 
all(colnames(pr.vcf@gt)[-1] == pr.env$sample) 
pop(df) <- pr.env$pop #set populationp
df$other <-pr.env
df.genpop <- adegenet::genind2genpop(df)
library(dplyr)
df.genpop$other <- pr.env[,c("pop","longitude", "latitude")] |> distinct()
all(rownames(df.genpop$tab) == df.genpop$other$pop)
Xcount <- adegenet::tab(df.genpop, NA.method="zero")
freq_mean <- adegenet::makefreq(df.genpop, missing="mean")  
library(stringr)
freq_Q <- freq_mean[,(!is.na(str_extract(colnames(freq_mean),".0$")))] 
freq_pr <- freq_Q

# RDA-Variance Partitioning
# Repeat variance partitioning analysis 1000 times and calculate the proportion of the neutral genetic structure, the proportion explained by climatic variables, the proportion jointly explained by the former two, and the proportion unexplained.
traits=4 
cycles=1000
Adj_Rs = matrix(nrow=cycles, ncol=traits)
for(r in 1:cycles){
  train= as.matrix(sample(1:205149, 20000))
  rda_vp <- varpart(freq_pr[,train], Variables[,c(2:4)], Variables[,c(7:8)])
  Adj_Rs[r,c(1:4)]<-rda_vp$part$indfract$Adj.R.squared
}
colMeans(Adj_Rs)
apply(Adj_Rs,2,sd) #sd

#pRDA to identify candidate adaptive SNPs
pRDA <- rda(freq_pr ~ bio15 + bio8 + bio12 + 
                 Condition(PC1 + PC2 + PC3),  Variables)
RsquareAdj(pRDA)
anova(pRDA,by="term")
screeplot(pRDA, main="Eigenvalues of constrained axes")
load.prda <- scores(pRDA, choices=c(1:3), display="species")
write.table(load.prda,file="./output_file/load_dbrda_3.txt", quote=F)

SNP_df <- as.data.frame(row.names(load.prda))
SNP_df$pos <- 1:nrow(SNP_df)
names(SNP_df) <- c("snp", "pos")

#standard deviation cutoff
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.prda[,1],3)    # if a conservative approach is required, set it to 3.5 to minimize false positives. 
cand2 <- outliers(load.prda[,2],3)
cand3 <- outliers(load.prda[,3],3)

ncand <- length(cand1) + length(cand2) + length(cand3) 
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)    # correlations of each candidate SNP with the predictors
foo <- matrix(nrow=nrow(cand), ncol=3)
colnames(foo) <- c("bio8","bio12","bio15")

for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- freq_pr[,nam]
  foo[i,] <- apply(Variables[,c(7,8,9)],2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)
length(cand$snp[duplicated(cand$snp)])
cand <- cand[!duplicated(cand$snp),]     # remove duplicate detections
head(cand)
 
cand.pos <- left_join(cand, SNP_df, by = "snp")

for (i in 1:length(cand.pos$snp)) {
  bar <- cand.pos[i,]
  cand.pos[i,8] <- names(which.max(abs(bar[4:6]))) # gives the predictor 
  cand.pos[i,9] <- max(abs(bar[4:6]))              # gives the correlation value
}
colnames(cand.pos)[8] <- "predictor"
colnames(cand.pos)[9] <- "correlation"

write.table(cand.pos[,7], "./output_file/cand.pos.prda.sites", quote = F, row.names = F,col.names = T)
write.table(cand.pos, "./output_file/cand.pos.prda.3axis", quote = F, row.names = F,col.names = T)
