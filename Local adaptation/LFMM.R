#LFMM
setwd("/Ro/GEA/LFMM")
.libPaths("D:/R/R-4.3.1/library")
library(LEA)
library(qvalue)
library(qqman)
library(qvalue)

project = NULL
clim.env <- read.table("./clim.points",header = T)
project = lfmm("./SNP_file/pr.lfmm", "./clim.env", K = 5, repetitions = 5, CPU = 10, iterations = 1000, burnin = 500, project = "new") #重复5次

#Merge 5 repitations and calculate the adjust p-value for each variable.
z.bio8 = z.scores(project, K = 5, d = 1)
z.bio8 <- apply(z.bio8, 1, median)
lambda.bio8 = median(z.bio8^2)/qchisq(0.5, df = 2)
p.bio8.adj = pchisq(z.bio8^2/lambda.bio8, df = 1, lower = FALSE)

z.bio12 = z.scores(project, K = 5, d = 3)
z.bio12 <- apply(z.bio12, 1, median)
lambda.bio12 = median(z.bio12^2)/qchisq(0.5, df = 1)
p.bio12.adj = pchisq(z.bio12^2/lambda.bio12, df = 1, lower = FALSE)

z.bio15 = z.scores(project, K = 5, d = 3)
z.bio15 <- apply(z.bio15, 1, median)
lambda.bio15 = median(z.bio15^2)/qchisq(0.5, df = 1)
p.bio15.adj = pchisq(z.bio15^2/lambda.bio15, df = 1, lower = FALSE)

q.bio8<-qvalue(p.bio8.adj)$qvalues
q.bio12<-qvalue(p.bio12.adj)$qvalues
q.bio15<-qvalue(p.bio15.adj)$qvalues

#Take bio8 as an example, check the QQ plot, histogram and Manhattan plot.
sum(q.bio8<0.05 & abs(z.bio8)>2)
qq(p.bio8.adj)
hist(p.bio8.adj, col = "lightblue",main = "p.bio8.adj")
plot(-log10(q.bio8), pch = 19, col = "grey", cex = .2, xlab = '')

lfmm.results <- cbind(z.bio8, q.bio8, z.bio12, q.bio12, z.bio15, q.bio15)
write.table(lfmm.results, "./results/lfmm.results", sep="\t", quote=F, row.names=F)

#output outliers
alpha <- 0.05
outliers_8 <- which(q.bio8 < alpha)
write.table(outliers_8,file ="./results/bio8_outliers",row.names = FALSE)
outliers_12 <- which(q.bio12 < alpha)
write.table(outliers_12,file ="./results/bio12_outliers",row.names = FALSE)
outliers_15 <- which(q.bio15< alpha)
write.table(outliers_15,file ="./results/bio15_outliers",row.names = FALSE)























