# RDA plot of climate-associated SNPs screened fromLFMM and RDA.
setwd("D:/R/workdir/RDA_GEA_GPA/")
library(vegan)
library("data.table")
library(vcfR)

pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/PR/input_file/clim_pheno_all.recode.vcf") 
###vcf bio8_12 pheno5
pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/LFMM/SNP_anno_struc/SNP_site8_12/bio_pheno_union.recode.vcf") 
###vcf lfrd_bio81215
pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/GF/results_prda_adap/dfmod_raster_81215/input_file/lf_rd_bio81215.recode.vcf")
pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/bio81215/PR/input_file/lf_rd_union.recode.vcf")

df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)] #属于随机提取2个等位基因中的一个，每一个SNP会有Chr01_92104.0/Chr01_92104.1两个互补的等位基因频率，我们只取其中一列，但是79个个体都是同一列。相比于vcftools --012，它统一取的是Chr01_92104.1这一列。只是注意这样不能说0是参考基因组位点，2是变异位点     

climdata <- read.table("D:/R/workdir/workspace/GEA/RDA/clim_data/clim_uncorr_6.env", header=T)
sample.coord <- read.table("D:/R/workdir/workspace/GEA/RDA/clim_data/sample.coord.txt", header=T, stringsAsFactors=F)
hyp_sam_clim <- cbind(sample.coord[,c(1)], climdata)
colnames(hyp_sam_clim)[1] <- "sample"
identical(rownames(geno), hyp_sam_clim[,1])

rda <- rda(geno ~ bio15 + bio8 + bio12  + bio2 + bio11 + bio17, 
           data = hyp_sam_clim, scale = T)

RsquareAdj(rda)
vif.cca(rda)



rda <- rda(geno ~ bio8 + bio12, 
           data = hyp_sam_clim, scale = T)
RsquareAdj(rda)
#$r.squared
#[1] 0.2135446
#$adj.r.squared
#[1] 0.1928484
vif.cca(rda)  #没有共线性
anova(rda,by="term")
#  Model: rda(formula = geno ~ bio8 + bio12, data = hyp_sam_clim, scale = T)
#  Df Variance      F Pr(>F)    
#  bio8      1   1708.3 10.072  0.001 ***
#  bio12     1   1791.7 10.564  0.001 ***
#  Residual 76  12890.0
summary(rda)

#######bio81215
rda <- rda(geno ~ bio8 + bio12 + bio15, 
           data = hyp_sam_clim, scale = T)
RsquareAdj(rda)
vif.cca(rda)
anova(rda,by="term")
summary(rda)

######绘图 8：6画幅 输出pdf "D:\R\workdir\workspace\GEA\GEA_GPA_converg\RDA_GEA_GPA\PR\output_file" 
library(ggrepel)
library(ggplot2)
ii <- summary(rda)
sp=as.data.frame(ii$species[,1:2])*5###提取相应变量坐标，乘以5是使图美观，不影响分析
st=as.data.frame(ii$sites[,1:2])###提取样方坐标，有两种模式，可根据自己数据探索：二选一即可
yz=as.data.frame(ii$biplot[,1:2])*10###提取解释变量坐标,*使箭头变长
#se=as.data.frame(ii$constraints[,1:2])

grp <- c("ML","QH","QH","QH","QH","QH","QH","QH","QH","QH","QH","QH","ML","ML","ML","ML","ML","ML","HL","HL","HL","HL","DB","DB","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","HL","HL","DB","DB","DB","DB","HL","HL","QH","QH","QH","QH","QH","QH","QH","QH","QH","QH","HL","QH","QH","QH","QH","HL","HL","DB","DB","DB","ML","ML","ML","NJ","NJ","NJ","NJ","NJ","HL","HL","HL")
#GEA_GPA_union no HL and QH: grp <- c("ML","ML","ML","ML","ML","ML","ML","DB","DB","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","DB","DB","DB","DB","DB","DB","DB","ML","ML","ML","NJ","NJ","NJ","NJ","NJ")

grp <- as.data.frame(grp)
colnames(grp) <- "group"
ggplot() +
  geom_point(data = st,aes(RDA1,RDA2,shape=grp$group,fill=grp$group,color = grp$group),size=7)+###此处修改圆圈大小
  scale_shape_manual(values = c(21,21,21,21,21))+
  scale_color_manual(labels = c("PRDB","PRML","PRQH","PRHL","PRNJ"),
                     breaks = c("PRDB","PRML","PRQH","PRHL","PRNJ"),
                     values = c("HL"="black", "ML"="black", "QH"="black", "NJ"="black", "DB"="black"),
                     guide = "legend",) +
  scale_fill_manual(values = c("HL"="#31646c", "ML"="#d49c87", "QH"="#96b89b", "NJ"="#20364f", "DB"="#ecd9cf"))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=20,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=1,colour = "68")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)),size = 6)+
  labs(x="RDA1 68.03%",y="RDA2 31.97%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+###此处修改
  theme_bw()+theme(panel.grid=element_blank(),
                   axis.text=element_text(size=16,color="black"),
                   legend.text=element_text(size=16),
                   axis.title.x=element_text(vjust=2, size=16),
                   axis.title.y=element_text(vjust=2, size=16))

ggsave("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/PR/output_file/PR_env.pdf", units="in", dpi=400, width=6, height=4.5, device="pdf")

###bio8_12 pheno5 outfile
#x="RDA1 68.03%",y="RDA2 31.97%"
ggsave("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/PR/output_file_bio8_12_pheno5/PR_env.pdf", units="in", dpi=400, width=6, height=4.5, device="pdf")


####################################################################
##
##pheno
setwd("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/")
library(vegan)
library("data.table")

library(vcfR)
pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/PR/input_file/clim_pheno_all.recode.vcf") #random select 100k SNP from intergenic region 
df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)]
geno <- geno[2:79,]

climdata <- read.table("D:/R/workdir/workspace/GEA/RDA/morph_data/morph.data",header = T)
sample.coord <- read.table("D:/R/workdir/workspace/GEA/RDA/clim_data/sample.coord.txt", header=T, stringsAsFactors=F)
sample.coord <- sample.coord[2:79,]

hyp_sam_clim <- cbind(sample.coord[,c(1)], climdata)
colnames(hyp_sam_clim)[1] <- "sample"
identical(rownames(geno), hyp_sam_clim[,1])

rda_pheno <- rda(geno ~ wei + len + tail + pes  + ear,
                 data = hyp_sam_clim, scale = T)

RsquareAdj(rda_pheno)
vif.cca(rda_pheno)


######绘图
library(ggrepel)
library(ggplot2)
ii <- summary(rda_pheno)
sp=as.data.frame(ii$species[,1:2])*5###提取相应变量坐标，乘以5是使图美观，不影响分析
st=as.data.frame(ii$sites[,1:2])###提取样方坐标，有两种模式，可根据自己数据探索：二选一即可
yz=as.data.frame(ii$biplot[,1:2])*10###提取解释变量坐标,*使箭头变长

grp <- c("QH","QH","QH","QH","QH","QH","QH","QH","QH","QH","QH","ML","ML","ML","ML","ML","ML","HL","HL","HL","HL","DB","DB","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","ML","HL","HL","DB","DB","DB","DB","HL","HL","QH","QH","QH","QH","QH","QH","QH","QH","QH","QH","HL","QH","QH","QH","QH","HL","HL","DB","DB","DB","ML","ML","ML","NJ","NJ","NJ","NJ","NJ","HL","HL","HL")

grp <- as.data.frame(grp)
colnames(grp) <- "group"
ggplot() +
  geom_point(data = st,aes(RDA1,RDA2,shape=grp$group,fill=grp$group,color = grp$group),size=7)+###此处修改
  scale_shape_manual(values = c(21,21,21,21,21))+
  scale_color_manual(labels = c("PRDB","PRML","PRQH","PRHL","PRNJ"),
                     breaks = c("PRDB","PRML","PRQH","PRHL","PRNJ"),
                     values = c("HL"="black", "ML"="black", "QH"="black", "NJ"="black", "DB"="black"),
                     guide = "legend",) +
  scale_fill_manual(values = c("HL"="#31646c", "ML"="#d49c87", "QH"="#96b89b", "NJ"="#20364f", "DB"="#ecd9cf"))+
  geom_segment(data = yz,aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(angle=20,length = unit(0.35,"cm"),
                             type = "closed"),linetype=1, size=1,colour = "68")+
  geom_text_repel(data = yz,aes(RDA1,RDA2,label=row.names(yz)),size = 6)+
  labs(x="RDA1 52.73%",y="RDA2 26.05%")+
  geom_hline(yintercept=0,linetype=3,size=1) + 
  geom_vline(xintercept=0,linetype=3,size=1)+
  guides(shape=guide_legend(title=NULL,color="black"),
         fill=guide_legend(title=NULL))+###此处修改
  theme_bw()+theme(panel.grid=element_blank(),
                   axis.text=element_text(size=16,color="black"),
                   legend.text=element_text(size=16),
                   axis.title.x=element_text(vjust=2, size=16),
                   axis.title.y =element_text(vjust=2, size=16))

###bio8_12_pheno5 outfile
#x="RDA1 44.47%",y="RDA2 31.97%"  pheno输出width修改为6.1
ggsave("D:/R/workdir/workspace/GEA/GEA_GPA_converg/RDA_GEA_GPA/PR/output_file_bio8_12_pheno5/PR_pheno.pdf", units="in", dpi=400, width=6, height=4.5, device="pdf")

