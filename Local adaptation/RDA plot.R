# RDA plot of climate-associated SNPs screened fromLFMM and RDA. take P. roborovskii as an example
setwd("D:/Ro/GEA/RDA_plot/")
.libPaths("D:/R/R-4.3.1/library")
library(vegan)
library("data.table")
library(vcfR)
library(ggrepel)
library(ggplot2)

# Input vcf file of climate-associated SNPs and calculate allele frequency
pr.vcf <- read.vcfR("./input_file/lfmm_rda.recode.vcf")
df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)] # Randomly select one allele
climdata <- read.table("./input_file/clim_data/clim.env", header=T)
sample.coord <- read.table("./input_file/clim_data/sample.coord.txt", header=T, stringsAsFactors=F) #sample name, sample group, and corrdinates
hyp_sam_clim <- cbind(sample.coord[,c(1)], climdata)
colnames(hyp_sam_clim)[1] <- "sample"
identical(rownames(geno), hyp_sam_clim[,1])

rda <- rda(geno ~ bio8 + bio12  + bio15, 
           data = hyp_sam_clim, scale = T)
RsquareAdj(rda)
vif.cca(rda)
anova(rda,by="term")
summary(rda)

# RDA plot
ii <- summary(rda)
sp=as.data.frame(ii$species[,1:2])*5 # Extract the variable coordinates and multiply by 5 for aesthetics
st=as.data.frame(ii$sites[,1:2]) # Extract the sampling coordinates
yz=as.data.frame(ii$biplot[,1:2])*10 # Extract explanatory variable coordinates and multiply by 10 to extend arrows
#se=as.data.frame(ii$constraints[,1:2])

grp <- sample.coord$group
grp <- as.data.frame(grp)
colnames(grp) <- "group"
ggplot() +
  geom_point(data = st,aes(RDA1,RDA2,shape=grp$group,fill=grp$group,color = grp$group),size=7)+ 
  scale_shape_manual(values = c(21,21,21,21,21))+
  scale_color_manual(labels = c("DB","ML","QH","HL","NJ"),
                     breaks = c("DB","ML","QH","HL","NJ"),
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
         fill=guide_legend(title=NULL))+
  theme_bw()+theme(panel.grid=element_blank(),
                   axis.text=element_text(size=16,color="black"),
                   legend.text=element_text(size=16),
                   axis.title.x=element_text(vjust=2, size=16),
                   axis.title.y=element_text(vjust=2, size=16))

ggsave("./output_file/PR_env.pdf", units="in", dpi=400, width=6, height=4.5, device="pdf")

