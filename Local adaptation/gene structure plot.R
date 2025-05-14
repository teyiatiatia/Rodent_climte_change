#1. vcf to pop freq







setwd("D:/R/workdir/workspace/GEA/GEA_GPA_converg/allele_freq/PR/")
.libPaths("D:/R/R-4.3.1/library")
library(vcfR)
library(adegenet)
library(hierfstat)

clim.points <-  read.table("D:/R/workdir/workspace/GEA/GF/clim.points",sep="\t",header=T)
pr.env<-clim.points
pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/GEA_GPA_converg/allele_freq/PR/gene.recode.vcf")
##bio8_12 foot_tail
pr.vcf <- read.vcfR("D:/R/workdir/workspace/GEA/LFMM/env_K/SNP_anno/SNP_site/bio81215/lf_rd_bio8_2.recode.vcf")


df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)]
all(colnames(pr.vcf@gt)[-1] == pr.env$sample) #检查文件名
pop(df) <- pr.env$pop #set populationp
df$pop
df$other <-pr.env
pop(df) <- pr.env$pop #set populationp
df.genpop <- adegenet::genind2genpop(df)
library(dplyr)
df.genpop$other <- pr.env[,c("pop","longitude", "latitude", "bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19")] |> distinct()
all(rownames(df.genpop$tab) == df.genpop$other$pop)
Xcount <- adegenet::tab(df.genpop, NA.method="zero")
freq_mean <- adegenet::makefreq(df.genpop, missing="mean")  #函数“makefreq”从“genpop”对象中提取具有等位基因频率的表
library(stringr)
freq_Q <- freq_mean[,(!is.na(str_extract(colnames(freq_mean),".0$")))] # 提取Q基因
freq_pr <- freq_Q

write.table(freq_pr, file = "D:/R/workdir/workspace/GEA/GEA_GPA_converg/allele_freq/PR/PR_pop_freq.txt", quote=F, sep = "\t")



####ggplot for gene structure

.libPaths("D:/R/R-4.3.1/library")
setwd("D:/R/workdir/workspace/GEA/LFMM/lfmm_prda_venn/output_set/positions/GEA_GPA_genes/")

library(ggh4x)
library(ggplot2)
library(ggbio)
library(GenomicRanges)


#######################################
setwd("D:/R/workdir/workspace/GEA/GEA_GPA_converg/gene_structure_map/PR/")
####SREBF1_Chr15_55598786
df<-read.table("./input_file/SREBF1_plot.txt",header=F)
df_snp <- read.table("./input_file/SREBF1_plot_snp.txt",header = T)

SREBF1<-GRanges("Chr15",IRanges(df$V4,end=df$V5,group=df$V3))
SREBF1_snp<-data.frame(xmin=unique(df_snp$x_location),
                     xmax=unique(df_snp$x_location),
                     ymin=0.6,
                     ymax=1.4)

pdf(file = "./output_file/SREBF1_Chr15_55598786.pdf",width = 10,height = 4)
autoplot(SREBF1,aes(fill=group),geom="alignment")+
  theme_bw()+
  scale_x_continuous(limits = c(55592900,55599000),
                     breaks = c(seq(55592900,55599000,by=3000)),
                     position = "top")+
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.length = unit(0.2,'cm'))+
  guides(x=guide_axis_truncated(trunc_lower = 55592900,
                                trunc_upper = 55599000))+
  scale_fill_manual(values = c("#92d050","lightgrey"))+
  theme(aspect.ratio = 0.1)+
  scale_y_continuous(limits = c(0,3))+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())+
  annotate(geom = "text",x=55592900,y=1,
           label="\n\nSREBF1")+
  coord_cartesian(clip="off")+
  geom_segment(data=SREBF1_snp,aes(x=xmin,xend=xmax,y=ymin,yend=ymax),
               color="red")
dev.off()
