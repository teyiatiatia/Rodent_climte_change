# Allele frequency and gene structure of climate-related SNPs
# 1. Allele frequency
setwd("D:/Ro/GEA/allele_freq/")
.libPaths("D:/R/R-4.3.1/library")
library(vcfR)
library(adegenet)
library(hierfstat)
library(dplyr)
library(xlsx)
library(ggplot2)
library(ggpubr)

clim.points <- read.table("./input_file/clim.points",sep="\t",header=T)
pr.env <- clim.points
pr.vcf <- read.vcfR("./input_file/lfmm_rda.recode.vcf")

df <-vcfR2genind(pr.vcf) #vcf2genind
geno <- df$tab
geno <- geno[,seq(1,ncol(geno)-1,2)]
all(colnames(pr.vcf@gt)[-1] == pr.env$sample) 
pop(df) <- pr.env$pop #set populationp
df$other <-pr.env
pop(df) <- pr.env$pop #set populationp
df.genpop <- adegenet::genind2genpop(df)
df.genpop$other <- pr.env[,c("pop","longitude", "latitude")] |> distinct()
all(rownames(df.genpop$tab) == df.genpop$other$pop)
Xcount <- adegenet::tab(df.genpop, NA.method="zero")
freq_mean <- adegenet::makefreq(df.genpop, missing="mean")  # Extract the allele frequency table
library(stringr)
freq_Q <- freq_mean[,(!is.na(str_extract(colnames(freq_mean),".0$")))] 
freq_pr <- freq_Q

write.table(freq_pr, file = "./output_file/PR_pop_freq.txt", quote=F, sep = "\t")

# Allele frequency table needs to have additional climate variables and population names included.
# Allele frequency plot for Chr15_55598786
data_1<-read.xlsx2("./output_file/PR_pop_freq.xlsx",
                   sheetName = "Sheet1",
                   header = TRUE,
                   as.data.frame = TRUE,
                   row.names = 1)
data_1$bio8<-as.numeric(data_1$bio8)
data_1$Chr15_55598786.1<-as.numeric(data_1$Chr15_55598786.1) 
ggplot(data_1,aes(x=bio8,y=Chr15_55598786.1))+
  ylab("Allele Frequency")+
  geom_point(size=5,alpha=0.7,color="black")+
  geom_smooth(method = "lm", formula = y~x, color = "black", fill = "gray")+
  stat_cor(method = "pearson", label.x = 16, label.y = 1.2) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text=element_text(size=18,color="black"),
    legend.text=element_text(size=18),
    axis.title.x=element_text(vjust=2, size=18),
    axis.title.y =element_text(vjust=2, size=18)
  )+geom_point(data = data_1, aes(x = bio8, y =Chr15_55598786.1,color=pop.1),size=5)+theme(legend.position = 'none')
ggsave("./output_file/bio8_Chr15_55598786.1_SREBF1.pdf",width = 5, height = 5,units = "in") 

# 2. ggplot for gene structure of gene SREBF1
.libPaths("D:/R/R-4.3.1/library")
setwd("D:/Ro/GEA/allele_freq/gene_struc/")
library(ggh4x)
library(ggplot2)
library(ggbio)
library(GenomicRanges)

df<-read.table("./input_file/SREBF1_plot.txt",header=F) # This gene annotation information
df_snp <- read.table("./input_file/SREBF1_plot_snp.txt",header = T) # The annotation information of this SNP

SREBF1 <- GRanges("Chr15",IRanges(df$V4,end=df$V5,group=df$V3))
SREBF1_snp <- data.frame(xmin=unique(df_snp$x_location),
                     xmax=unique(df_snp$x_location),
                     ymin=0.6,
                     ymax=1.4)

pdf(file = "./output_file/SREBF1_Chr15_55598786.pdf",width = 10,height = 4)
autoplot(SREBF1,aes(fill = group),geom="alignment")+
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
