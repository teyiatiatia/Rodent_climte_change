# Windows
julia
using Pkg
using RCall  #download RCall
Pkg.build("RCall")
Pkg.add("Circuitscape")

# R
setwd("/Ro/GEA/ResistanceGA/")
.libPaths("D:/R/R-4.3.1/library")

install.packages("devtools")
devtools::install_github("wpeterman/ResistanceGA")

# Calling Julia platform in R
JULIA_HOME <- "D:/software/julia/Julia-1.10.2/bin"
JuliaCall::julia_setup(JULIA_HOME)

library(ResistanceGA)
library(raster)
library(parallel)
library(doParallel)
library(sf)

# Take three variables as an example. Variables are plotted on the same coordinate system with the same range in ASC format.
elevslop <- raster("D:/inputfile/elev_slop.asc")
ai <- raster("D:/inputfile/ai.asc")
bio12 <- raster("D:/inputfile/bio12.asc")
#图层可以裁剪到固定区域
extent <- c(+70, +130, 25, 55)
elevslop <- crop(elevslop, extent)
ai <- crop(ai, extent)
bio12 <- crop(bio12, extent)

#将图层变量叠加在一起
all_surface <- stack(elevslop,ai,bio12)

#载入采样点坐标信息和样点之间的遗传距离，遗传距离是样点之间的fst值矩阵，可以在hierfstat R包中利用vcf文件生成。
samples<-read.table("D:/inputfile/samples.txt",header=T,row.names = 1, stringsAsFactors=F)
sample.locales <- SpatialPoints(samples[, c(1, 2)])
fst <- read.table("D:/inputfile/popFst.txt", header=T,row.names = 1, stringsAsFactors=F)

#多线程
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

#单个图层依次分析，会对每一个图层优化生成阻碍值图层，并且预估每个图层的置信度
write.dir <- "D:/outputfile/SS/"
jl.inputs <- jl.prep(n.Pops = length(sample.locales),
                     response = fst[lower.tri(fst)],
                     CS_Point.File = sample.locales,
                     JULIA_HOME = JULIA_HOME,
                     run_test = F)
GA.inputs <- GA.prep(ASCII.dir = all_surface,
                     Results.dir = write.dir,
                     min.cat = 0,
                     max.cat = 1000,
                     max.cont = 1000,
                     method = "AIC",
                     parallel=cl,
                     quiet = TRUE)
SS <- SS_optim(jl.inputs = jl.inputs,
               GA.inputs = GA.inputs)

#多个图层联合优化，会评估每个图层的贡献值，生成联合的阻碍值图层
write.dir <- "D:/outputfile/MS/"
jl.inputs <- jl.prep(n.Pops = length(sample.locales),
                     response = fst[lower.tri(fst)],
                     CS_Point.File = sample.locales,
                     JULIA_HOME = JULIA_HOME,
                     run_test = F)
GA.inputs <- GA.prep(ASCII.dir = all_surface,
                     Results.dir = write.dir,
                     min.cat = 0,
                     max.cat = 1000,
                     max.cont = 1000,
                     method = "AIC",
                     parallel=cl,
                     quiet = TRUE)
MS <- MS_optim(jl.inputs = jl.inputs,
               GA.inputs = GA.inputs)

#自由组合单个或多个图层，生成单个或者多个图层的阻碍值图层。自由组合的范围在 max.combination 中设置，max.combination = c(2,3)指的是自由组合环境图层中的2个或者3个图层联合优化。也可以设置为max.combination = 3，意味着将进行单个图层优化，2个图层自由组合联合优化，3个图层自由组合联合优化
#同时进行bootstrap分析，会给每一个单个图层生成的阻碍图层或多个图层生成的联合阻碍图层的贡献值和可信度排名
write.dir <- "D:/outputfile/all_comb/"
CS.inputs <- jl.prep(n.Pops = length(sample.locales),
                     response = fst[lower.tri(fst)],
                     CS_Point.File = sample.locales,
                     JULIA_HOME = JULIA_HOME,
                     run_test = F)

GA.inputs <- GA.prep(ASCII.dir = all_surface,
                     Results.dir = 'all.comb',
                     min.cat = 0,
                     max.cat = 1000,
                     max.cont = 1000,
                     method = "AIC",
                     parallel=cl,
                     quiet = TRUE)
all_comb <- all_comb(jl.inputs = CS.inputs,
                     GA.inputs = GA.inputs, 
                     results.dir = write.dir,
                     max.combination = c(2,3),
                     iters = 1000,
                     replicate = 1,
                     sample.prop = 0.75,
                     nlm = FALSE,
                     dist_mod = TRUE,
                     null_mod = TRUE)

#如果需要保存结果，可以先写出
save(SS,file= "D:/outputfile/SS/Results/SS.rda")
#再读入结果
SS<-load("D:/outputfile/SS/Results/SS.rda ")                      
#如果想单独对单个表面运算结果以及多个表面联合生成的结果进行额外的bootstrap分析
mat.list <- c(SS$cd,MS$cd)
k <- rbind(SS$k,MS$k)
fst_boots <- fst

colnames(fst_boots) <- NULL
rownames(fst_boots) <- NULL
response <- fst_boots

# Run bootstrap，obs = 18为样点数
AIC.boots <- Resist.boot(mod.names = names(mat.list),
                         dist.mat = mat.list,
                         n.parameters = k[,2],
                         sample.prop = 0.75,
                         iters = 1000,
                         obs = 18, 
                         genetic.mat = response)
write.csv(AIC.boots, file = paste0(write.dir,"./Results/boots.1000iter.csv"))




