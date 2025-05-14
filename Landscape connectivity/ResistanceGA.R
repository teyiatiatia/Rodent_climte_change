#R code for assessing landscape connectivity, taking three variables as an example in P. roborovskii
setwd("/Ro/GEA/ResistanceGA/")
.libPaths("D:/R/R-4.3.1/library")

# Calling Julia platform in R
JULIA_HOME <- "D:/software/julia/Julia-1.10.2/bin"
JuliaCall::julia_setup(JULIA_HOME)

install.packages("devtools")
devtools::install_github("wpeterman/ResistanceGA")
library(ResistanceGA)
library(raster)
library(parallel)
library(doParallel)
library(sf)

# Variables are plotted on the same coordinate system with an identical range, using the ASC format. 
elevslop <- raster("D:/inputfile/elev_slop.asc")
ai <- raster("D:/inputfile/ai.asc")
bio12 <- raster("D:/inputfile/bio12.asc")

extent <- c(+70, +130, 25, 55)
elevslop <- crop(elevslop, extent)
ai <- crop(ai, extent)
bio12 <- crop(bio12, extent)
# Stack the variables
all_surface <- stack(elevslop,ai,bio12)

# The coordinate of sample locations and genetic dstance matrix (Fst).
samples<-read.table("D:/inputfile/samples.txt",header=T,row.names = 1, stringsAsFactors=F)
sample.locales <- SpatialPoints(samples[, c(1, 2)])
fst <- read.table("D:/inputfile/popFst.txt", header=T,row.names = 1, stringsAsFactors=F)

# Multithreading
cl <- makePSOCKcluster(6)
registerDoParallel(cl)

# Each layer is analyzed sequentially to generate an optimized resistance layer and estimate its confidence level.
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

# Multiple layers are jointly optimized by evaluating each layer's contribution to generate a combined resistance layer.
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

#optimize single layer or combined layers to generate resistance layer. The combianation range is defined in max.combination. Perform bootstrap analysis to rank the contribution and credibility of each resistance layer.
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

# Save the result
save(SS,file = "D:/outputfile/SS/Results/SS.rda")
save(MS,file = "D:/outputfile/MS/MS.rda")
save(all.comb,file = "D:/outputfile/all_comb/all_comb.rda")
                    
# Conduct additional bootstrap analysis on the combined the SS and MS results.
mat.list <- c(SS$cd,MS$cd)
k <- rbind(SS$k,MS$k)
fst_boots <- fst

colnames(fst_boots) <- NULL
rownames(fst_boots) <- NULL
response <- fst_boots

# Run bootstrap, obs = 18 is the populations
AIC.boots <- Resist.boot(mod.names = names(mat.list),
                         dist.mat = mat.list,
                         n.parameters = k[,2],
                         sample.prop = 0.75,
                         iters = 1000,
                         obs = 18, 
                         genetic.mat = response)
write.csv(AIC.boots, file = paste0(write.dir,"./Results/boots.1000iter.csv"))




