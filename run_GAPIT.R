
source("http://www.bioconductor.org/biocLite.R")
biocLite("multtest", suppressUpdates=TRUE, lib="~/Rlib")

library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

load("../simulated_phenos_genos.rda")

param <- read.table("params.txt", stringsAsFactors=FALSE)
trait <- param[1,1]

Y <- cbind(allPhenos[,c("taxa", trait)])

myGD <- GD
myGM <- GM
myY <- Y

rm(GD, GM, Y); gc()

myK <- read.csv("../GAPIT.Kin.VanRaden.csv", stringsAsFactors=FALSE, header=FALSE)

myY[,trait] <- as.numeric(myY[,trait])
head(myY)

set.seed(345)

foo <- GAPIT(Y=myY,
             GD=myGD,
             GM=myGM,
             KI=myK,
             group.from=9999,
             group.to=9999,
             file.output=FALSE
             )

out <- data.frame(foo$GWAS,
                  SNP.effect=foo$effect.snp,
                  heritability=foo$h2)

gwasFile <- paste0(trait, "_gwas.txt")
write.table(out, gwasFile, row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)


