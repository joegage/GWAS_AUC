source("simPheno.R")

genoFile <- "~/Dropbox/papers/tassel_GWAS_paper/analyses/widiv/genos_phenos_farmCPU.rda"

load(genoFile)
rm(allPhenos); gc()
GM$chrom <- as.character(GM$chrom)
GM$pos <- as.character(GM$pos)

h2 <- seq(0.1, 1, 0.1)
nCausativeValues <- c(10,100, 1000)
nReps <- 10

GDfreqs <- colMeans(GD[,-1])
GDfreqs <- apply(cbind(GDfreqs, 1-GDfreqs), 1, min)
keep <- GDfreqs > 0.02 & GM$chrom %in% 1:10
print(table(keep))
GD <- GD[,c(TRUE, keep)]
GM <- GM[keep, ]


set.seed(345)
for(nCausative in nCausativeValues){
    causativeSNPs <- sample(2:ncol(GD), nCausative)
    genotypes <- t(as.matrix(GD[,causativeSNPs]))
    alleleFreqs <- rowMeans(genotypes)

    simulated <- simPheno(h2, nCausative, nrow(GD), alleleFreqs, genotypes, nReps)
    simulated$effects <- data.frame(SNP=colnames(GD)[causativeSNPs],
                                    chr=GM[causativeSNPs, "chrom"],
                                    pos=GM[causativeSNPs, "pos"],
                                    effect=simulated$effects)
    params <- expand.grid(1:nReps, h2*100)
    colnames(simulated$pheno) <- paste0("c", nCausative, "h", params[,2], "r", params[,1])

    assign(paste0("c", nCausative), simulated)
}

allPhenos <- data.frame(taxa=GD$taxa, c10$pheno, c100$pheno, c1000$pheno)

save("allPhenos", "GD", "GM", file="simulated_phenos_genos.rda")

write.table(allPhenos, "all_simulated_phenos.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)

for(nCausative in nCausativeValues){
    fileName <- paste0("c", nCausative, "_effects.txt")
    outDat <- get(paste0("c", nCausative))$effects
    write.table(outDat, fileName, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}
