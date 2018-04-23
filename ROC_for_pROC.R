args <- commandArgs(TRUE)

nC <- args[1]

library(data.table)
library(ggplot2)
library(plyr)
library(pROC)

########################################
## Read files and calculate true/false positive rates


pattern <- paste0("c",nC,"h[0-9]*r[0-9]*_gwas.txt")
files <- list.files("./", pattern=pattern, full.names=TRUE, recursive=TRUE)
files <- files[!grepl("h100", files)] ## Remove sims with h2=100%

toPlot <- NULL
AUC <- NULL
ROC <- list()
i <- 1

for(fileName in files){
    print(paste0("Starting ", fileName))

    nCausative <- as.numeric(gsub(".*/c([0-9]*)h([0-9]*)r([0-9]*)_gwas.txt", "\\1", fileName))
    print(paste0("Causative: ", nCausative))
    h2 <- as.numeric(gsub(".*/c([0-9]*)h([0-9]*)r([0-9]*)_gwas.txt", "\\2", fileName)) / 100
    print(paste0("H^2: ", h2))
    r <- as.numeric(gsub(".*/c([0-9]*)h([0-9]*)r([0-9]*)_gwas.txt", "\\3", fileName))
    print(paste0("Rep: ", r))

    results <- as.data.frame(fread(fileName, header=TRUE, stringsAsFactors=FALSE))
    realFile <- paste0("c", nCausative, "_effects.txt")
    real <- read.table(realFile, header=TRUE, stringsAsFactors=FALSE, sep="\t")

    results$logP <- -log10(results$P.value)
    results$Causative <- results$SNP %in% real$SNP

    thisROC <- roc(results$Causative, results$logP, algorithm=2, quiet=FALSE, levels=c("FALSE", "TRUE"), direction="<")
    thisAUC <- as.numeric(thisROC$auc)

    toPlot <- rbind(toPlot,
                    data.frame(TPR=thisROC$sensitivities,
                               FPR=(1-thisROC$specificities),
                               Threshold=thisROC$thresholds, 
                               threshNum=1:length(thisROC$sensitivities),
                               h2realized=results$heritability[1],
                               h2sim=h2))

    AUC <- c(AUC, thisAUC)
    ROC[[i]] <- thisROC
    i <- i+1

}

ROCplot <- ggplot(toPlot, aes(FPR, TPR, group=h2realized, color=h2sim)) +
    geom_line() +
    geom_abline(slope=1, intercept=0, col="lightgray", linetype=2) +
    xlim(0,1) + ylim(0,1) +
    coord_fixed() +
    labs(main=paste0(nC, " Causative Loci")) +
    theme_classic()
ggsave(paste0("ROC_c", nC, ".png"), ROCplot, width=5, height=5)

summarizeFun <- function(dat){
    TPR <- mean(dat$TPR)
    FPR <- mean(dat$FPR)
    n <- nrow(dat)
    TPR.sd <- sd(dat$TPR)
    FPR.sd <- sd(dat$FPR)
    Thresh <- mean(dat$Threshold)
    Thresh.sd <- sd(dat$Threshold)
    c("TPR"=TPR, "FPR"=FPR, "Threshold"=Thresh, "TPR.sd"=TPR.sd, "FPR.sd"=FPR.sd, "Threshold.sd"=Thresh.sd, "n"=n)
}
avg <- ddply(toPlot, .(h2sim, threshNum), summarizeFun )

test <- matrix(nrow=length(files), ncol=length(files))
allNames <- gsub(".*/(c[0-9]*h[0-9]*r[0-9]*)_gwas.txt", "\\1", files)
colnames(test) <- allNames
rownames(test) <- allNames
for(i in 1:(length(files)-1)){
    for(j in (i+1):length(files)){
        thisTest <- roc.test(ROC[[i]], ROC[[j]], method="delong")
        print(names(thisTest$statistic))
        test[i,j] <- thisTest$statistic
        test[j,i] <- thisTest$p.value
    }
}

assign(paste0("toPlot", nC), toPlot)
assign(paste0("AUC", nC), AUC)
assign(paste0("tests", nC), test)
assign(paste0("files", nC), files)
assign(paste0("avg", nC), avg)


print("Saving...")
toKeep <- c("toPlot", "AUC", "test", "files", "avg")
outFile <- paste0("ROC_results_", nC, ".rda")
save(list=toKeep, file=outFile)

print("Done.")
