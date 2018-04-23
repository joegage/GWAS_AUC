library(gplots)
library(ggplot2)
library(cowplot)
library(cobs)
library(ggrepel)

tassels <- read.table("tassel_trait_heritabilities.txt", header=TRUE, stringsAsFactors=FALSE)
rownames(tassels) <- tassels$trait

manual <- tassels[c("TL", "SL", "BN", "weight"), "H2"]
img <- tassels[c("TL_pred", "SL_pred", "BN_pred", "weight_pred"), "H2"]

hDiffs <- data.frame(Trait=c("TL", "SL", "BN", "TW"),
                     diff=abs(manual-img),
                     color=c("dodgerblue3", "darkred", "darkgreen", "darkorchid4"))
hDiffs <- hDiffs[c(2:4, 1), ] ## reorder to make plots clearer

allCompare <- NULL

for(nC in c(10, 100, 1000)){
    datFile <- paste0("ROC_results_", nC, ".rda")
    load(datFile)

    h <- gsub(".*h([0-9]*).*", "\\1", colnames(test))
    r <- gsub(".*r([0-9]*)", "\\1", colnames(test))

    Z <- test
    Z[lower.tri(Z, diag=TRUE)] <- NA
    Z <- -Z
    P <- test
    P[upper.tri(P, diag=TRUE)] <- NA
    P <- -log10(P)

    diff <- matrix(nrow=90, ncol=90)
    h <- as.numeric(h)
    for(i in 1:nrow(diff)){
        for(j in 1:ncol(diff)){
            diff[i,j] <- abs(h[i] - h[j]) / 100
        }
    }

    ## table(diff)
    
    compare <- data.frame(diff=as.numeric(diff), Z=as.numeric(Z), P=as.numeric(P))
    allCompare <- rbind(allCompare, data.frame(compare, NCL=nC))

    ## Plots of Z and P changing with difference in h2
    pPlot <- ggplot(compare, aes(diff, P)) +
        geom_point() +
        ylim(0, 10) +
        stat_summary(fun.y="median", color="red", size=3, geom="point") +
        labs(x="Difference in Heritability", y=expression(-log[10](p))) +
        theme_classic()

    nullDist <- compare[compare$diff == 0 & !is.na(compare$Z),]
    ggplot(nullDist, aes(sample=Z)) +
        geom_qq() +
        geom_abline(slope=1, intercept=0)
    ggsave(paste0("qq_null_", nC, ".png"))
    significant <- quantile(nullDist$Z, c(0.025, 0.975))

    reg <- lm(Z~diff, compare)
    hDiffs$Z <- predict(reg, hDiffs)
    hDiffs$Zhigh <- significant[2]
    hDiffs$Zlow <- significant[1]

    empiricalP <- function(x){
        null <- nullDist$Z
        sum(abs(null) > abs(x)) / length(null)
    }

    hDiffs$pVals <- sapply(1:nrow(hDiffs), function(i) empiricalP(hDiffs$Z[i]) )
    
    zPlot <- ggplot(compare, aes(diff, Z)) +
        geom_hline(yintercept=c(significant), linetype=2, color="gray") +
        geom_point(color="gray", size=0.5) +
        ylim(-2.5, 8) +
        stat_summary(fun.y="median", color="darkgray", size=10, geom="point", pch="_") +
        scale_color_manual(values=c(TL="dodgerblue3", SL="darkred", BN="darkgreen", TW="darkorchid4")) +
        labs(x=ifelse(nC==1000, "Difference in Heritability", ""), y="Z-score") +
        theme_classic()

    write.table(hDiffs[,c("Trait", "Zlow", "Z", "Zhigh", "pVals")],
                paste0("zScores", nC, ".txt"),
                row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

    assign(paste0("z", nC), zPlot)
    assign(paste0("p", nC), pPlot)

    ## Plot averaged ROC curves
    plotAvg <- NULL
    nPoint <- 25
    surface <- NULL
    for(h2loop in (1:9)/10){
        tmp <- avg[avg$h2sim == h2loop, ]

        FPR <- tmp$FPR[seq(0, 1, 1/nPoint)^(1/3)*nrow(tmp)+1]
        FPR <- c(1, FPR, 0)
        TPR <- tmp$TPR[seq(0, 1, 1/nPoint)^(1/3)*nrow(tmp)+1]
        TPR <- c(1, TPR, 0)

        Thresh <- tmp$Threshold[seq(0, 1, 1/nPoint)^(1/3)*nrow(tmp)+1]
        Thresh <- c(0, Thresh, 1)       
        
        plotAvg <- rbind(plotAvg,
                         data.frame(TPR=TPR,
                                    FPR=FPR,
                                    h2=h2loop)
                         )
    }


    cols <- colorRampPalette(c("cadetblue3", "darkblue"))(9)
    if(nC == 1000){
    avgPlot <- ggplot(plotAvg, aes(FPR, TPR, group=h2, color=factor(h2))) +
        geom_smooth(se=FALSE) +
        geom_abline(intercept=0, slope=1, col="gray") +
        xlim(0,1) + ylim(0, 1) +
        scale_color_manual(expression(h^2), values=cols) +
        coord_fixed() +
        theme_classic()
    assign(paste0("avgPlot", nC), avgPlot)
    } else{
        avgPlot <- ggplot(plotAvg, aes(FPR, TPR, group=h2, color=factor(h2))) +
            geom_line() +
            geom_abline(intercept=0, slope=1, col="gray") +
            xlim(0,1) + ylim(0, 1) +
            scale_color_manual(expression(h^2), values=cols, guide=FALSE) +
            coord_fixed() +
            theme_classic()
        assign(paste0("avgPlot", nC), avgPlot)

    }
    
}

together <- plot_grid(z10 + theme(legend.position="none"),
                      z100 + theme(legend.position="none"),
                      z1000 + theme(legend.position="none"),
                      ncol=1,
                      labels=c("a","b","c"))
legendZ <- get_legend(z10 + theme(legend.position="bottom"))
Zplot <- plot_grid(together, legendZ, ncol=1, rel_heights=c(10,1))
save_plot("Zs_and_Ps.png", Zplot, ncol=1, base_height=6, base_width=7)

rocPlots <- plot_grid(avgPlot10, avgPlot100, avgPlot1000 + theme(legend.position="none"), align='vh', hjust=-1, nrow=1, labels=c("a", "b", "c"))
legendR <- get_legend(avgPlot1000)
rocPlotsOut <- plot_grid(rocPlots, legendR, rel_widths=c(3, .3))
save_plot("averaged_ROCs.png", rocPlotsOut, base_width=11.5, base_height=4.5)
rocPlotsVertical <- plot_grid(avgPlot10, avgPlot100, avgPlot1000 + theme(legend.position="none"), align='vh', hjust=-1, nrow=3)
save_plot("vertical_averaged_ROCs.png", rocPlotsVertical, base_width=4.5, base_height=11.5)
