library(ggplot2)
library(cowplot)

phenos <- read.table("~/Dropbox/papers/tassel_GWAS_paper/analyses/widiv/allTasselTraits_2013_through_2015_BLUPs.txt", header=TRUE, stringsAsFactors=FALSE, sep="\t")

phenos$TL_pred <- phenos$TL_pred / 10
phenos$SL_pred <- phenos$SL_pred/ 10

cors <- cor(phenos[,-1], use="complete.obs")

TLcor <- paste0(" r = ", round(cors["TL", "TL_pred"], 2))
TL <- ggplot(phenos, aes(TL, TL_pred)) +
    geom_point(color="dodgerblue3", alpha=0.5) +
    geom_text(aes(x=-Inf, y=Inf), label=TLcor, hjust=0, vjust=1) +
    labs(x="TL", y="TLp") +
    theme_classic()
TL

SLcor <- paste0(" r = ", round(cors["SL", "SL_pred"], 2))
SL <- ggplot(phenos, aes(SL, SL_pred)) +
    geom_text(aes(x=-Inf, y=Inf), label=SLcor, hjust=0, vjust=1) +
    geom_point(color="darkred", alpha=0.5) +
    labs(x="SL", y="SLp") +
    theme_classic()
SL

BNcor <- paste0(" r = ", round(cors["BN", "BN_pred"], 2))
BN <- ggplot(phenos, aes(BN, BN_pred)) +
    geom_text(aes(x=-Inf, y=Inf), label=BNcor, hjust=0, vjust=1) +
    geom_point(color="darkgreen", alpha=0.5) +
    labs(x="BN", y="BNp") +
    theme_classic()
BN

TWcor <- paste0(" r = ", round(cors["weight", "weight_pred"], 2))
TW <- ggplot(phenos, aes(weight, weight_pred)) +
    geom_point(color="darkorchid4", alpha=0.5) +
    geom_text(aes(x=-Inf, y=Inf), label=TWcor, hjust=0, vjust=1) +
    labs(x="TW", y="TWp") +
    theme_classic()
TW

all <- plot_grid(TL, SL, BN, TW, ncol=2, labels=c("a", "b", "c", "d"))
save_plot("fig1_correlations.tiff", all, ncol=2, base_width=3, base_height=6)
save_plot("fig1_correlations.pdf", all, ncol=2, base_width=3, base_height=6)
