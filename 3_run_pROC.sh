for N in 10 100 1000; do
    nohup Rscript ROC_for_pROC.R $N &> ROC_for_pROC_n$N.log &
done
