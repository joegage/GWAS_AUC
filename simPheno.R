## The correlation between two traits X and Y is theoretically equal to
## the mean of their heritabilities.

## Use the correlation between Img and Manual as the dependent variable for
## The number of common loci identified?


## h2: scalar between 0 and 1
## nCausative: scalar indicating number of loci that contribute to the trait
## popSize: scalar indicating number of individuals
## alleleFreqs: vector with length nCausative containing allele frequencies
## genotypes: matrix of dimension nCausative x popSize with genotype calls [0,2]

simPheno <- function(h2, nCausative, popSize, alleleFreqs=NULL, genotypes=NULL, nReps){

    effectSize <- rnorm(nCausative, 0, 10)
    
    ## These could be replaced by actual data:
    if(is.null(alleleFreqs) & is.null(genotypes)){
        print("Simulating genotypes")
        alleleFreqs <- runif(nCausative, 0.05, .5)
        genotypes <- rbinom(popSize*nCausative, 1, rep(alleleFreqs, each=popSize))   
        genotypes <- matrix(genotypes, nCausative, byrow=TRUE) *2
    }

    w <- (genotypes - 2*alleleFreqs) / sqrt(2*alleleFreqs*(1-alleleFreqs))
    g <- colSums(w * effectSize)

    varG <- var(g)
    varE <- (varG / h2) - varG

    resid <- rnorm(popSize*length(varE)*nReps, 0, rep(sqrt(varE), each=nReps))
    resid <- matrix(resid, length(varE)*nReps)

    pheno <- t(resid) + g

    return(list(effects=effectSize, g=g, pheno=pheno))

}
