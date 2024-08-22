### This script is to do the genomic selection
library(rrBLUP)
library(cvTools)
## Convert the vcf data to the {-1,0,1} format which made by plink with ped and map
# Denote the function
parse.GBS <- function(x) {
  unique.x <- unique(x)
  alleles <- c("0/0","1/1")
  y <- rep(0,length(x))
  y[which(x==alleles[1])] <- -1
  y[which(x==alleles[2])] <- 1
  y[which(x=="./.")] <- NA
  return(y)
}
# read the vcf data 
raw_vcf <- read.table("3K.BJ.SLB.GT.vcf", header = TRUE, sep = '\t',stringsAsFactors = F,comment.char = "")
# convert the vcf to {-1,0,1}
marker <- apply(raw_vcf[,-c(1:9)], 1, parse.GBS)
# impute the missing
marker_impute <- A.mat(marker, impute.method = "mean", max.missing = 1, return.imputed = T)
#marker_impute <- A.mat(marker, impute.method = "EM", max.missing = 1, return.imputed = T)
myGD1 <- marker_impute$imputed
# change name 
dimnames(myGD1)[[1]] <- colnames(raw_vcf)[10:ncol(raw_vcf)]
dimnames(myGD1)[[2]] <- raw_vcf$ID
#write.table(myGD1, "3K.BJ.SLB.GT.numeric.txt", sep = '\t', quote = F, row.names = T, col.names = NA)

## begin to do GS (mix.solve)
# read data
pheno <- read.table("1K.NE.NLB.phenotype.txt",sep = '\t',header = T, stringsAsFactors = F)
geno_ms <- myGD1 # genotype_ms
geno_kin <- A.mat(myGD1) # genotype_kin
# set the relative parameters
k <- 5 # fold
l <- 20 # repeat time for one pop size
pop <- c(200,400,600,800,1000,nrow(pheno)) # pop=200,400,600,800,1000,1121 do cycle
# set the output file (two model)
corr_ms <- matrix(NA, nrow = k*l, ncol = length(pop), dimnames = list(c(1:(k*l)), pop))
R_sqr_ms <- matrix(NA, nrow = k*l, ncol = length(pop), dimnames = list(c(1:(k*l)), pop))
corr_kin <- matrix(NA, nrow = k*l, ncol = length(pop), dimnames = list(c(1:(k*l)), pop))
R_sqr_kin <- matrix(NA, nrow = k*l, ncol = length(pop), dimnames = list(c(1:(k*l)), pop))
# begin the prediction
pop.idx <- 0 # set the pop.idx for the column of output
for (p in pop){
  # cycle for pop size
  set.seed(1234+p) # fix the data
  fold <- cvFolds(p, K = 5) # 5 fold and no repeat
  ix <- 0 # ix represent the row of output
  pop.idx <- pop.idx + 1
  for (i in 1:l){
    # cycle for repeat
    set.seed(108+i) # fix the data
    pop_sample <- sample(c(1:nrow(pheno)), size = p, replace = FALSE) # set the training sample
    pop_geno_ms <- geno_ms[pop_sample,] # extract the genotype of ms for pop
    pop_geno_kin <- geno_kin[pop_sample, pop_sample] # extract the genotype of kin for pop
    pop_pheno <- pheno[pop_sample,] # extract the phenotype for pop
    for (j in 1:k){
      # cycle for fold
      ix <- ix + 1 # cycle once, ix +1 ,each row of the corr represent one prediction
      test_sample <- fold[[4]][fold[[5]]==j,1] # extract the testing set
      ### mix.solved ###
      train_geno_ms <- pop_geno_ms[-test_sample,] # genotype of training set of ms
      train_pheno_ms <- pop_pheno[-test_sample] # phenotype of training set of ms
      test_geno_ms <- pop_geno_ms[test_sample,] # genotype of testing set of ms
      test_pheno_ms <- pop_pheno[test_sample] # phenotype of testing set of ms
      model <- mixed.solve(train_pheno_ms, Z = train_geno_ms, K = NULL, SE = FALSE, return.Hinv = FALSE) # build model
      GEBV <- as.matrix(test_geno_ms) %*% model$u # calculate the GEBV
      y.est_ms <- GEBV[,1] + model$beta # estimated phenotype
      y.obs_ms <- test_pheno_ms # observed phenotype
      corr_ms[ix,pop.idx] <- cor(y.est_ms, y.obs_ms) # calculate the corr
      R_sqr_ms[ix,pop.idx] <- 1-(var(y.obs_ms-y.est_ms)/var(y.obs_ms)) # calculate the R sqrt
      ### kin.blup ###
      yNA <- pop_pheno # set the phenotype
      yNA[test_sample] <- NA # convert the test sample phenotype to NA
      data <- data.frame(y=yNA, gid=rownames(pop_geno_kin)) # create a dataframe
      ans <- kin.blup(data = data, geno = "gid", pheno = "y", K = pop_geno_kin) # kin.blup prediction
      y.est_kin <- ans$g[test_sample] # estimated phenotype
      y.obs_kin <- pop_pheno[test_sample] # observed phenotype
      corr_kin[ix,pop.idx] <- cor(y.est_kin, y.obs_kin) # calculate the corr
      R_sqr_kin[ix,pop.idx] <- 1-(var(y.obs_kin-y.est_kin)/var(y.obs_kin)) # calculate the R sqrt
    }
  }
}
# write the output
write.table(corr_ms, "1K.NE.NLB.corr.ms.txt", sep = '\t', quote = F, row.names = T)
write.table(R_sqr_ms, "1K.NE.NLB.rsqr.ms.txt", sep = '\t', quote = F, row.names = T)
write.table(corr_kin, "1K.NE.NLB.corr.kin.txt", sep = '\t', quote = F, row.names = T)
write.table(R_sqr_kin, "1K.NE.NLB.rsqr.kin.txt", sep = '\t', quote = F, row.names = T)