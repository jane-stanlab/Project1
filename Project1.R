# Project 1 -- revised

# important -- this makes sure our runs are consistent and reproducible
set.seed(0)  #make your result somewhat predictable

file <- "TCGA_breast_cancer_LumA_vs_Basal_PAM50.tsv"
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold=5
  
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

LumA <- data[data$sample_id %in% first10,header2=='Luminal A']
Basal <- data[data$sample_id %in% first10,header2=='Basal-like']


# define function cross_valid so we can rerun the cross validataion with various parameters
cross_validation <- function (nfold, alg="centroid") {

  # split each cancer type samples into nfold groups
  LumA_groups <- split(sample(colnames(LumA)), 1+(seq_along(colnames(LumA)) %% nfold)) #since the sample function returns random samples, which is not good for debugging, we use seed function at the beggining
  Basal_groups <- split(sample(colnames(Basal)), 1+(seq_along(colnames(Basal)) %% nfold))
  
  result <- array()
 
  # iterate from 1 to nfold groups -- to choose test group
  for (test_group in 1:nfold) {
   
    # return all samples in the chosen test group
    testLumA <- LumA[,colnames(LumA) %in% unlist(LumA_groups[test_group])] #This line gives you the test data for LumA
    testBasal <- Basal[,colnames(Basal) %in% unlist(Basal_groups[test_group])]
   
    # return all samples *not* in the chosen test group 
    trainingLumA <- LumA[,!(colnames(LumA) %in% unlist(LumA_groups[test_group]))]
    trainingBasal <- Basal[,!(colnames(Basal) %in% unlist(Basal_groups[test_group]))]
   
    # compute centroid for each cancer type -- mean for each gene based on all samples
    # note -- rows are gene
    centroidLumA <- rowMeans(trainingLumA) #####glm in logistic regression
    centroidBasal <- rowMeans(trainingBasal)
   
    # For each sample in the test set decide whether it will be classified
    # distance from centroid Lum A: sum(abs(x-centroidLumA))
    # distance from centroid Basal: sum(abs(x-centroidBasal))
    # distance is a sum of distances over all genes 
    # misclassification if when the distance is greater from centroid associated with known result
    misclassifiedLumA <- sum(sapply(testLumA, function(x) { sum(abs(x-centroidLumA))>sum(abs(x-centroidBasal)) })) ####predict in logistic regression
    misclassifiedBasal <- sum(sapply(testBasal, function(x) { sum(abs(x-centroidLumA))<sum(abs(x-centroidBasal)) }))
    
    result[test_group] <- (misclassifiedLumA+misclassifiedBasal)/(ncol(testLumA)+ncol(testBasal))
 }
 
 c(mean(result), sd(result))
}

x<-data.frame(three=cross_validation(3), five=cross_validation(5), ten=cross_validation(10))
# rownames(x) <- c('mean','sd')
print(x) #thw first row is mean, the second row is sd

# #Logistic regression(by Jane) also should do the for loop:see the comments ####for reference
# trainingAB <- rbind(cbind(data.frame(t(trainingLumA)),cancer=0),cbind(data.frame(t(trainingBasal)),cancer=1))
# testLumA0 <- data.frame(t(testLumA))
# model <- glm(cancer~.,data=trainingAB,family=binomial)
# model
# p1<-predict(model, newdata= testLumA0, type="response")
# ifelse(p1<0.5, "0", "1")
# misclassifiedLumA <- sum(ifelse(p1<0.5,0,1)) #labeles as 0, 1:misclassified
# misclassifiedBasal <-sum(ifelse(p2>0.5,0,1)) #basal labeled as one, 0:misclassified
# misclassifiedLumA #michael:0
# testBasal0 <- data.frame(t(testBasal))
# p2<-predict(model, newdata= testBasal0, type="response")
# misclassifiedBasal#michael:2
