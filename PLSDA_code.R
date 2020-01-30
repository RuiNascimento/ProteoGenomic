library(mixOmics)
data(srbct)
X <- srbct$gene
Y <- srbct$class 
summary(Y) ## class summary


MyResult.splsda <- splsda(X, Y, keepX = c(50,50)) # 1 Run the method
plotIndiv(MyResult.splsda)                        # 2 Plot the samples (coloured by classes automatically)
plotVar(MyResult.splsda)                          # 3 Plot the variables




var <- apply(SummarizedExperiment::assay(dep), 1, sd)
df <- SummarizedExperiment::assay(dep)[order(var, decreasing = TRUE),]

df6 <- df[,c(7,8,9,16,17,18)]
df12 <- df[,c(1,2,3,10,11,12)]
df24 <- df[,c(4,5,6,13,14,15)]

MyResult.splsda.6h <- splsda(t(df6), Y = c("I", "I", "I", "M", "M", "M"))
MyResult.splsda.12h <- splsda(t(df12), Y = c("I", "I", "I", "M", "M", "M"))
MyResult.splsda.24h <- splsda(t(df24), Y = c("I", "I", "I", "M", "M", "M"))

plotIndiv(MyResult.splsda.6h, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title = 'sPLS-DA at 6hpi',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

plotIndiv(MyResult.splsda.12h, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title = 'sPLS-DA at 12hpi',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

plotIndiv(MyResult.splsda.24h, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title = 'sPLS-DA at 24hpi',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

MyResult.splsda.all <- splsda(t(df), Y = c("12I", "12I", "12I", "24I", "24I", "24I", "6I", "6I", "6I", "12M", "12M", "12M", "24M", "24M", "24M", "6M", "6M", "6M"))
plotIndiv(MyResult.splsda.all, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title = 'sPLS-DA of all Samples',
          X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
