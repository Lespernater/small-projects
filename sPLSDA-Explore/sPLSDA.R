#### Load necessary libraries ####

rm(list = ls()) # Make sure starting with a clean environment
library(knitr)
library(mixOmics)
library(devtools)
library(factoextra)

#### Load breast cancer dataset ####

data("breast.TCGA") 

X <- breast.TCGA$data.train$mirna # use miRNA expression data of 184 miRNA as X matrix
Y <- breast.TCGA$data.train$subtype # use tumour subtype as the Y matrix

head(X[,1:5]) # Columns are features (miRNAs) and rows are samples (individuals)


#### PCA ####

pca.miRNA <- prcomp(X, center = TRUE, scale = TRUE)

PCs = pca.miRNA$x # Gives us 150 PCs

# Scree plot
fviz_eig(pca.miRNA, geom = "bar", main="X Variance Explained by PCA Components", xlab = "Component", ylab = "Variance explained (%)") 

# Visualize clustering 
ggplot(as.data.frame(PCs), aes(y=PC2, x=PC1, colour=Y)) + geom_point(shape=19, size=1) + stat_ellipse() + theme_bw()

# Examine loading vectors
PCAloadings = pca.miRNA$rotation
head(PCAloadings[1:5,1:5]) # Columns are the loading vectors

#### sPLSDA ####

# Initial model
splsda.breast <- splsda(X = X, Y = Y, ncomp = 8) 

# Evaluate model performance
perf.splsda.breast <- perf(splsda.breast, folds = 10, nrepeat = 50, 
                           dist = "mahalanobis.dist") 

plot(perf.splsda.breast, criterion = 'Q2.total') # Plot the perf object

#### Tuning the model ####

list.keepX <- c(seq(5, 120, 5))

# Optimize over possible keepX values
tune.splsda.breast <- tune.splsda(X, Y, ncomp = 8,
                                  test.keepX = list.keepX,
                                  nrepeat = 10, folds = 10)

plot(tune.splsda.breast) # mixOmics makes visualization easy

# Optimizing ncomp and keepX
optimal.keepX <- tune.splsda.breast$choice.keepX 
optimal.ncomp <-  length(optimal.keepX) # extract optimal number of components

final.splsda.breast <- splsda(X, Y, ncomp = optimal.ncomp, 
                              keepX = optimal.keepX) 

perf.final.splsda.breast <- perf(final.splsda.breast, dist = "mahalanobis.dist", 
                                 folds = 10, nrepeat = 50) 

#### Visualizations ####

# Plot loading vectors
plotLoadings(final.splsda.breast, comp = 1, title = "Loadings for Component 1")
plotLoadings(final.splsda.breast, comp = 2, title = "Loadings for Component 2")
plotLoadings(final.splsda.breast, comp = 3, title = "Loadings for Component 3")

# Plot stability of features
plot(perf.final.splsda.breast$features$stable$comp1, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = 'Stability of Comp 1', las =2,
     xlim = c(0, 184),
     ylim = c(0, 1))

plot(perf.final.splsda.breast$features$stable$comp2, type = 'h',
     ylab = 'Stability',
     xlab = 'Features',
     main = 'Stability of Comp 2', las =2,
     xlim = c(0, 184), 
     ylim = c(0, 1))

#### Compare PCA and sPLSDA ####

pca.miRNA = pca(X, center = TRUE, scale = TRUE)

plotIndiv(final.splsda.breast, ind.names = FALSE, 
          rep.space = "X-variate",
          group = breast.TCGA$data.train$subtype, # colour by group
          col.per.group = color.mixo(1:3), 
          legend = TRUE, legend.title = 'Subtype', title = 'sPLSDA comp 1 - 2')

# Let's compare plots
plotIndiv(pca.miRNA, comp = c(1, 2), ind.names = FALSE,
          group = breast.TCGA$data.train$subtype, 
          legend = TRUE, legend.title = 'Subtype', title = 'PCA comp 1 - 2')

## Some more complex visualizations
# col.tox <- color.mixo(as.numeric(as.factor(c(1,3)))) # create set of colours
# library(rgl) # we'll use this to make a 3D plot 
# plotIndiv(final.splsda.breast, ind.names = FALSE, rep.space = "XY-variate", axes.box = "both", style = '3d')

# Some plots that only really appropriate when using multiblock.splsda

#### Prediction using Model ####

predict.model = predict(final.splsda.breast, breast.TCGA$data.test$mirna)

confusion.mat <- get.confusion_matrix(
  truth = breast.TCGA$data.test$subtype, 
  predicted = predict.model$MajorityVote$mahalanobis.dist[,5])
kable(confusion.mat)

#### Multiomic Modelling (DIABLO) ####

X1 <- breast.TCGA$data.train$mirna
X2 <- breast.TCGA$data.train$mrna  
X3 <- breast.TCGA$data.train$protein
X_all <- list(mirna = X1, mrna = X2, protein = X3)
Y_all <- breast.TCGA$data.train$subtype # use tumour subtype as the outcome variable

list.keepX = list(mirna = c(16, 17), mrna = c(18,5), protein = c(5,5)) 
result.sparse.diablo.breast <-  block.splsda(X_all, Y_all, keepX = list.keepX, ncomp = 5) 

# Plot contribution of each modality to each component 1
plotLoadings(result.sparse.diablo.breast, ncomp = 1)

plotIndiv(result.sparse.diablo.breast, var.names = FALSE) # plot the samples

# Plotting
plotVar(result.sparse.diablo.breast, cex = c(2,2,2), var.names = FALSE) # plot the variables

circosPlot(result.sparse.diablo.breast, cutoff = 0.7)

#### Prediction using Model ####
predict.model.diablo = predict(result.sparse.diablo.breast, list(mirna = breast.TCGA$data.test$mirna, mrna =  breast.TCGA$data.test$mrna), protein = NULL) 

# Note protein data is missing for prediction and that is OK. 
# This is explained in warning below

plotDiablo(result.sparse.diablo.breast)

confusion.mat_diablo <- get.confusion_matrix(
  truth = breast.TCGA$data.test$subtype, 
  predicted = predict.model.diablo$WeightedVote$mahalanobis.dist[,5])
kable(confusion.mat_diablo)

# And we didn't even properly tune this model! Imagine the possibilities...
