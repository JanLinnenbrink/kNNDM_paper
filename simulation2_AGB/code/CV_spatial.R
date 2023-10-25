# *****************************************************************************
# R Script implementing spatial  k-fold cross-validation.  
# Modified from https://doi.org/10.5281/zenodo.6514923
# *****************************************************************************

# ****** load required libraries *******
library(ranger)
library(parallel)
library(caret)
library(CAST)
source("simulation2_AGB/code/knndm_W.R")


# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "data"
outfolder <- "results"

startseed <- 1234567
n_samp    <- 100  # number of sample replicates (for each design)
cores <- 20 # number of cores for parallel computing

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/spatial")))
  dir.create(paste0(outfolder, "/spatial"))


# ************ FUNCTIONS ***************

predfun <- function(object, newdata){
  pred <- predict(object, newdata)
  pred[[1]]
}


spatialCV <- function(smpl, number, seed){

  fname <- paste0("AGB", "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  
  
  set.seed(seed)
  RMSE=W=R2=MAE <- numeric()
  
  coords <- data.frame(x=AGBdata$xcoord, y=AGBdata$ycoord)
  coords$clust <- stats::kmeans(x=coords, centers = 10)$cluster
  folds_spatial <- CAST::CreateSpacetimeFolds(coords, spacevar="clust",k=10)
  trControl <- trainControl(method="cv", index=folds_spatial$index, savePredictions=TRUE)
  
  AGBdata$glc2017 <- as.factor(AGBdata$glc2017)
  mtry <- 4
  pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
  RFmodel <- caret::train(agb~., AGBdata, respect.unordered.factors=TRUE, method = "ranger",
                          tuneGrid=pgrid, num.trees=500,trControl = trControl)
  
  
  err_spat <- global_validation(RFmodel)
  
  RMSE   <- error_spat[[1]]
  R2 <- err_spat[[2]]
  MAE <- error_spat[[3]]
  
  # calculate W statistic
  load(file.path(infolder, "ppoints.Rdata"))
  
  Gij <- c(FNN::knnx.dist(query = sf::st_coordinates(ppoints)[,1:2],
                          data = coords[,1:2], k = 1))
  
  
  W_spatial <- distclust_proj(coords, coords$clust)
  W_spatial <- twosamples::wass_stat(W_spatial, Gij) 
  
  W <- W_spatial
  
  seed <- seed + 1
  
  
  
  fname <-  paste0("AGB", "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out <- file.path(outfolder, "spatial", fname)
  save(RMSE, W, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    spatialCV(smpl, i, startseed)
  }
}, mc.cores = cores)
