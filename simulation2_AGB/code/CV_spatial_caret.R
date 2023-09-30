# *****************************************************************************
# R Script implementing spatial  k-fold cross-validation.  
# Modified from doi
# *****************************************************************************

# ****** load required libraries *******
#.libPaths("~/r_packages/")

setwd("simulation2_AGB/")
library(ranger)
library(parallel)
library(caret)
library(CAST)
source("simulation1_virtualSpecies/code/knndm_W.R")


# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "data"
outfolder <- "results"

startseed <- 1234567
n_CV      <- 3  # number of cross validation replications
n_samp    <- 100  # number of sample replicates (for each design)
cores <- 20 # number of cores for parallel computing

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/spatial_caret_5k_W")))
  dir.create(paste0(outfolder, "/spatial_caret_5k_W"))


# ************ FUNCTIONS ***************

predfun <- function(object, newdata){
  pred <- predict(object, newdata)
  pred[[1]]
}


spatialCV <- function(smpl, number, variate, seed){

  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  for(i_CV in 1:n_CV){
    
    set.seed(seed)
    RMSE=W <- numeric(n_CV)
    
    coords <- data.frame(x=AGBdata$xcoord, y=AGBdata$ycoord)
    coords$clust <- stats::kmeans(x=coords, centers = 10)$cluster
    folds_spatial <- CAST::CreateSpacetimeFolds(coords, spacevar="clust",k=10)
    trControl <- trainControl(method="cv", index=folds_spatial$index, savePredictions=TRUE)
    
    AGBdata$glc2017 <- as.factor(AGBdata$glc2017)
    mtry <- 4
    pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
    RFmodel <- caret::train(agb~., AGBdata, respect.unordered.factors=TRUE, method = "ranger",
                            tuneGrid=pgrid, num.trees=500,trControl = trControl)
      

    error_spat <- global_validation(RFmodel)
    
    RMSE[i_CV]   <- error_spat[[1]]
    
    # calculate W statistic
    load(file.path(infolder, "ppoints.Rdata"))
    
    Gij <- c(FNN::knnx.dist(query = sf::st_coordinates(ppoints)[,1:2],
                            data = coords[,1:2], k = 1))
    
    
    W_spatial <- distclust_proj(coords, coords$clust)
    W_spatial <- twosamples::wass_stat(W_spatial, Gij) 
    
    W[i_CV] <- W_spatial
    
    seed <- seed + 1
  }
  
  
  fname <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out <- file.path(outfolder, "spatial_caret_5k_W", fname)
  save(RMSE, W, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    spatialCV(smpl, i, "AGB", startseed)
  }
}, mc.cores = cores)
