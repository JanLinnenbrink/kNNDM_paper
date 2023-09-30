# *****************************************************************************
# R Script implementing random  k-fold cross-validation.  
# Modified from doi
# *****************************************************************************

# ****** load required libraries *******
#.libPaths("~/r_packages/")

setwd("simulation2_AGB/")
library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)


# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "data"
outfolder <- "results"
startseed <- 1234567
n_CV   <- 3  # number of cross validation replications
n_samp <- 100 # number of sample replicates (for each design)
cores <- 20 # number of cores for parallel computing

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/random_caret_5k")))
  dir.create(paste0(outfolder, "/random_caret_5k"))


# ************ FUNCTIONS ***************
knndmCV <- function(smpl, number, variate, seed){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"random_caret_5k", fname)
  
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  RMSE=R2=MAE <- numeric(n_CV)
  
  pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
  pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))

  for(i_CV in 1:n_CV){
    
    SSR <- 0
    SST <- 0
    
    ntrees <- 500
    
    # RMSE random caret
    set.seed(seed)
    mtry <- floor(sqrt(ncol(AGBdata[,-1])))
    pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
    RFmodel_rand <- caret::train(agb~., AGBdata,respect.unordered.factors=TRUE, method = "ranger",
                                 tuneGrid=pgrid,num.trees=500,
                                 trControl = trainControl(method = "cv", savePredictions = "final"))
      

    err_rand_caret <- global_validation(RFmodel_rand)
    RMSE[i_CV] <- err_rand_caret[[1]]
    R2[i_CV] <- err_rand_caret[[2]]
    MAE[i_CV] <- err_rand_caret[[3]]
    
  } # loop over i_CV
  save(RMSE,R2,MAE, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
  }
}, mc.cores = cores)

