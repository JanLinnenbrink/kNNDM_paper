# *****************************************************************************
# R Script implementing random  k-fold cross-validation.  
# Modified from https://doi.org/10.5281/zenodo.6514923
# *****************************************************************************

# ****** load required libraries *******
#.libPaths("~/r_packages/")

library(ranger)
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)


# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")

infolder <- "simulation2_AGB/data"
outfolder <- "simulation2_AGB/results"
startseed <- 1234567
n_samp <- 100 # number of sample replicates (for each design)
cores <- 10 # number of cores for parallel computing

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/random")))
  dir.create(paste0(outfolder, "/random"))


# ************ FUNCTIONS ***************
knndmCV <- function(smpl, number, variate, seed){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"random", fname)
  
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  RMSE=R2=MAE <- numeric()
  
  pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
  pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))

  
  # RMSE random caret
  set.seed(seed)
  mtry <- floor(sqrt(ncol(AGBdata[,-1])))
  pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
  RFmodel_rand <- caret::train(agb~., AGBdata,respect.unordered.factors=TRUE, method = "ranger",
                               tuneGrid=pgrid,num.trees=500,
                               trControl = trainControl(method = "cv", savePredictions = "final"))
  
  
  err_rand_caret <- global_validation(RFmodel_rand)
  RMSE <- err_rand_caret[[1]]
  
  save(RMSE, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
  }
}, mc.cores = cores)

