# *****************************************************************************
# R Script implementing kNNDM  k-fold cross-validation.  
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
n_samp <- 100  # number of sample replicates
cores <- 20 # number of cores for parallel computing

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/knndm_caret_5k_3")))
  dir.create(paste0(outfolder, "/knndm_caret_5k_3"))


# ************ FUNCTIONS ***************
knndmCV <- function(smpl, number, variate, seed){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"knndm_caret_5k_3", fname)
  
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
  pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
  
  n <- length(pts_df$x)
  RMSE=R2=MAE <- numeric(n_CV)
  
  for(i_CV in 1:n_CV){
    
    set.seed(seed)
    knndm <- knndm(pts_sf, ppoints = ppoints, k = 10, maxp = 0.5)
    trControl = trainControl(method = "cv", savePredictions = "final",
                             index=knndm$indx_train)
    
    set.seed(seed)
    AGBdata$glc2017 <- as.factor(AGBdata$glc2017)
    mtry <- floor(sqrt(ncol(AGBdata[,-1])))
    pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
    RFmodel <- caret::train(agb~., AGBdata, repect.unordered.factors=TRUE, method = "ranger",
                            tuneGrid=pgrid, num.trees=500,
                            trControl = trControl)
      
    knndm_error <- global_validation(RFmodel)
    RMSE[i_CV] <- knndm_error[[1]]
    R2[i_CV] <- knndm_error[[2]]
    MAE[i_CV] <- knndm_error[[3]]      
    
    seed <- seed + 1
  } # loop over i_CV
  
  save(RMSE, R2, MAE, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  for(smpl in samples) {
    knndmCV(smpl, i, "AGB", startseed)
  }
}, mc.cores = cores)

