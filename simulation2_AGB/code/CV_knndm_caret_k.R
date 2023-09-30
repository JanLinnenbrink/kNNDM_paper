# *****************************************************************************
# R Script testing the influence of different k on the outcome of
# kNNDM k-fold cross-validation.  
# Modified from doi
# *****************************************************************************

# ****** load required libraries *******
#.libPaths("~/r_packages/")

setwd("simulation2_AGB/")
library(CAST)
library(sf)
library(raster)
library(caret)
library(parallel)

# ************ GLOBALS ***************
infolder <- "data"
outfolder <- "results"
startseed <- 1234567
n_samp <- 100  # number of sample replicates
cores <- 20 # number of cores for parallel computing
num_k <- seq(2,20,by=2) # tested numbers of k

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/knndm_k_cstrong")))
  dir.create(paste0(outfolder, "/knndm_k_cstrong"))

# ************ FUNCTIONS ***************


knndmCV <- function(smpl, number, variate, seed, num_k){
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"knndm_k_cstrong", fname)
  
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
  pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))

  n <- length(pts_df$x)
  
  
  k_err <- data.frame(RMSE=NaN, W=NaN, k=num_k)
  for(k in num_k){
    
    set.seed(seed)
    kout <- knndm(pts_sf, ppoints = ppoints, k = k, maxp = (1/k)+0.1)
    trControl = trainControl(method = "cv", savePredictions = "final",
                             index=kout$indx_train)
    
    AGBdata$glc2017 <- as.factor(AGBdata$glc2017)
    mtry <- floor(sqrt(ncol(AGBdata[,-1])))
    pgrid <- expand.grid(splitrule="variance",min.node.size=5,mtry=mtry)
    RFmodel <- caret::train(agb~., AGBdata, repect.unordered.factors=TRUE, method = "ranger",
                            tuneGrid=pgrid, num.trees=500,
                            trControl = trControl)
    
    knndm_error <- global_validation(RFmodel)
    k_err[k_err$k==k,"RMSE"] <- knndm_error[[1]]
    k_err[k_err$k==k,"W"] <- kout$W
    
    seed <- seed + 1
  }
  
  save(k_err, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  knndmCV("clusterStrong", i, "AGB", startseed, num_k)
}, mc.cores = cores)

