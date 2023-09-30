# *****************************************************************************
# R Script for calculating the W statistic and RMSE associated with each
# split in the kNNDM function.
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

source("simulation1_virtualSpecies/code/knndm_W.R")

# ************ GLOBALS ***************

infolder <- "data"
outfolder <- "results"
startseed <- 1234567
n_samp <- 100  # number of sample replicates (for each design)
cores <- 30 # number of cores for parallel computing

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/knndm_W_AGB_clStr")))
  dir.create(paste0(outfolder, "/knndm_W_AGB_clStr"))


# ************ FUNCTIONS ***************


knndmCV <- function(smpl, number, variate, seed){
  
  
  fname  <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out  <- file.path(outfolder,"knndm_W_AGB_cl_Str", fname)
  fname <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  f_in <- file.path(infolder,smpl,fname)
  load(f_in)
  
  message(f_out) 
  # load ppoints
  load(file.path(infolder, "ppoints.Rdata"))
  
  
  pts_df <- data.frame(x=AGBdata$xcoord * 1000, y=AGBdata$ycoord * 1000)
  pts_sf <- st_as_sf(pts_df, coords = c("x","y"), crs = st_crs(ppoints))
  
  
  n <- length(pts_df$x)
  
  set.seed(seed)
  kndm_out <- knndmW(pts_sf, ppoints = ppoints, k = 2, maxp = 0.8)
  
  # Validate with knndm
  form = agb~.; data=AGBdata
  
  kndm_folds <- kndm_out$clusters
  W <- kndm_out$W
  message("ksplit completed") 
  if (class(kndm_folds)=="list") {
    kndm_stats <- lapply(seq_along(kndm_folds), function(i) {
      fdf <- data.frame(f=kndm_folds[[i]])
      ctrl <- CAST::CreateSpacetimeFolds(fdf, spacevar=1, k = max(fdf))
      kndm_ctrl <- trainControl(method="cv",
                                index=ctrl$index,
                                indexOut=ctrl$indexOut,
                                savePredictions=TRUE)
      kndm_out_mod <- suppressWarnings( # train() can't compute R2
        train(form, data, respect.unordered.factors=TRUE, method = "ranger",
              trControl=kndm_ctrl))
      kndm_stats <- kndm_out_mod$pred |> 
        dplyr::summarise(RMSE = sqrt(mean((obs-pred)^2)),
                         MAE = mean(abs(obs-pred)),
                         R2 = cor(obs, pred)^2, W=W[[i]])
      names(kndm_stats) <- paste0(names(kndm_stats), "_kndm")
      kndm_stats
    })
    kndm_stats <- do.call(rbind, kndm_stats)
  } else {
    fdf <- data.frame(f=kndm_folds)
    ctrl <- CAST::CreateSpacetimeFolds(fdf, spacevar=1, k = max(fdf))
    kndm_ctrl <- trainControl(method="cv",
                              index=ctrl$index,
                              indexOut=ctrl$indexOut,
                              savePredictions=TRUE)
    kndm_out_mod <- suppressWarnings( # train() can't compute R2
      train(form, data, respect.unordered.factors=TRUE, method = "ranger",
            trControl=kndm_ctrl))
    kndm_stats <- kndm_out_mod$pred |> 
      dplyr::summarise(RMSE = sqrt(mean((obs-pred)^2)),
                       MAE = mean(abs(obs-pred)),
                       R2 = cor(obs, pred)^2, W=W)
    names(kndm_stats) <- paste0(names(kndm_stats), "_kndm")
  }
  save(kndm_stats, file=f_out)
  
}


# ************ CALL THE FUNCTIONS ************ 
mclapply(seq(n_samp), function(i) {
  
  knndmCV("clusterStrong", i, "AGB", startseed)
  
}, mc.cores = cores)

