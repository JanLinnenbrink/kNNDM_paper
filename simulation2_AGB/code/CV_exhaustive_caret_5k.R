# *****************************************************************************
# R Script implementing exhaustive validation, which is used as the reference.  
# Modified from doi
# Due to the long run-time, this script calculates the reference for 1 realization 
# of each sampling design, and thus this script needs to be run 100 times.
# *****************************************************************************

# ****** load required libraries *******
#Sys.sleep(round(runif(1, min = 1, max = 240)))

#.libPaths("~/r_packages/")

setwd("simulation2_AGB/")

library(caret)
library(CAST)
library(terra)
library(parallel)

# ************ GLOBALS ***************
samples   <- c("clusterMedium", "clusterStrong", "clusterGapped", "regular", 
               "simpleRandom")
infolder <- "data"
outfolder <- "results"
startseed <- 1234567

# create outfolders if they don't exist
if(!dir.exists(outfolder))
  dir.create(outfolder)

if(!dir.exists(paste0(outfolder, "/exhaustive_caret_5k")))
  dir.create(paste0(outfolder, "/exhaustive_caret_5k"))

if(!file.exists(file.path(outfolder, "exhaustive_caret_5k", "runs.csv"))) {
  write.csv(data.frame("runs"=1), file.path(outfolder, "exhaustive_caret_5k", "runs.csv"), row.names = FALSE)
}

csv_file <- file.path(outfolder, "exhaustive_caret_5k", "runs.csv")
runs <- read.csv(csv_file)
lastIndex <- runs[nrow(runs),1]
thisIndex <- lastIndex + 1
print(paste0("this Index is: ", thisIndex))
runs[thisIndex,1] <- thisIndex
write.csv(runs, file = csv_file, row.names = FALSE)

# download data from https://doi.org/10.5281/zenodo.6513429
# ****** load input raster data ******
msk <- rast(file.path(infolder, "TOTmask.tif"))
AGBstack <- rast(file.path(infolder, "AGBstack.tif"))

AGBstack$glc2017 <- terra::as.factor(AGBstack$glc2017)

AGt <- rast(file.path(infolder, "agb.tif"))
AGB <- mask(AGt, msk, filename="tmpagb.tif", overwrite=T)
rm(AGt, OCt)

# ************ FUNCTIONS ***************
rmsefu <- function(ref, pred){
  residsq <- (ref - pred)^2
  sqrt(global(residsq, "mean", na.rm=T)[[1]])
}

exhaustive <- function(smpl, number, variate, seed){
  
  message(number)
  
  fname1 <- paste0(variate, "data", sprintf("%03d", number), ".Rdata")
  fname2 <- paste0(variate, smpl, sprintf("%03d", number), ".tif")
  
  f_in  <- file.path(infolder,smpl,fname1)
  f_out <- file.path(outfolder, "exhaustive_caret_5k", fname2)
  
  load(f_in)
  
  set.seed(seed)  
  
  mtry <- 4
  pgrid <- expand.grid(mtry=mtry, splitrule="variance", min.node.size=5)
  RFmodel <- caret::train(agb~., AGBdata[,!names(AGBdata) %in% c("ID")],
                          respect.unordered.factors=TRUE, method = "ranger",
                          tuneGrid=pgrid, num.trees=500)
  map  <- terra::predict(AGBstack, RFmodel, filename=f_out, overwrite=TRUE, na.rm=TRUE)
  RMSE <- rmsefu(AGB, map); message(paste("RMSE:", RMSE))
  
  fname <-  paste0(variate, "_", smpl, sprintf("%03d", number), ".Rdata")
  f_out2 <- file.path(outfolder,"exhaustive_caret_5k", fname)
  message(f_out2)
  save(RMSE, file=f_out2)
  file.remove(f_out)
}


# ************ CALL THE FUNCTIONS ************ 
lapply(samples, function(smpl) {
  i <- thisIndex 
  exhaustive(smpl, i, "AGB", startseed)
})

