# Libraries and utils ----
library("tidyverse")
library("raster")
library("terra")
library("CAST")
library("sf")
library("caret")
library("parallel")
library("doParallel")
library("pbapply")
source("code/sim_utils.R")

# No need for proj4 warnings
options("rgdal_show_exportToProj4_warnings"="none")

#' Sample simulation: virtual species.
#' @details
#' Simulates a series of sampling points for the simulation.
#' @param nsamples Integer. Number of samples to simulate.
#' @param dsamples Character. Spatial distribution of the samples. 5 are
#' possible: sregular, wregular, random, wclust, sclust.
#' @param sarea sf/sfc polygon where samples will be simulated.
sim2_samples <- function(nsamples, dsamples, sarea){


  if(dsamples=="sregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 40000)
  }else if(dsamples=="wregular"){
    simpoints <- jitterreg_sample(sarea, nsamples, 80000)
  }else if(dsamples=="random"){
    simpoints <- st_sample(sarea, nsamples)
  }else if(dsamples=="wclust"){
    simpoints <- clustered_sample(sarea, nsamples, 25, 80000)
  }else if(dsamples=="sclust"){
    simpoints <- clustered_sample(sarea, nsamples, 10, 80000)
  }

  simpoints <- st_sf(geometry=simpoints)
  simpoints
}


#' Fits a RF model and evaluates it using several methods (species simulation).
#' @details
#' Fits a RF model and evaluates it using random 10-fold CV, spatial 10-fold CV,
#' NNDM LOO CV, kNNDM 10-fold CV and true errors.
#' @param form String. Model formula.
#' @param spatial_index list. Indices for spatial CV.
#' @param knndm_indexTrain list. Indices for kNNDM 10-fold CV training data.
#' @param knndm_indexTest list. Indices for kNNDM 10-fold CV test data.
#' @param nndm_indexTrain list. Indices for NNDM LOO CV test data.
#' @param nndm_indexTest list. Indices for NNDM LOO CV training data.
#' @param pgrid Data frame. Parameter rand of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param surfdf Data frame. Surface data.
fitval_rf_species <- function(form,
                              knndm_indexTrain,
                              knndm_indexTest,
                              pgrid, traindf,
                              surfdf) {


  # Validate with random CV and compute metrics
  r_cntrl <- trainControl(method="CV", savePredictions=TRUE)
  rand_mod <- train(form, data=traindf, method="rf",
                   trControl=r_cntrl, tuneGrid=pgrid, ntree=100)
  err_stats_rand <- global_validation(rand_mod)
  rand_stats <- data.frame(RMSE = err_stats_rand[[1]],
              MAE = err_stats_rand[[3]],
              R2 = err_stats_rand[[2]])
  names(rand_stats) <- paste0(names(rand_stats), "_rand")

  # Compute CV statistics in surface
  surfdf$preds <- predict(rand_mod, newdata=surfdf)
  surf_stats <- surfdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(surf_stats) <- paste0(names(surf_stats), "_surf")

 
  # Validate with knndm and compute CV statistics
  knndm_cntrl <- trainControl(method="cv",
                                index=knndm_indexTrain,
                                indexOut=knndm_indexTest,
                                savePredictions=TRUE)
  knndm_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=knndm_cntrl, tuneGrid=pgrid, ntree=100))
  err_stats_knndm <- global_validation(knndm_mod)
  knndm_stats <- data.frame(RMSE = err_stats_knndm[[1]],
              MAE = err_stats_knndm[[3]],
              R2 = err_stats_knndm[[2]])
  names(knndm_stats) <- paste0(names(knndm_stats), "_knndm")


  # Tidy and return results
  data.frame(surf_stats,knndm_stats)
}


#' Simulation function.
#' @details
#' The function takes a virtual species simulated landscape for the Iberian
#' peninsula using bioclim data, simulates sampling points
#' and fits a RF and evaluates it using spatial 10-fold CV, random 10-fold CV,
#' kNNDM 10-fold CV and NNDM LOO CV; as well as the true error.
#' @param rgrid sf or sfc point object. Prediction grid.
#' @param rstack SpatRaster object. Contains the landscape data and the
#' simulated outcome.
#' @param sampling_area sf or sfc polygon object. Sampling area.
#' @param sample_dist String or vector or string. Distribution of the sampling
#' points. 5 are possible: "sregular", "wregular", "random", "wclust","sclust".
sim_species <- function(rgrid, rstack, sampling_area, kseq = seq(2,20,2),
                        sample_dist=c("sregular", "wregular", "random",
                                      "wclust","sclust")){
  
  # sample prediction points from the prediction area
  ppoints <- st_sample(sampling_area, 1000, type="regular")

  # Initiate results object and fixed information for all models
  res_k <- data.frame()
  res <- list()
  grid_data <- as.data.frame(terra::extract(rstack, terra::vect(rgrid)))
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  pgrid<- data.frame(mtry=6)
  
  
  grid <- expand.grid(sample_dist, kseq)
  
  res <- mapply(function(dist_it, k) {
    
    # Simulate sampling points according to parameters and constraints
    train_points <- sim2_samples(100, dist_it, sampling_area) |> st_transform(st_crs(ppoints))
    
    # Get training and surface data for modelling and validation
    train_data <- terra::extract(rstack, train_points)
    
    # Estimate outcome range
    train_points$outcome <- train_data$outcome
    # Define folds based on kNNDM
    folds_knndm <- CAST::knndm(train_points, ppoints = ppoints, k = k, maxp = (1/k)+0.4)
    
    #### Model fitting and validation
    mod <- fitval_rf_species(form,
                             folds_knndm$indx_train,
                             folds_knndm$indx_test,
                             pgrid, train_data,
                             grid_data)
    
    #### Compute W statistics for knndmCV
    W_knndmCV <- folds_knndm["W"][[1]]
    
    #### store results
    W_df <- data.frame(W_knndmCV = W_knndmCV)
    
    # Store results of the number of k
    res_k_it <- cbind(data.frame(k=k, dsample=dist_it, stringsAsFactors = FALSE),
                      mod, W_df)
    res_k_it
  }, grid$Var1, grid$Var2, SIMPLIFY = FALSE)

  
  res <- do.call(rbind.data.frame, res)
  row.names(res) <- NULL
  res
  
}
