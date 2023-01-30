# Libraries and utils ----
library("tidyverse")
library("raster")
library("terra")
library("CAST")
library("sf")
library("caret")
library("gstat")
library("virtualspecies")
library("parallel")
library("doParallel")
library("pbapply")
setwd("C:/0_Msc_Loek/M7_Fernerkundung/NNDMpaper-main/")
source("./code/sim_utils.R")

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


# helper function: calculate AOA percentages
AOA_perc <- function(rstack, model, indexTrain=NULL, indexTest=NULL) {
  AOA <- aoa(rstack, model=model, CVtrain = indexTrain, CVtest = indexTest)$AOA
  AOAn <- terra::freq(AOA)[,2:3]
  overall_cells <- sum(AOAn$count)
  if(nrow(AOAn) == 1) {
    list(data.frame(inside_AOA = 100, outside_AOA = 0), AOA)
  } else {
    AOAp <- AOAn$count / overall_cells * 100
    list(data.frame(inside_AOA = AOAp[[2]], outside_AOA = AOAp[[1]]), AOA)
  }
}

# helper function: calculate RMSE inside AOA
AOA_err <- function(rstack, AOA, rand_mod) {
  rstack_curr <- rstack
  rstack_curr[AOA<1] <- NA
  rgrid_curr <- as.points(rstack_curr)
  surfdf_curr <- as.data.frame(terra::extract(rstack_curr, rgrid_curr))
  surfdf_curr$preds <- predict(rand_mod, newdata=surfdf_curr)
  surf_curr_stats <- surfdf_curr %>%
    summarise(RMSE_AOA = sqrt(mean((outcome-preds)^2)),
              MAE_AOA = mean(abs(outcome-preds)),
              R2_AOA = cor(outcome, preds)^2)
  names(surf_curr_stats) <- paste0(names(surf_curr_stats), "_surf")
  surf_curr_stats
}


#' Fits a RF model and evaluates it using several methods (species simulation).
#' @details
#' Fits a RF model and evaluates it using LOO CV, bLOO CV, NNDM LOO CV and true errors.
#' @param form String. Model formula.
#' @param folds_rand list. Indices for rand CV
#' training data.
#' @param folds_kndm list. Indices for kndm CV
#' data.
#' @param ndmout_indexTrain list. Indices for NNDM LOO CV (outcome range) training
#' data.
#' @param ndmout_indexTest list. Indices for NNDM LOO CV (outcome range) test data.
#' @param ndmres_indexTrain list. Indices for NNDM LOO CV (residual range) training
#' data.
#' @param ndmres_indexTest list. Indices for NNDM LOO CV (residual range) test data.
#' @param prand Data frame. Parameter rand of the model.
#' @param traindf Data frame. Training data to fit the model.
#' @param surfdf Data frame. Surface data.
fitval_rf_species <- function(form,
                              spatial_index,
                              kndm_indexTrain,
                              kndm_indexTest,
                              nndm_indexTrain,
                              nndm_indexTest,
                              prand, traindf,
                              surfdf,train_points,
                              rstack) {


  # Validate with random CV and compute metrics
  print("random CV")
  r_cntrl <- trainControl(method="CV", savePredictions=TRUE)
  rand_mod <- train(form, data=traindf, method="rf",
                   trControl=r_cntrl, tunerand=prand, ntree=100)
  err_stats_rand <- global_validation(rand_mod)
  rand_stats <- data.frame(RMSE = err_stats_rand[[1]],
              MAE = err_stats_rand[[3]],
              R2 = err_stats_rand[[2]])
  AOA_rand <- AOA_perc(rstack,rand_mod)
  AOA_rand_err <- AOA_err(rstack, AOA_rand[[2]],rand_mod)
  rand_stats <- cbind(rand_stats, AOA_rand_err)
  rand_stats$in_AOA <- AOA_rand[[1]]$inside_AOA
  names(rand_stats) <- paste0(names(rand_stats), "_rand")

  # Compute CV statistics in surface
  print("surface CV")
  surfdf$preds <- predict(rand_mod, newdata=surfdf)
  surf_stats <- surfdf %>%
    summarise(RMSE = sqrt(mean((outcome-preds)^2)),
              MAE = mean(abs(outcome-preds)),
              R2 = cor(outcome, preds)^2)
  names(surf_stats) <- paste0(names(surf_stats), "_surf")

  # Validate with 10-fold cluster CV and compute CV statistics
  print("spatial CV")
  spatial_cntrl <- trainControl(method="cv", index=spatial_index,
                                 savePredictions=TRUE)
  spatial_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=spatial_cntrl, tunerand=prand, ntree=100))
  err_stats_spatial <- global_validation(spatial_mod)
  spatial_stats <- data.frame(RMSE = err_stats_spatial[[1]],
              MAE = err_stats_spatial[[3]],
              R2 = err_stats_spatial[[2]])
  AOA_spatial <- AOA_perc(rstack, spatial_mod)
  AOA_spatial_err <- AOA_err(rstack,AOA_spatial[[2]],rand_mod)
  spatial_stats <- cbind(spatial_stats, AOA_spatial_err)
  spatial_stats$in_AOA <- AOA_spatial[[1]]$inside_AOA
  names(spatial_stats) <- paste0(names(spatial_stats), "_spatial")

  # Validate with knndm and compute CV statistics
  print("knndm CV")
  kndm_cntrl <- trainControl(method="cv",
                                index=kndm_indexTrain,
                                indexOut=kndm_indexTest,
                                savePredictions=TRUE)
  time_kndm_mod <- system.time(
    kndm_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=kndm_cntrl, tunerand=prand, ntree=100))
    )[[3]]
  err_stats_kndm <- global_validation(kndm_mod)
  kndm_stats <- data.frame(RMSE = err_stats_kndm[[1]],
              MAE = err_stats_kndm[[3]],
              R2 = err_stats_kndm[[2]])
  AOA_kndm <- AOA_perc(rstack,kndm_mod, indexTrain = kndm_indexTrain, indexTest = kndm_indexTest)
  AOA_kndm_err <- AOA_err(rstack,AOA_kndm[[2]],rand_mod)
  kndm_stats <- cbind(kndm_stats, AOA_kndm_err)
  kndm_stats$in_AOA <- AOA_kndm[[1]]$inside_AOA
  names(kndm_stats) <- paste0(names(kndm_stats), "_kndm")

  # Validate with NNDM LOO and compute CV statistics (outcome range)
  print("NNDM CV")
  nndm_cntrl <- trainControl(method="cv",
                                index=nndm_indexTrain,
                                indexOut=nndm_indexTest,
                                savePredictions=TRUE)
  time_ndm_mod <- system.time(
    nndm_mod <- suppressWarnings( # train() can't compute R2
    train(form, data=traindf, method="rf",
          trControl=nndm_cntrl, tunerand=prand, ntree=100))
  )[[3]]
  err_stats_nndm <- global_validation(nndm_mod)
  nndm_stats <- data.frame(RMSE = err_stats_nndm[[1]],
              MAE = err_stats_nndm[[3]],
              R2 = err_stats_nndm[[2]])
  AOA_ndm <- AOA_perc(rstack,nndm_mod, indexTrain = nndm_indexTrain, indexTest = nndm_indexTest)
  AOA_ndm_err <- AOA_err(rstack,AOA_ndm[[2]],rand_mod)
  nndm_stats <- cbind(nndm_stats, AOA_ndm_err)
  nndm_stats$in_AOA <- AOA_ndm[[1]]$inside_AOA
  names(nndm_stats) <- paste0(names(nndm_stats), "_nndm")

  # Tidy and return results
  data.frame(rand_stats,surf_stats,spatial_stats,kndm_stats,nndm_stats,
             time_kndm_mod, time_ndm_mod)
}


#' Simulation function 2: virtual species.
#' @details
#' The function takes a virtual species simulated landscape for the Iberian
#' peninsula using bioclim data, simulates sampling points
#' and fits a RF and evaluates it using spatial 10-fold CV, random 10-fold CV,
#' kNNDM CV and NNDM LOO CV; as well as the true error.
#' @param rrand sf or sfc point object. Prediction rand.
#' @param rstack terra SpatRaster object. Contains the landscape data and the
#' simulated outcome.
#' @param sampling_area sf or sfc polygon object. Sampling area.
#' @param sample_dist String or vector or string. Distribution of the sampling
#' points. 5 are possible: "sregular", "wregular", "random", "wclust","sclust".
sim_species <- function(rrand, rstack, sampling_area,
                        sample_dist=c("sregular", "wregular", "random",
                                      "wclust","sclust")){
  
  # sample prediction points from the prediction area
  ppoints <- st_sample(sampling_area, 1000, type="regular")

  # Initiate results object and fixed information for all models
  res <- data.frame()
  rand_data <- as.data.frame(terra::extract(rstack, terra::vect(rrand)))
  form <- as.formula(paste0("outcome~", paste0("bio", 1:19, collapse="+")))
  prand <- data.frame(mtry=6)

  i <- 0
  # Start sampling loop
  for(dist_it in sample_dist){
    i=i+1
    # Simulate sampling points according to parameters and constraints
    train_points <- sim2_samples(100, dist_it, sampling_area) |> st_transform(st_crs(ppoints))

    # Get training and surface data for modelling and validation
    train_data <- terra::extract(rstack, train_points)

    # cluster points based on kmeans
    coords <- sf::st_coordinates(train_points) |> as.matrix()
    clust <- stats::kmeans(x=coords, centers = 10)
    pts_id <- cbind(train_points, clust$cluster)
    pts_id$rand <- sample(c(1:10),nrow(pts_id), replace = TRUE)

    #### Folds based on outcome range
    # Estimate outcome range
    train_points$outcome <- train_data$outcome
    # Define random folds
    folds_random <- CAST::CreateSpacetimeFolds(pts_id, spacevar="rand", k=10)
    # Define folds based on 10-fold cluster CV
    folds_spatial <- CAST::CreateSpacetimeFolds(pts_id, spacevar="clust.cluster",k=10)
    # Define folds based on kNNDM
    time_kndm <- system.time(folds_kndm <- CAST::knndm(train_points, ppoints = ppoints))[[3]]
    # Define folds based on NNDM
    time_ndm <- system.time(folds_ndm <- CAST::nndm(train_points, ppoints = ppoints, min_train =  0.5))[[3]]
    #### Model fitting and validation
    mod <- fitval_rf_species(form,
                             folds_spatial$index,
                             folds_kndm$indx_train,
                             folds_kndm$indx_test,
                             folds_ndm$indx_train,
                             folds_ndm$indx_test,
                             prand, train_data,
                             rand_data,train_points,rstack)

    # Store results of the iteration
    res_it <- cbind(data.frame(dsample=dist_it, time_ndm=time_ndm, time_kndm=time_kndm,
                               stringsAsFactors = FALSE), mod)
    res <- bind_rows(res, res_it)
  }

  row.names(res) <- NULL
  res
}
