# First two from Ludwig et al. (2023): Assessing and improving the transferability of current global spatial prediction models 
# https://github.com/LOEK-RS/global_applicability/tree/main/sla

train_model = function(modelname, training_samples, predictors, response, folds, hyperparameter){
  
  training_samples = training_samples %>% 
    dplyr::select(all_of(c(predictors, response))) %>% 
    st_drop_geometry()
  
  i = fold2index(folds)
  
  
  set.seed(7353)
  rfmodel = caret::train(x = training_samples %>% dplyr::select(all_of(predictors)),
                         y = training_samples %>% dplyr::pull(response),
                         method = "ranger",
                         tuneGrid = hyperparameter,
                         num.trees = 300,
                         trControl = trainControl(method = "cv", number = length(unique(folds)),
                                                  index = i$index, indexOut = i$indexOut,
                                                  savePredictions = "final"),
                         importance = "impurity")
  
  return(rfmodel)
  
}

fold2index = function(fold){
  
  fold = data.frame(fold = fold)
  
  indOut = fold %>% dplyr::group_by(fold) %>%
    attr('groups') %>% dplyr::pull(.rows)
  
  ind = purrr::map(seq(length(indOut)), function(x){
    s = seq(nrow(fold))
    s = s[!s %in% indOut[[x]]]
    return(s)
  })
  return(
    list(
      index = ind,
      indexOut = indOut
    )
  )
  
}


dist_rand <- function(x, folds) {
  
  distclust <- function(distm, folds){
    alldist <- rep(NA, length(folds))
    for(f in unique(folds)){
      alldist[f == folds] <- apply(distm[f == folds, f != folds, drop=FALSE], 1, min)
    }
    alldist
  }
  
  distmat <- sf::st_distance(x)
  units(distmat) <- NULL
  diag(distmat) <- NA
  Gj <- apply(distmat, 1, function(y) min(y, na.rm=TRUE))
  Gij <- sf::st_distance(ppoints, x)
  units(Gij) <- NULL
  Gij <- apply(Gij, 1, min)
  Gjstar <- distclust(distmat, folds)
  
  W <- twosamples::wass_stat(Gjstar, Gij)
  
  r <- list(Gj=Gj, Gij=Gij, Gjstar=Gjstar, W=W)
  class(r) <- c("knndm", "list")
  r
  
}
