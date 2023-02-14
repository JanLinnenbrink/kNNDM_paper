#
#Reproduce study specific leaf area
setwd("C:/git/kNNDM_paper/")

library(caret)
library(CAST)
library(sf)
library(terra)
library(dplyr)
library("cowplot")
source("./code/case_study/reproduce_utils.R")

# load training points
sla <- st_read("data/case_study/sla.gpkg") 
trainDat <- st_drop_geometry(sla) |> 
  dplyr::select(-coord.ID)

# load prediction points
ppoints <- st_read("data/case_study/ppoints.gpkg")

# train model
hyperparameter = expand.grid(mtry = 3,
                             splitrule = "variance",
                             min.node.size = 5)

predictor_names <- c("B5","EVIstd","B6","EVImax","alt","B2", "bio3",
                     "bio19","B4","bio4","bio9","B7","NDWIstd", "B3","bio6")

# run kNNDM to obtain folds
library("future")
library("future.apply")
plan(multisession, workers = 10)

k_s = 2:10
sla_kndm_k <- future_lapply(k_s, function(x) knndm(sla,ppoints=ppoints,k=x,maxp=0.8,clustering="hierarchical"))
k_sel_idx <- which.min(unlist(lapply(sla_kndm_k, function(x) x$W)))
sla_kndm <- sla_kndm_k[[k_sel_idx]]
sla_kndm$W
plot.knndm(sla_kndm)
saveRDS(sla_kndm,"results/sla_knndm.RDS")

sla_kndm <- knndm(sla,ppoints=ppoints,maxp=0.8,clustering="hierarchical")
saveRDS(sla_kndm,"results/sla_knndm.RDS")

# reproduce model with those folds
model_knndm <- train_model(training_samples=sla, predictors = predictor_names,
                         response = "SLA",folds = sla_kndm$clusters,
                         hyperparameter = hyperparameter)
err_knndm <- global_validation(model_knndm)
saveRDS(err_knndm,"results/err_knndm.RDS")
saveRDS(model_knndm, "results/sla_model_knndm.RDS")

# reproduce model with random CV split
sla_rand <- sample(1:10, nrow(sla), replace=TRUE)
rand_dist <- dist_rand(st_as_sfc(sla), folds=sla_rand)
saveRDS(rand_dist, "results/rand_dist.RDS")
model_rand <- train_model(training_samples=sla, predictors = predictor_names, response = "SLA",
                          folds = sla_rand, hyperparameter = hyperparameter)
err_rand <- global_validation(model_rand)
saveRDS(sla_rand, "results/folds_rand.RDS")
saveRDS(err_rand, "results/err_rand.RDS")
saveRDS(model_rand, "results/sla_model_random.RDS")

