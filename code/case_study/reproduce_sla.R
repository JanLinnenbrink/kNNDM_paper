#
#Reproduce study specific leaf area
setwd("C:/git/kNNDM_paper/")

library(caret)
library(CAST)
library(sf)
library(terra)
library(dplyr)
source("./code/reproduce_sla/reproduce_functions.R")

# load training points
sla <- st_read("./data/reproduce_sla/sla.gpkg") 
trainDat <- st_drop_geometry(sla) |> 
  dplyr::select(-coord.ID)

# train model
hyperparameter = expand.grid(mtry = 3,
                             splitrule = "variance",
                             min.node.size = 5)


r <- c("system.index","Total_Number")
predictor_names=names(trainDat)
predictor_names <- setdiff(predictor_names,r)

sla_kndm <- knndm(sla, co.ee,k=4, maxp=0.5, clustering = "hierarchical")
plot(sla_kndm)

model_spat = train_model("reproduced_knndmcv", sla, predictors = predictor_names,
                         response = "Total_Number",folds = sla_kndm$clusters,
                         hyperparameter = hyperparameter)
err_spat <- global_validation(model_spat)
k <- 2:10
sla_rand <- sample(1:10, nrow(sla), replace=TRUE)
model_rand <- train_model("reproduced_random",sla, predictors = predictor_names, response = "Total_Number",
                          folds = sla_rand, hyperparameter = hyperparameter)
err_rand <- global_validation(model_rand)

dist_general <- plot_geodist(sla,co.ee,
                             sampling="Fibonacci",
                             unit="km",
                             showPlot = FALSE)

dist_general <- dist_general$dist
dist_general$CV <- NaN
dist_general <- rbind(dist_general, dist_general)
dist_general[1:(nrow(dist_general)/2),"CV"] <- "random" 
dist_general[(nrow(dist_general)/2+1):nrow(dist_general),"CV"] <- "kNNDM" 


dist_rand <- plot_geodist(sla,co.ee,
                          sampling="Fibonacci",
                          cvfolds=sla_rand,
                          unit="km",
                          showPlot = FALSE)

dist_rand <- dist_rand$distances[dist_rand$distances$what=="CV-distances",]
dist_rand$CV <- "random"

dist_clstr <- plot_geodist(sla,co.ee,
                           sampling="Fibonacci",
                           cvfolds=sla_knndm$clusters,
                           unit="km",
                           showPlot = FALSE)
dist_clstr <- dist_clstr$distances[dist_clstr$distances$what=="CV-distances",]
dist_clstr$CV <- "kNNDM"

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

dist_all <- rbind(dist_general, dist_rand, dist_clstr)

dist_all$CV <- forcats::fct_relevel(dist_all$CV, "random","kNNDM")

(sla_dens <- ggplot() +
  stat_density(data=dist_all[dist_all$what%in% c("sample-to-sample", "sample-to-prediction"),],
               aes(x=dist, fill=what),
               geom="area",position="identity") +
  stat_density(data=dist_all[dist_all$what=="CV-distances",],aes(x=dist, color="CV-distances"), 
               linewidth=1.6, geom="line",position="identity") +
  scale_color_manual("", values="black") +
  scale_fill_manual("what",values=c(gg_color_hue(2))) +
  scale_x_log10() +
  xlab("geographic distance [km]")+
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", aspect.ratio=0.6) +
  facet_wrap(~CV, scales="fixed") +
  theme(strip.text.x = element_blank(),
        panel.spacing = unit(10, "lines")))

rmse_labs <- data.frame(
  label = c(paste0("RMSE = ", round(err_rand[[1]])), paste0("RMSE = ", round(err_spat[[1]]))),
  CV   = factor(c("random", "kNNDM")),
  x     = c(18^3, 18^3),
  y     = c(1.1, 1.1)
)
sla_dens <- sla_dens + geom_text(
  data    = rmse_labs,
  mapping = aes(x = x, y = y, label = label)
)

sla <- rbind(sla,sla)
sla$CV <- NaN
sla[1:(nrow(sla)/2),"CV"] <- sla_rand
sla[(nrow(sla)/2+1):nrow(sla),"CV"] <- sla_knndm$clusters
sla[1:(nrow(sla)/2),"what"] <- "random"
sla[(nrow(sla)/2+1):nrow(sla),"what"] <- "kNNDM"
sla$what <- forcats::fct_relevel(sla$what, "random","kNNDM")

(sla_map <- ggplot() +
    geom_sf(data=st_union(co.ee)) +
    geom_sf(data=sla, size=0.2, color=sla$CV) +
    theme_minimal()+
    facet_wrap(~what) +
    theme(strip.text.x = element_text(face="bold", size = 14)))


bp <- plot_grid(sla_map,sla_dens,nrow=2)
ggsave("figures/sla_plots.pdf",bp)
