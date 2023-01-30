library(terra)
library(sf)
library(CAST)
library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(tidyr)
library(forcats)
library(ggthemes)
library(cowplot)
library(caret)

#library(ggspatial)

setwd("C:/0_Msc_Loek/M7_Fernerkundung/NNDMpaper-main/example/")
#predictors <- terra::rast("gridded_data_for_AQbench.nc") 
#crs(predictors) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

#plot(predictors[["max_nightlight_25km"]])
#names(predictors)

stations <- read.csv("AQbench_dataset.csv") |> 
  st_as_sf(coords=c("lon","lat"), crs="epsg:4326")

predictors <- names(stations)[4:35]
response <- "o3_average_values"

ee <- st_crs("+proj=eqearth")
co <- ne_countries(returnclass = "sf")
co.ee <- st_transform(co, ee)


stations <- stations |> 
  st_transform(st_crs(ee)) |> 
  mutate(spatialCV = case_when(htap_region=="NAM"~1,
                               htap_region=="EUR"~2,
                               htap_region=="EAS"~3),
         randomCV = sample(1:4, nrow(stations), replace=TRUE))


knndm_folds <- CAST::knndm(tpoints=stations, modeldomain=co.ee,
                     clustering = "hierarchical", maxp = 0.8, k=3)

stations$knndmCV <- knndm_folds$clusters

stations_CV <- stations |> 
  select(spatialCV, randomCV, knndmCV) |> 
  pivot_longer(!geometry) |> 
  mutate(name=fct_relevel(name, "randomCV","spatialCV","knndmCV"))


stations_CV <- stations |> 
  select(spatialCV, randomCV, knndmCV) 

stations_plot <- sample_n(stations_CV,1000)


space_right <- -0.3
rCV <- ggplot() +
    geom_sf(data=co.ee$geometry) +
  geom_sf(data=stations_plot,size=1,aes(color=as.factor(randomCV)),shape = "+") +
  scale_color_colorblind("CV fold") +
  theme_minimal(base_size=15) +
    labs(title="random CV") + 
  theme(legend.position = "None",
        plot.title = element_text(face = "bold", size = (15),hjust = 0.5))
kCV <- ggplot() +
    geom_sf(data=co.ee$geometry) +
    geom_sf(data=stations_plot,size=1,aes(color=as.factor(knndmCV)), shape = "+") +
    scale_color_colorblind() +
    labs(title="kNNDM CV") + 
    theme_minimal(base_size = 15) +
    theme(legend.position = "None",
          plot.title = element_text(face = "bold", size = (15),hjust = 0.5))

legend_stations <- get_legend(rCV +
                             guides(color = guide_legend(override.aes = list(size=5))) +
                             theme(legend.position = "right"))
top_row <- plot_grid(rCV,NULL,kCV, NULL, legend_stations, NULL,nrow=1, rel_widths = c(1,0,1, -0.3, 1, space_right))
random_geodist <- plot_geodist(stations_CV,modeldomain=co.ee, unit = "km",
                               sampling="Fibonacci", cvfolds=stations_CV$randomCV)
knndm_geodist <- plot_geodist(stations_CV,modeldomain=co.ee, unit = "km",
                               sampling="Fibonacci", cvfolds=stations_CV$knndmCV)

rgeod <- random_geodist$plot + scale_x_log10() + xlab("")+ theme_minimal(base_size=15) +
  theme(legend.position = "None", aspect.ratio = 0.6) 
kgeod <- knndm_geodist$plot + scale_x_log10() + xlab("")+ theme_minimal(base_size=15) +
  theme(legend.position = "None", aspect.ratio = 0.6) 

legend_geodist <- get_legend(rgeod +
                                guides(color = guide_legend(nrow = 1,override.aes = list(size=5))) +
                                theme(legend.position = "right"))
bottom_row <- plot_grid(rgeod,NULL,kgeod, NULL, legend_geodist, NULL,nrow=1, rel_widths = c(1,0,1, -0.3, 1, space_right))
both_pl <- plot_grid(NULL, top_row, NULL, bottom_row, nrow=4, rel_heights = c(-0.3,1,-0.2,0.5))
ggsave("CV_dist.pdf",both_pl)


(o3_map <- ggplot() +
    geom_sf(data=co.ee) +
    geom_sf(data=stations, size=0.5, aes(color=o3_average_values)) +
    scale_color_viridis_c(expression(paste(O[3], " [ppb]")), limits=c(15,40)) +
    theme_minimal())

ggsave("station_map.pdf", o3_map)

library(doParallel)
cl <- makeCluster(3)
registerDoParallel(cl)

trainDat <- st_drop_geometry(stations)
ctrl <- trainControl(method="cv",index = knndm_folds$clusters)
mod <- ffs(trainDat[,predictors], trainDat[,response], globalval=TRUE, trControl = ctrl)


