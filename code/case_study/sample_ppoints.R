# sampling of global prediction points.
# Prediction raster needed, can be obtained from https://isp.uv.es/code/try.html

setwd("C:/git/kNNDM_paper/")

library("terra")

sla_pred <- rast("data/SLA_1km_v1.tif") |> 
  terra::aggregate(10)
sla_pred[sla_pred<=0] <- NA
ppoints <- spatSample(sla_pred, 10000, "regular", values=FALSE,
                      na.rm=TRUE,as.points=TRUE) |> 
  sample(1000) 

plot(sla_pred)
plot(ppoints, add=TRUE)

writeVector(ppoints,"data/case_study/ppoints.gpkg", overwrite=TRUE)
