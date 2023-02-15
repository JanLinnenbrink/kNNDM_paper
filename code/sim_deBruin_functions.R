# Prepare simulation results as in de Bruin et al.
library(sf)
library(raster)
library(ggpubr)

setwd("C:/0_Msc_Loek/Z_Palma/deBruin_add_nndm/")

infolder  <- "./CVresults"
outfolder <- "./material"

mets <- c("exhaustive", "heteroscedastic", "intensity", "modelbased","random", "spatial", "knndm_sample", "knndm_sample_4")

colnms <- c("method", "variate", "design", "number", "RMSE", "MEC", "time","WS")
outtab <- data.frame(matrix(NA, 0, 8))
names(outtab) <- colnms

for(m in mets){
  p <- file.path(infolder, m)
  f_ins <- list.files(p, glob2rx("???_*.Rdata"))
  for(f_in in f_ins){
    lchar <- nchar(f_in)
    variate <- substr(f_in, 1, 3)
    design <- substr(f_in, 5, lchar-9)
    number <- as.numeric(substr(f_in, lchar-8, lchar-6))
    load(file.path(p, f_in))
    if(m == "modelbased" | m == "heteroscedastic"){
      MEC    <- mean(MECs)
      RMSE   <- mean(RMSEs)
    } else{
      if(length(MEC) > 1){
        MEC  <- mean(MEC)
        RMSE <- mean(RMSE)
      }
      
    }
    if (m == "exhaustive") WS=time = 0 
    newrow <- data.frame(method = m, variate = variate, design = design,
                         number = number, RMSE = RMSE, MEC = MEC, time=time,WS=WS)
    outtab <- rbind(outtab, newrow)
  }
}

outtab$methodID <- with(outtab, ifelse(method=="exhaustive",0,1))


# relative RMSE & MEC
outtab$rRMSE <- NA
outtab$rMEC  <- NA

# some rows missing in exhaustive method for OCS data
numbers <- outtab[outtab$variate == "OCS" & outtab$method=="exhaustive" & outtab$design=="clusterGapped",]$number
outtab <- outtab[outtab$number %in% numbers, ]

for(variate in c("AGB", "OCS")){
  for(design in unique(outtab$design)){
    for(number in numbers){
      idx1 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$number == number)
      idx2 <- which(outtab$design == design & outtab$variate == variate & 
                      outtab$methodID == 0 & outtab$number == number)
      
      outtab$rRMSE[idx1] <- 100 * (outtab$RMSE[idx1] - outtab$RMSE[idx2])/
        outtab$RMSE[idx2]
      outtab$rMEC[idx1] <- 100 * (outtab$MEC[idx1] - outtab$MEC[idx2])/
        outtab$MEC[idx2]
    }
  }
}

write.csv(outtab, file.path(outfolder, "outtab.csv"))

