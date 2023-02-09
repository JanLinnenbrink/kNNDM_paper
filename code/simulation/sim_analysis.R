#-----------------------------------------------------------#
#             Simulation analysis 2: virtual species        #
#-----------------------------------------------------------#
#.libPaths("/home/j/jlinnenb/r_packages/")
library("parallel")
library("doParallel")
library("pbapply")
setwd("C:/git/kNNDM_paper/")

# Load utils, functions, and define number of iterations
source("code/simulation/sim_functions.R")
nsim <- 100
pboptions(type = "timer")


# Read data
spoly <- st_read(dsn="data/simulation/species_vdata.gpkg", layer="sampling_polygon")
wclim <- rast("data/simulation/species_stack.grd")
wgrid <- st_read(dsn="data/simulation/species_vdata.gpkg", layer="landscape_grid")

# Prepare parallelization
print(paste0("Process started with ", 11, " cores."))
cl <- makeCluster(detectCores()-1)
registerDoParallel(cl)

# Launch simulation
set.seed(1234)
sims <- pbreplicate(nsim, sim_species(wgrid, wclim, spoly), simplify=FALSE)

# We're done
stopCluster(cl)
rm("cl")
write_csv(do.call(rbind, sims), "results/simulation/sim_res.csv")
