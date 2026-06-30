library(grinnell)
library(sf)
library(terra)
library(geodata)
library(spThin)


# 1.0 Climate data (WorldClim at 5 arc minutes resolution)

# 1.1 Download WorldClim bioclimatic variables

dir.create("ClimateData")

wrclim <- worldclim_global(var = "bio", res = 2.5, path = "ClimateData")

# read tif files if needed
if(!exists("wrclim")) {
  wcPath <- paste0("ClimateData/climate/wc2.1_2.5m")
  wcFiles <- list.files(wcPath, pattern = "\\.tif$", full.names = T)
  wrclim <- rast(wcFiles)
}


# 1.2 Crop wrclim data for America extent

# define tropical extent
Ne_ext <- c(-120, -20, -60, 40)

# crop wrclim data for Nearctic
Ne_wrclim <- crop(wrclim, ext(Ne_ext))
plot(Ne_wrclim[[1]]) # looks ok

# 1.3 Crop wrclim data for Holarctic extent

# define Holarctic extent
Ho_ext <- c(-170,-20, -60, 84)

# crop wrclim data for Holarctic
Ho_wrclim <- crop(wrclim, ext(Ho_ext))
plot(Ho_wrclim[[1]]) # ok

# 2.0 Occurrence data

# 2.1 Read cleaned and thinned occurrence data. First for one species, then adapt the code to run over loops. We chose Bombus bimaculatus as an example species to work with initially.
data <-read.csv("ENMs/Atta sexdens/occurrences.csv")

# 2.2 Additional thinning - B. bimaculatus (and some other common ones) present thousands of points, a lot of them biased toward some regions with a lot of iNat observations. This issue is resulting in poor simulations of M  that do not recover Florida, where this species actually occur. So, apply a very severe thinning (500 km or 1000 km ) and see what happens.

occ50 <- thin(loc.data = data, lat.col = "dec_lat", long.col = "dec_long",   
              spec.col = "valid_species_name", thin.par = 5, reps = 5,
              locs.thinned.list.return = FALSE, write.files = TRUE,
              max.files = 1, out.dir = "ENMs/Atta sexdens/",
              out.base = "50thinned")

occ50 <-read.csv("ENMs/Atta sexdens/50thinned_thin1.csv")

# XX. just some ideas for the loop: define Nearctic and Holarctic extents, then, for each species, check species' geographical extent. If the species' geographical extent fell within Nearctic, then use variables cropped to Nearctic; else, use Holarctic variables. 

sp_ext <- c(min(data$dec_long), max(data$dec_long),
            min(data$dec_lat), max(data$dec_lat))
sp_ext <- ext(sp_ext) # spatExtent

# 3. Parameters for M simulations
dir.create("ENMs/Atta sexdens/M_simulation")

ks <- seq(from = 1, to = 5, by = 1) # Kernel spreads
dispersal.events <- c(200, 500) # events
max.dispersers <- 4 # numbers of dispersers per event

# list to store all Ms
m_ls <- list()
p = 0 

# Select G variables for M simulations

# If sp_ext is within Nearctic extent, then select Ne_wrclim as G; else, select Ho_wrclim as G
if (xmin(sp_ext) >= xmin(Ne_wrclim) &&
    xmax(sp_ext) <= xmax(Ne_wrclim) &&
    ymin(sp_ext) >= ymin(Ne_wrclim) &&
    ymax(sp_ext) <= ymax(Ne_wrclim)) {
  
  G_vars <- Ne_wrclim
  
} else {
  
  G_vars <- Ho_wrclim
  
}

#Perform the M simulation
# Loop for all parameterizations
for (j in 1:length(ks)){
  for (n in 1:length(dispersal.events)){
    for (d in 1:length(max.dispersers)){
      
      # create output directory
      out_dir <- paste("ENMs/Atta sexdens/M_simulation/M_simulation", 
                       gsub("\\.", "-", ks[j]),
                       "SD", dispersal.events[n], "events",max.dispersers[d], "rep10", sep = "_")
      p = p + 1
      m_ls[[p]] <- M_simulationR(data = occ50, 
                                 current_variables = G_vars,
                                 starting_proportion = 0.5,
                                 max_dispersers = max.dispersers[d],
                                 replicates = 10,
                                 suitability_threshold = 1,
                                 dispersal_events = dispersal.events[n],
                                 kernel_spread = ks[j],
                                 dispersal_kernel = "normal",
                                 write_all_scenarios = T,
                                 set_seed = 21,
                                 output_directory = out_dir, 
                                 overwrite = TRUE)
    }
  }
}

# 1.5 SD, 100 events, resulted the most reliable one

