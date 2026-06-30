library(ks)
library(terra)
library(kuenm2)


# 1.1 Read occurrences
data <-read.csv("D:/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020.csv") # unique (species, lat, lon, month, year) records
data <- subset (data, genus_name_pub == "Atta")
data <- subset(data,select = c("valid_species_name","dec_long","dec_lat"))
data$dec_long <- as.numeric(data$dec_long)
data$dec_lat <- as.numeric(data$dec_lat)


# env
if(!exists("wrclim")) {
  wcPath <- paste0("ClimateData/climate/wc2.1_2.5m")
  wcFiles <- list.files(wcPath, pattern = "\\.tif$", full.names = T)
  wrclim <- rast(wcFiles)
}
plot(wrclim)
# define Holarctic extent
Ne_ext <- c(-120, -20, -60, 40)

# crop wrclim data for Holarctic
base_raster <- crop(wrclim[[1]], ext(Ne_ext))
plot(base_raster)

# points
occ_data<- initial_cleaning(data, "valid_species_name", "dec_long", "dec_lat",
                            other_columns = NULL, keep_all_columns = TRUE,
                            sort_columns = TRUE, remove_na = TRUE, remove_empty = TRUE,
                            remove_duplicates = FALSE, by_decimal_precision = TRUE,
                            decimal_precision = 0, longitude_precision = NULL,
                            latitude_precision = NULL)
data <- occ_data[,c(2,3)] # only lat, lon
points <- vect(data, geom = c("dec_long", "dec_lat"), crs = crs(base_raster),
               keepgeom = T)

points <- geom(points)[, c("x", "y")]

#
active_cell_coords_projected <- crds(base_raster, na.rm = TRUE)

kde_result_vec <- kde(
  x = points,
   #H = Hpi(point_matrix_projected), # Optional: calculate optimal bandwidth in meters
  eval.points = active_cell_coords_projected
)

density_values_vec <- kde_result_vec$estimate

base_raster[!is.na(values(base_raster))] <- density_values_vec
base_raster <- (base_raster-minmax(base_raster)[1])/minmax(base_raster)[2]
# base_raster[base_raster == minmax(base_raster)[1]] <- 0
plot(base_raster)
minmax(base_raster)
dir.create("BiasFile")
writeRaster(base_raster, "BiasFile/bias_raster.tif")

#
cell_samp <- terra::as.data.frame(base_raster, na.rm = TRUE, cells = TRUE)

bias_value <- cell_samp[, 2]
cell_samp <- cell_samp$cell

cell_samp1 <- sample(cell_samp, size = 1000, replace = FALSE,
                     prob = bias_value)

bg_var <- terra::extract(x = base_raster, y = cell_samp1, xy = TRUE)
points(bg_var[,1:2])
