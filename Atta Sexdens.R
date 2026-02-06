library(kuenm2)
library(terra)
library(geodata)


#tratamento da planilha
data <-read.csv("D:/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020/GABI_Data_Release1.0_18012020.csv") 
data <- subset (data, valid_species_name == "Atta.sexdens")
data <- subset(data,select = c("valid_species_name","dec_long","dec_lat"))
data$dec_long <- as.numeric(data$dec_long)
data$dec_lat <- as.numeric(data$dec_lat)
#carregar as variaveis climaticas
bio01 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif")
bio02 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_2.tif")
bio03 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_3.tif")
bio04 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_4.tif")
bio05 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_5.tif")
bio06 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_6.tif")
bio07 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_7.tif")
bio08 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_8.tif")
bio09 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_9.tif")
bio10 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_10.tif")
bio11 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_11.tif")
bio12 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_12.tif")
bio13 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_13.tif")
bio14 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_14.tif")
bio15 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_15.tif")
bio16 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_16.tif")
bio17 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_17.tif")
bio18 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_18.tif")
bio19 <- rast("D:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_19.tif")
var1 <- c(bio01,bio02,bio03,bio04,bio05,bio06,bio07,bio10,bio11,bio12,bio13,bio14,bio15,bio16,bio17)
names(var1)
numeros <- setdiff(1:19, c(8, 9, 18, 19))

# Cria os nomes: "bio_1", "bio_2" ... "bio_7", "bio_10" ... "bio_17"
new_names <- paste0("bio_", numeros)
names(var1) <- new_names
print(names(var1))

occ_data<- initial_cleaning(data, "valid_species_name", "dec_long", "dec_lat",
                            other_columns = NULL, keep_all_columns = TRUE,
                            sort_columns = TRUE, remove_na = TRUE, remove_empty = TRUE,
                            remove_duplicates = TRUE, by_decimal_precision = TRUE,
                            decimal_precision = 0, longitude_precision = NULL,
                            latitude_precision = NULL)

Ne_ext <- c(-120, -20, -60, 40)

# crop wrclim data for Holarctic
var_crop <- crop(var1, ext(Ne_ext))
writeRaster(var_crop, "varcrop.tif", overwrite=TRUE)
G_curr <- crop(var1,Ne_ext)
dir.create("ENMs/Atta sexdens/Maxnet", recursive = TRUE)
dir.create("ENMs/Atta sexdens/PCA_Mvars", recursive = TRUE)
biasfile = rast("BiasFile/bias_raster.tif")

bias_file <- resample(x = biasfile, 
                               y = var_crop,
                               method = "bilinear")


sp <- prepare_data(algorithm = "maxnet", occ = occ_data , x = "dec_long", y = "dec_lat", raster_variables = var_crop, species = "Atta.sexdens",
                   n_background = 400, features = c("l", "q", "p", "lq", "lqp"),
                   r_multiplier = c(0.1, 0.25, 0.5, 0.75, 1, 2,3,4,5), partition_method = "kfolds",
                   n_partitions = 4, train_proportion = 0.7,
                   seed = 42,
                   do_pca = TRUE,
                   min_explained = 1,
                   min_number = 4,
                   center = TRUE,
                   scale = TRUE,
                   write_pca = TRUE,
                   bias_file = bias_file,
                   bias_effect = "direct",
                   pca_directory = "ENMs/Atta sexdens/PCA_Mvars",
                   write_file = TRUE,
                   file_name = "ENMs/Atta sexdens/Maxnet/data_preparation")
print(sp)



calibrado <- calibration(data = sp, error_considered = 5,
                         remove_concave = TRUE,
                         proc_for_all = FALSE, 
                         omission_rate = 5, 
                         delta_aic = 2,
                         allow_tolerance = TRUE, 
                         tolerance = 0.01,
                         addsamplestobackground = TRUE, 
                         write_summary = T, 
                         output_directory = "ENMs/Atta sexdens/Maxnet/calibration",
                         skip_existing_models = FALSE, 
                         return_all_results = TRUE,
                         parallel = TRUE, 
                         ncores = 10, 
                         progress_bar = TRUE,
                         verbose = TRUE)
print(calibrado)
partition_response_curves(calibrado, 256, n = 100,
                          averages_from = "pr_bg", col = "darkblue",
                          ylim = NULL, las = 1, parallel = FALSE,
                          ncores = NULL)
encaixado <- fit_selected(calibration_results = calibrado,
                          n_replicates = 10, 
                          write_models = F,
                          parallel = TRUE, 
                          ncores = 6,
                          progress_bar = TRUE, 
                          verbose = TRUE, 
                          seed = 42)

selected_table <- encaixado$selected_models
selected_table$species <- "Atta sexdens"
dir.create("ENMs/Atta sexdens/Maxnet/maxnet_selected_models")
write.csv(selected_table, 
          "ENMs/Atta sexdens/Maxnet/maxnet_selected_models/selected_models_eval.csv", 
          row.names = FALSE)

predict_present <- predict_selected(models = encaixado, 
                                    mask = NULL, 
                                    write_files = TRUE,
                                    out_dir = "ENMs/Atta sexdens/Maxnet/Final_models_current",
                                    consensus_per_model = TRUE, 
                                    consensus_general = TRUE,
                                    consensus = c("median", "range", "mean", "stdev"),
                                    extrapolation_type = "E",
                                    type = "cloglog", 
                                    overwrite = TRUE, 
                                    progress_bar = TRUE,
                                    new_variables = G_curr )

saveRDS(predict_present, "ENMs/Atta sexdens/Maxnet/predict_present.rds")
predict_present <- readRDS("ENMs/Atta sexdens/Maxnet/predict_present.rds")

# 7.0 Response curves

dir.create("ENMs/Atta sexdens/Maxnet/Response_curves")

# 7.1  mean curve
resp_curves <- all_response_curves(
  models = encaixado,
  predictors = G_curr,
  out_dir = "ENMs/Atta sexdens/Maxnet/Response_curves"
)

par(mfrow = c(1,2)) #Set grid of plot
response_curve(models = encaixado, variable = "PC1") # quite nice curve
response_curve(models = encaixado, variable = "PC2") # don't know why it's not showing

geodata_dir <- file.path("ENMs/Future_worldclim")
dir.create(geodata_dir)
#Define GCMs, SSPs and time periods
gcms <- c("ACCESS-CM2", "MIROC6")
ssps <- c("126","245","585")
periods <- c("2081-2100")
#Create a grid of combination of periods, ssps and gcms
g <- expand.grid("period" = periods, "ssps" = ssps, "gcms" = gcms)
g #Each line is a specific scenario for future
#Loop to download variables for each scenario
lapply(1:nrow(g), function(i){
  cmip6_world(model = g$gcms[i], 
              ssp = g$ssps[i], 
              time = g$period[i], 
              var = "bioc", 
              res =  2.5, path = geodata_dir)})
list.files(geodata_dir, recursive = TRUE)
in_dir <- file.path(geodata_dir,"climate")
out_dir_future <- file.path("ENMs/Future_worldclim")
organize_future_worldclim(input_dir = in_dir, #Path to the raw variables from WorldClim
                          output_dir = out_dir_future,
                          variables = c("bio_1","bio_2","bio_3","bio_4","bio_5","bio_6","bio_7","bio_10","bio_11","bio_12","bio_13","bio_14","bio_15","bio_16","bio_17"),
                          name_format = "bio_", mask = m_ext,overwrite = TRUE) 


projetado <- prepare_projection(models = encaixado,
                                present_dir = "ENMs/Var", #Directory with present-day variables
                                past_dir = NULL, #NULL because we won't project to the past
                                past_period = NULL, #NULL because we won't project to the past
                                past_gcm = NULL, #NULL because we won't project to the past
                                future_dir = out_dir_future, #Directory with future variables
                                future_period = c("2081-2100"),
                                future_pscen = c("ssp126","ssp245", "ssp585"),
                                future_gcm = c("ACCESS-CM2", "MIROC6"))

out.dir <- "ENMs/Atta sexdens/Maxnet/Projections"
dir.create(out.dir)
projetado_selecionado <- project_selected(models = encaixado, 
                                          projection_data = projetado,
                                          out_dir = out.dir, 
                                          progress_bar = TRUE, overwrite = TRUE)
print(projetado_selecionado)

p_median <- import_projections(projection = projetado_selecionado, consensus = "median")
#Plot all scenarios
plot(p_median, cex.main = 0.8)
p_mean <- import_projections(projection = projetado_selecionado, consensus = "mean")
plot(p_mean, cex.main = 0.8)

# 9.0 Binarization

# 9.1 Current projection

# Open General consensus
GenCons_curr <- rast("ENMs/Atta sexdens/Maxnet/Final_models_current/General_consensus.tif")

# Extract suitability values
val_training_presence <- extract(GenCons_curr$median, occ_data[2:3])

# Sort the minimum training value
train_pres <- val_training_presence$median
train_pres_sort <- sort(train_pres)

# Select the 5 percent minimum training presence
threshold <- 5
thres_min_train <- train_pres_sort[ceiling(length(occ_data[, 1]) * threshold / 100) + 1]

# Binarization
curr_bin <- (GenCons_curr$median>= thres_min_train) * 1

plot(curr_bin)
writeRaster(curr_bin, "ENMs/Atta sexdens/Maxnet/Final_models_current/General_consensus_bin.tif",overwrite=TRUE)

# 9.2 Future projection

# Open General consensus
GenCons_fut <- rast("ENMs/Atta sexdens/Maxnet/Projections/Future/2081-2100/ssp245/ACCESS-CM2/General_consensus.tif")

# Binarization
fut_bin <- (GenCons_fut$median>= thres_min_train) * 1

plot(fut_bin)
writeRaster(fut_bin, "ENMs/Atta sexdens/Maxnet/Projections/Future/2081-2100/ssp245/ACCESS-CM2/General_consensus_bin.tif",overwrite=TRUE)


# 10.0 MOP 

mop.dir <- "ENMs/Atta sexdens/Maxnet/MOP"
dir.create(mop.dir)

mop <- projection_mop(
  data = encaixado, 
  projection_data = projetado,
  na_in_range = TRUE,
  out_dir = mop.dir,
  type = "basic",
  parallel = FALSE,  
  progress_bar = TRUE,
  overwrite = TRUE)

# PLOT
plot(rast("ENMs/Atta sexdens/Maxnet/MOP/Present/Present_mopbasic.tif"))
plot(rast("ENMs/Atta sexdens/Maxnet/MOP/Future/2081-2100/ssp245/ACCESS-CM2_mopbasic.tif"))
