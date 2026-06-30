library(terra)
presente <- rast("ENMs/Atta sexdens/Maxnet/Final_models_current/General_consensus_bin.tif")
futuro <- rast( "ENMs/Atta sexdens/Maxnet/Projections/Future/2081-2100/ssp245/ACCESS-CM2/General_consensus_bin.tif")
M <- vect("ENMs/Atta sexdens/M_simulation/M_simulation_2_SD_1000_events_4_rep10/accessible_area_M.shp")
plot(M)
presente <- crop(presente, M, mask = T) 
futuro <- crop(futuro, M, mask = T) 
plot(futuro)
plot(presente)
resultado <-  (presente) + (futuro*2)
cls <- data.frame(id = 0:3, 
                  categoria = c("Ausente em ambos", "Trailing edge", "Leading edge", "Estabiliade"))
levels(resultado) <- cls

# 4. Plotar
plot(resultado, col=c("lightgray", "red", "green", "yellow"), plg=list(title="Classes"))


     