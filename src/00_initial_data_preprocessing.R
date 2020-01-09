library(tidyverse)
library(sp)
library(openxlsx)
library(rgdal)
library(proj4)
library(raster)
library(geosphere)
library(sf)
library(magrittr)

### Reading and preprocessing MRCI data
#### Preparation: Decimal separator changed from ',' to '.'. NA values identifier '--' changed to '-9999' 
data<- read.xlsx("data/initial/Indicatori_Intero_territorio_nazionale_2018.xlsx",  sheet = 1) %>%
  filter(DZCOM != 'Mappano')

str_vec <- apply(data[,8:397], 2, function(x) sum(is.na(as.numeric(x))))
data[, which(str_vec == 0) + 7] <- lapply(data[, which(str_vec == 0) + 7], as.numeric)
data[data == -9999] <- NA
rm(str_vec)

### Reading and preprocessing Map Of Italy spatial df:

italy <- readOGR('data/initial/map_of_italy_istat/Com01012018_g_WGS84.shp')
italy$COMUNE <- iconv(from = "UTF-8", to = 'ISO-8859-1', italy$COMUNE)



### Merging MRCI and  Map Of Italy spatial df:
italy$COMUNE[italy$PRO_COM_T %in% setdiff( italy$PRO_COM_T,data$PROCOM)]
data$DZCOM[ which(data$PROCOM %in% setdiff( data$PROCOM,italy$PRO_COM_T))]

### Barbaro Mossano = Barbaro Vicentino + Mossano 
### Borgo Veneto =  Megliadino San Fidenzio, Saletto and Santa Margherita d'Adige
### Fiumicello Villa Vicentina = Fiumicello and Villa Vicentina.
### Treppo Ligosullo = Ligosullo and Treppo Carnico.
### Corigliano-Rossano = Corigliano Calabro and Rossano.
### Torre de' Busi code has been updated

com_upd_l <- list(
  list('Barbarano Mossano', c('Barbarano Vicentino', 'Mossano') ),
  list('Borgo Veneto',  c('Megliadino San Fidenzio', 'Saletto',  "Santa Margherita d'Adige")),
  list('Fiumicello Villa Vicentina', c('Fiumicello', 'Villa Vicentina')),
  list('Treppo Ligosullo',c( 'Ligosullo', 'Treppo Carnico')),
  list("Corigliano-Rossano", c("Corigliano Calabro","Rossano")),
  list("Torre de' Busi", "Torre de' Busi")
)

italy$PROCOM <- as.character(italy$PRO_COM_T)

italy_g_t <- italy_g <- italy[!italy$COMUNE %in% unlist(sapply(1:length(com_upd_l), function(i) com_upd_l[[i]][[2]])), ]
for (i in 1:length(com_upd_l)) {
  italy_com <- com_upd_l[[i]][[2]]
  italy_t <- raster::aggregate(rbind(italy[italy$COMUNE %in% italy_com,]))
  id_ <- (italy_g@polygons %>% length()) -1
  italy_g <- spChFIDs(italy_g, as.character(0:id_))
  italy_t <- spChFIDs(italy_t, as.character(id_ + 1))
  italy_g <- maptools::spRbind(italy_g, italy_t)
}


df_for_it_sp_pol <- data.frame(PROCOM = c(as.character(italy_g_t$PRO_COM_T),
                                          sapply(1:length(com_upd_l), function(i) data$PROCOM[data$DZCOM == com_upd_l[[i]][[1]]]))) 



italy_sp <- SpatialPolygonsDataFrame(italy_g, df_for_it_sp_pol, match.ID = F)
italy_sp@data %<>% left_join(data %>% dplyr::select(PROCOM, DZCOM, DZPRO, DZREG))

italy <- italy_sp

rm(italy_sp, italy_g, italy_g_t, df_for_it_sp_pol)


italy <- spTransform(italy, CRS("+proj=longlat +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"), inverse = T)
italy_fortified <- fortify(italy) %>%
  left_join(italy@data %>%
              mutate(id = as.character(0:(nrow(.) - 1))     )
  )


italy_sf <- sf::st_as_sf(italy_fortified, coords = c("long", "lat")) %>% 
  group_by(id, piece) %>% 
  summarize(do_union=FALSE) %>%
  st_cast("POLYGON") %>% 
  ungroup()


id_median <- italy_fortified %>% 
  dplyr::select(id, PROCOM) %>% 
  distinct() %>% 
  left_join(data %>% 
              dplyr::select(AGMAX_50, PROCOM)) %>% 
  filter(AGMAX_50 > 0.05) %>% 
  dplyr::select(id) %>% 
  pull()


id_moderate <- italy_fortified %>% 
  dplyr::select(id, PROCOM) %>% 
  distinct() %>% 
  left_join(data %>% 
              dplyr::select(AGMAX_50, PROCOM)) %>% 
  filter(AGMAX_50 > 0.15) %>% 
  dplyr::select(id) %>% 
  pull()


id_high <- italy_fortified %>% 
  dplyr::select(id, PROCOM) %>% 
  distinct() %>% 
  left_join(data %>% 
              dplyr::select(AGMAX_50, PROCOM)) %>% 
  filter(AGMAX_50 > 0.25) %>% 
  dplyr::select(id) %>% 
  pull()


italy_ag_median_borders <- st_union(italy_sf[italy_sf$id %in% id_median,]) %>% as(., 'Spatial') %>% fortify()
italy_ag_moderate_borders <- st_union(italy_sf[italy_sf$id %in% id_moderate,]) %>% as(., 'Spatial') %>% fortify()
italy_ag_high_borders <- st_union(italy_sf[italy_sf$id %in% id_high,])%>% as(., 'Spatial')  %>% fortify()


italy_fortified <- italy_fortified %>% left_join(data %>% dplyr::select(PROCOM, DZREG, DZPRO, CODPRO))

saveRDS(list(data, italy, italy_fortified, italy_ag_moderate_borders,
             italy_ag_high_borders,
             italy_ag_median_borders), file ='data/created/all_data_2018.RDS')


rm(list = ls())