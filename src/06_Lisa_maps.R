# !diagnostics off
source('src/01_preprocessed_data_loading.R')
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
library(rgdal)
library(stringr)
library(magrittr)
library(compositions)
library(tidyr)
library(fda)
library(fields)


  pic_path <- paste0("results/pics/", get_filename(), "/")

  create_folder_pics(get_filename())


italy_fortified_ <- italy_fortified %>% left_join(data %>% select(PROCOM,AGMAX_50, DZCOM))

p <- ggplot(italy_fortified_) + # no backgroundcolor
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = AGMAX_50
  ), alpha = 0.9 
  ,size = 0.0015, color = 'black'
  ) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  #scale_fill_viridis_c(begin = 0.25, end = 0.85, option = 'inferno')+
  scale_fill_gradientn('ag[max]', colors =tim.colors(63)) +
  geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
  geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
  geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
  scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                      , "ag[max] > 0.15" = "dashed",
                                      'ag[max] > 0.25'= 'solid')
  ) +
  guides(linetype = guide_legend(order = 1))
p

save_MoI('MoI_a(g)_big_1')

data <- data %>% arrange(PROCOM)


italy_ <- italy %>% subset(., DZCOM != 'Mappano')
italy_ <- italy_[order(italy_$PROCOM),]

italy_@data %<>% left_join(data %>%
                             select(c('PROCOM', IVSM)))
italy_nb <- poly2nb(italy_, queen = TRUE)
italy_lw <- nb2listw(italy_nb, style = "B", zero.policy = T)  



W  <- as(italy_lw, "symmetricMatrix")
W  <- as.matrix(W)
W  <-(W/rowSums(W))
W[which(is.na(W))] <- 0

nr_ <- nrow(data)
neigh <- rowSums(W >0)

library(future)
library(tictoc)
future::plan(multiprocess)
tic()
perm_list <- furrr::future_map(.x = 1:nr_, .f = function(x) replicate(10000, sample(setdiff(1:nr_, x), neigh[x])))
toc()                  


#saveRDS(perm_list, 'perm_list1.RDS', compress = F)




#perm_list <- readRDS('perm_list1.RDS')
i <- 1

tic()
qt_IVSM <- sapply(1:7953, function(i) { 
  qt_ <- 1/2
  if (is.null(dim(perm_list[[i]]))== F) {
    if (dim(perm_list[[i]]) == 1)    perm_list[[i]] <- matrix(perm_list, nrow = 1)
    qt_ <- ecdf( apply(perm_list[[i]], 2, function(x) italy_$IVSM[x]) %>% colMeans())((W[i,]*italy_$IVSM) %>% sum())
  }
  return(qt_)
})
toc()


morans_I_IVSM <- ape::Moran.I(italy_$IVSM, W, na.rm = T)

italy_$IVSM_qt <- 'Not Significant'
italy_$IVSM_qt[qt_IVSM > 0.975] <- 'High'
italy_$IVSM_qt[qt_IVSM < 0.025] <- 'Low'

italy_$IVSM_cl <- 'Not Significant'
italy_$IVSM_cl[italy_$IVSM >= mean(italy_$IVSM)] <- 'High'
italy_$IVSM_cl[italy_$IVSM < mean(italy_$IVSM)] <- 'Low'

italy_$IVSM_LISA <- paste0(italy_$IVSM_cl,'-',italy_$IVSM_qt)
italy_$IVSM_LISA[grepl('Not Sign', italy_$IVSM_LISA) == T] <- 'Not Significant'
italy_$IVSM_LISA <- factor(italy_$IVSM_LISA, levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))
italy_$IVSM_LISA[1]

italy_$AG <- data$AGMAX_50
italy_$AG_gr <- ifelse(italy_$AG >= 0.15, 'Severe', 'Mild')


italy_$IVSM_LISA_2 <- as.character(italy_$IVSM_LISA)
italy_$IVSM_LISA_2[grepl('High-L|Low-H', italy_$IVSM_LISA_2) == T] <- 'Not Significant'
italy_$IVSM_LISA_2[grepl('High-H', italy_$IVSM_LISA_2) == T] <- 'High'
italy_$IVSM_LISA_2[grepl('Low-L', italy_$IVSM_LISA_2) == T] <- 'Low'

italy_$IVSM_LISA_AG <- paste0(italy_$IVSM_LISA_2,'-',italy_$AG_gr)
italy_$IVSM_LISA_AG[grepl('Not Sign', italy_$IVSM_LISA_AG) == T] <- 'NA'

italy_$IVSM_LISA_AG <- factor(italy_$IVSM_LISA_AG,  levels = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))


italy_fortified_ <- italy_fortified %>% left_join(italy_@data)


# LISA: IVSM

p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, IVSM_LISA)) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = IVSM_LISA
  ), alpha = 0.9, size = 0.0015, color = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val,na.translate = F) +
  geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
  geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
  geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
  scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                      , "ag[max] > 0.15" = "dashed",
                                      'ag[max] > 0.25'= 'solid')
  ) +
  guides(linetype = guide_legend(order = 1))
print(p)

save_MoI(name = paste0("LISA Map_IVSM and contours"))


# LISA: IVSM stratified by seismic hazard

col_val <- c("red2","yellow1", "dodgerblue2", "springgreen3")
p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, IVSM_LISA_AG) ) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = IVSM_LISA_AG
  ), alpha = 0.9, size = 0.0015, col = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val, na.translate = T, na.value = 'gray', breaks = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))
print(p)
save_MoI(name = paste0("LISA Map_IVSM_stratified by seismic hazard"))

####







### LISA + Stratified

####Var_perc
italy_@data %<>% left_join(data %>% select(PROCOM, VAR_PERC, ETA_Q3))

tic()
qt_VAR_PERC <- sapply(1:7953, function(i) { 
  qt_ <- 1/2
  if (is.null(dim(perm_list[[i]]))== F) {
    if (dim(perm_list[[i]]) == 1)    perm_list[[i]] <- matrix(perm_list, nrow = 1)
    qt_ <- ecdf( apply(perm_list[[i]], 2, function(x) italy_$VAR_PERC[x]) %>% colMeans())((W[i,]*italy_$VAR_PERC) %>% sum())
  }
  return(qt_)
})
toc()


morans_I_VAR_PERC <- ape::Moran.I(italy_$VAR_PERC, W, na.rm = T)

italy_$VAR_PERC_qt <- 'Not Significant'
italy_$VAR_PERC_qt[qt_VAR_PERC > 0.975] <- 'High'
italy_$VAR_PERC_qt[qt_VAR_PERC < 0.025] <- 'Low'

italy_$VAR_PERC_cl <- 'Not Significant'
italy_$VAR_PERC_cl[italy_$VAR_PERC >= mean(italy_$VAR_PERC)] <- 'High'
italy_$VAR_PERC_cl[italy_$VAR_PERC < mean(italy_$VAR_PERC)] <- 'Low'

italy_$VAR_PERC_LISA <- paste0(italy_$VAR_PERC_cl,'-',italy_$VAR_PERC_qt)
italy_$VAR_PERC_LISA[grepl('Not Sign', italy_$VAR_PERC_LISA) == T] <- 'Not Significant'
italy_$VAR_PERC_LISA <- factor(italy_$VAR_PERC_LISA, levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))
italy_$VAR_PERC_LISA[1]

italy_$AG <- data$AGMAX_50
italy_$AG_gr <- ifelse(italy_$AG >= 0.15, 'Severe', 'Mild')


italy_$VAR_PERC_LISA_2 <- as.character(italy_$VAR_PERC_LISA)
italy_$VAR_PERC_LISA_2[grepl('High-L|Low-H', italy_$VAR_PERC_LISA_2) == T] <- 'Not Significant'
italy_$VAR_PERC_LISA_2[grepl('High-H', italy_$VAR_PERC_LISA_2) == T] <- 'High'
italy_$VAR_PERC_LISA_2[grepl('Low-L', italy_$VAR_PERC_LISA_2) == T] <- 'Low'

italy_$VAR_PERC_LISA_AG <- paste0(italy_$VAR_PERC_LISA_2,'-',italy_$AG_gr)
italy_$VAR_PERC_LISA_AG[grepl('Not Sign', italy_$VAR_PERC_LISA_AG) == T] <- 'NA'
#italy_$VAR_PERC_LISA_AG[grepl('High-High', italy_$VAR_PERC_LISA_AG) == T] <- 'Not Significant'

italy_$VAR_PERC_LISA_AG <- factor(italy_$VAR_PERC_LISA_AG,  levels = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))


italy_fortified_ <- italy_fortified %>% left_join(italy_@data)



# LISA: VAR_PERC

p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, VAR_PERC_LISA)) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = VAR_PERC_LISA
  ), alpha = 0.9, size = 0.0015, col = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val,na.translate = F) +
  geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
  geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
  geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
  scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                      , "ag[max] > 0.15" = "dashed",
                                      'ag[max] > 0.25'= 'solid')
  ) +
  guides(linetype = guide_legend(order = 1))
print(p)

save_MoI(name = paste0("LISA Map_VAR_PERC and contours"))


# LISA: VAR_PERC stratified by seismic hazard

col_val <- c("red2","yellow1", "dodgerblue2", "springgreen3")
p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, VAR_PERC_LISA_AG) ) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = VAR_PERC_LISA_AG
  ), alpha = 0.9, size = 0.0015, col = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val, na.translate = T, na.value = 'gray', breaks = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))
print(p)
save_MoI(name = paste0("LISA Map_VAR_PERC_stratified by seismic hazard"))


#### ETA_Q3


tic()
qt_ETA_Q3 <- sapply(1:7953, function(i) { 
  qt_ <- 1/2
  if (is.null(dim(perm_list[[i]]))== F) {
    if (dim(perm_list[[i]]) == 1)    perm_list[[i]] <- matrix(perm_list, nrow = 1)
    qt_ <- ecdf( apply(perm_list[[i]], 2, function(x) italy_$ETA_Q3[x]) %>% colMeans())((W[i,]*italy_$ETA_Q3) %>% sum())
  }
  return(qt_)
})
toc()


morans_I_ETA_Q3 <- ape::Moran.I(italy_$ETA_Q3, W, na.rm = T)

italy_$ETA_Q3_qt <- 'Not Significant'
italy_$ETA_Q3_qt[qt_ETA_Q3 > 0.975] <- 'High'
italy_$ETA_Q3_qt[qt_ETA_Q3 < 0.025] <- 'Low'

italy_$ETA_Q3_cl <- 'Not Significant'
italy_$ETA_Q3_cl[italy_$ETA_Q3 >= mean(italy_$ETA_Q3)] <- 'High'
italy_$ETA_Q3_cl[italy_$ETA_Q3 < mean(italy_$ETA_Q3)] <- 'Low'

italy_$ETA_Q3_LISA <- paste0(italy_$ETA_Q3_cl,'-',italy_$ETA_Q3_qt)
italy_$ETA_Q3_LISA[grepl('Not Sign', italy_$ETA_Q3_LISA) == T] <- 'Not Significant'
italy_$ETA_Q3_LISA <- factor(italy_$ETA_Q3_LISA, levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))
italy_$ETA_Q3_LISA[1]

italy_$AG <- data$AGMAX_50
italy_$AG_gr <- ifelse(italy_$AG >= 0.15, 'Severe', 'Mild')


italy_$ETA_Q3_LISA_2 <- as.character(italy_$ETA_Q3_LISA)
italy_$ETA_Q3_LISA_2[grepl('High-L|Low-H', italy_$ETA_Q3_LISA_2) == T] <- 'Not Significant'
italy_$ETA_Q3_LISA_2[grepl('High-H', italy_$ETA_Q3_LISA_2) == T] <- 'High'
italy_$ETA_Q3_LISA_2[grepl('Low-L', italy_$ETA_Q3_LISA_2) == T] <- 'Low'

italy_$ETA_Q3_LISA_AG <- paste0(italy_$ETA_Q3_LISA_2,'-',italy_$AG_gr)
italy_$ETA_Q3_LISA_AG[grepl('Not Sign', italy_$ETA_Q3_LISA_AG) == T] <- 'NA'
#italy_$ETA_Q3_LISA_AG[grepl('High-High', italy_$ETA_Q3_LISA_AG) == T] <- 'Not Significant'

italy_$ETA_Q3_LISA_AG <- factor(italy_$ETA_Q3_LISA_AG,  levels = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))


italy_fortified_ <- italy_fortified %>% left_join(italy_@data)



# LISA: ETA_Q3
p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, ETA_Q3_LISA)) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = ETA_Q3_LISA
  ), alpha = 0.9, size = 0.0015, col='black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val,na.translate = F)+
  geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
  geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
  geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
  scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                      , "ag[max] > 0.15" = "dashed",
                                      'ag[max] > 0.25'= 'solid')
  ) +
  guides(linetype = guide_legend(order = 1))
print(p)

save_MoI(name = paste0("LISA Map_ETA_Q3 and contours"))


# LISA: ETA_Q3 stratified by seismic hazard

col_val <- c("red2","yellow1", "dodgerblue2", "springgreen3")
p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, ETA_Q3_LISA_AG) ) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = ETA_Q3_LISA_AG
  ), alpha = 0.9, size = 0.0015, col = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val, na.translate = T, na.value = 'gray', breaks = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))
print(p)
save_MoI(name = paste0("LISA Map_ETA_Q3_stratified by seismic hazard"))

### BUILDING_PC1

data$Building_PC <-  data %>% select(paste0('E', 8:16)) %>% compositions::clr() %>% princomp(.) %>% .$scores %>% .[,1]

#### 
italy_@data %<>% left_join(data %>% select(PROCOM, Building_PC))


tic()
qt_Building_PC <- sapply(1:7953, function(i) { 
  qt_ <- 1/2
  if (is.null(dim(perm_list[[i]]))== F) {
    if (dim(perm_list[[i]]) == 1)    perm_list[[i]] <- matrix(perm_list, nrow = 1)
    qt_ <- ecdf( apply(perm_list[[i]], 2, function(x) italy_$Building_PC[x]) %>% colMeans())((W[i,]*italy_$Building_PC) %>% sum())
  }
  return(qt_)
})
toc()


morans_I_Building_PC <- ape::Moran.I(italy_$Building_PC, W, na.rm = T)

qt_df <- data.frame(qt_IVSM,
                    qt_ETA_Q3,
                    qt_VAR_PERC,
                    qt_Building_PC)

saveRDS(qt_df, file = 'qt_df.RDS')


italy_$Building_PC_qt <- 'Not Significant'
italy_$Building_PC_qt[qt_Building_PC > 0.975] <- 'High'
italy_$Building_PC_qt[qt_Building_PC < 0.025] <- 'Low'

italy_$Building_PC_cl <- 'Not Significant'
italy_$Building_PC_cl[italy_$Building_PC >= mean(italy_$Building_PC)] <- 'High'
italy_$Building_PC_cl[italy_$Building_PC < mean(italy_$Building_PC)] <- 'Low'

italy_$Building_PC_LISA <- paste0(italy_$Building_PC_cl,'-',italy_$Building_PC_qt)
italy_$Building_PC_LISA[grepl('Not Sign', italy_$Building_PC_LISA) == T] <- 'Not Significant'
italy_$Building_PC_LISA <- factor(italy_$Building_PC_LISA, levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))
italy_$Building_PC_LISA[1]

italy_$AG <- data$AGMAX_50
italy_$AG_gr <- ifelse(italy_$AG >= 0.15, 'Severe', 'Mild')


italy_$Building_PC_LISA_2 <- as.character(italy_$Building_PC_LISA)
italy_$Building_PC_LISA_2[grepl('High-L|Low-H', italy_$Building_PC_LISA_2) == T] <- 'Not Significant'
italy_$Building_PC_LISA_2[grepl('High-H', italy_$Building_PC_LISA_2) == T] <- 'High'
italy_$Building_PC_LISA_2[grepl('Low-L', italy_$Building_PC_LISA_2) == T] <- 'Low'

italy_$Building_PC_LISA_AG <- paste0(italy_$Building_PC_LISA_2,'-',italy_$AG_gr)
italy_$Building_PC_LISA_AG[grepl('Not Sign', italy_$Building_PC_LISA_AG) == T] <- 'NA'

italy_$Building_PC_LISA_AG <- factor(italy_$Building_PC_LISA_AG,  levels = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))

italy_fortified_ <- italy_fortified %>% left_join(italy_@data)

# LISA: Building_PC

col_val <- c("red2","yellow1", "dodgerblue2", "springgreen3", "gray")

p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, Building_PC_LISA)) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = Building_PC_LISA
  ), alpha = 0.9, size = 0.0015, col = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val,na.translate = F) +
  geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
  geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
  geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
  scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                      , "ag[max] > 0.15" = "dashed",
                                      'ag[max] > 0.25'= 'solid')
  ) +
  guides(linetype = guide_legend(order = 1))

print(p)

save_MoI(name = paste0("LISA Map_Building_PC and contours"))

# LISA: Building_PC stratified by seismic hazard

col_val <- c("red2","yellow1", "dodgerblue2", "springgreen3")
p <- ggplot(italy_fortified_ %>% dplyr::select(long, lat, group, Building_PC_LISA_AG) ) +
  theme_minimal() +
  geom_polygon(aes(
    x = long, y = lat,
    group = group,
    fill = Building_PC_LISA_AG
  ), alpha = 0.9, size = 0.0015, col = 'black') +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18),
    plot.title = element_text(size = 18)
  ) +
  xlab("longitude") +
  ylab("latitude") +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2)) +
  coord_map(projection = "lambert", parameters = c(lat0 = 35, lat1 = 49)) +
  scale_fill_manual('',values = col_val, na.translate = T, na.value = 'gray', breaks = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))


print(p)
save_MoI(name = paste0("LISA Map_Building_PC1_stratified by seismic hazard"))
