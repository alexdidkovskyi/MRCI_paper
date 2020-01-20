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


### Map of Italy colored acc. to ag[max]

p <- MoI_p(data_ = italy_fortified_,
           fill_ = 'AGMAX_50',
           m_col = 'black',
           sc_fill_value = tim.colors(63),
           sc_fill_title = 'ag[max]',
           zone_b = T, 
           col_t = 'g_n')
p

save_MoI('MoI_a(g)')




### Lisa Maps. Spatal data preparation
data <- data %>% arrange(PROCOM)
data$Building_PC <-  data %>% select(paste0('E', 8:16)) %>% compositions::clr() %>% princomp(.) %>% .$scores %>% .[,1]


italy_ <- italy %>% subset(., DZCOM != 'Mappano')
italy_ <- italy_[order(italy_$PROCOM),]

italy_@data %<>% left_join(data %>%
                             select(c('PROCOM', IVSM,VAR_PERC, ETA_Q3, Building_PC)))
italy_$AG <- data$AGMAX_50
italy_$AG_gr <- ifelse(italy_$AG >= 0.15, 'Severe', 'Mild')

italy_nb <- poly2nb(italy_, queen = TRUE)
italy_lw <- nb2listw(italy_nb, style = "B", zero.policy = T)  



W  <- as(italy_lw, "symmetricMatrix")
W  <- as.matrix(W)
W  <-(W/rowSums(W))
W[which(is.na(W))] <- 0

nr_ <- nrow(data)
neigh <- rowSums(W >0)



### Municipalities permutation
library(future)
library(tictoc)
future::plan(multiprocess)
tic()
perm_list <- furrr::future_map(.x = 1:nr_, .f = function(x) replicate(100, sample(setdiff(1:nr_, x), neigh[x])))
toc()                  




for (var in c('IVSM', 'Building_PC','VAR_PERC','ETA_Q3')) {
  temp_v <- italy_@data %>% select(var) %>% pull()
  tic()
  temp <- sapply(1:7953, function(i) { 
    qt_ <- 1/2
    if (is.null(dim(perm_list[[i]]))== T) {
      qt_ <- ecdf(temp_v)((W[i,]*temp_v) %>% sum()) 
    } else {
      
      qt_ <- ecdf( apply(perm_list[[i]], 2, function(x) temp_v[x]) %>% colMeans())((W[i,]*temp_v) %>% sum())
      
    }
    return(qt_)
  })
  toc()
  assign(paste0('qt_',var), temp)
  
  ### Moran's I
  assign(paste0('morans_I_', var), ape::Moran.I(temp_v, W, na.rm = T))
  
  ### Create Lisa map
  italy_$temp_qt <- 'Not Significant'
  italy_$temp_qt[qt_IVSM > 0.975] <- 'High'
  italy_$temp_qt[qt_IVSM < 0.025] <- 'Low'
  
  italy_$temp_cl <- 'Not Significant'
  italy_$temp_cl[ temp_v >= mean( temp_v)] <- 'High'
  italy_$temp_cl[ temp_v < mean( temp_v)] <- 'Low'
  
  italy_$temp_LISA <- paste0(italy_$temp_cl,'-',italy_$temp_qt)
  italy_$temp_LISA[grepl('Not Sign', italy_$temp_LISA) == T] <- 'Not Significant'
  italy_$temp_LISA <- factor(italy_$temp_LISA, levels = c("High-High", "High-Low", "Low-Low", "Low-High", "Not Significant"))
   
  
  italy_$temp_LISA_2 <- as.character(italy_$temp_LISA)
  italy_$temp_LISA_2[grepl('High-L|Low-H', italy_$temp_LISA_2) == T] <- 'Not Significant'
  italy_$temp_LISA_2[grepl('High-H', italy_$temp_LISA_2) == T] <- 'High'
  italy_$temp_LISA_2[grepl('Low-L', italy_$temp_LISA_2) == T] <- 'Low'
  
  italy_$temp_LISA_AG <- paste0(italy_$temp_LISA_2,'-',italy_$AG_gr)
  italy_$temp_LISA_AG[grepl('Not Sign', italy_$temp_LISA_AG) == T] <- 'NA'
  
  italy_$temp_LISA_AG <- factor(italy_$temp_LISA_AG,  levels = c("High-Severe", "High-Mild", "Low-Mild", "Low-Severe"))
  italy_@data %<>% dplyr::rename(!!paste0(var,'_LISA') := temp_LISA,
                          !!paste0(var, '_LISA_2') := temp_LISA_2,
                          !!paste0(var,'_LISA_AG') := temp_LISA_AG,
                          !!paste0(var,'_qt') := temp_qt,
                          !!paste0(var, '_cl') := temp_cl)
  

  italy_fortified_ <- italy_fortified %>% 
    left_join(italy_@data %>% 
                select(PROCOM, DZCOM, DZPRO, DZREG, 
                       paste0(var,"_LISA",c('','_AG'))
                       )
              )
  
  
  
  ### LISA: IVSM
   p <- MoI_p(data_ = italy_fortified_ %>% dplyr::select(long, lat, group, paste0(var,'_LISA')),
             fill_ = !!paste0(var,'_LISA'),
             m_col = 'black', 
             sc_fill_title = '', 
             sc_fill_value = c("red2","yellow1", "dodgerblue2", "springgreen3", "gray"),
             zone_b = T, col_t = 'm')
  
  print(p)
  save_MoI(name = paste0("LISA Map_",var," and contours"))
  

  ### LISA: IVSM stratified by seismic hazard
  
  p <- MoI_p(data_ = italy_fortified_ %>% dplyr::select(long, lat, group, paste0(var,'_LISA_AG')),
             fill_ = !!paste0(var,'_LISA_AG'),
             m_col = 'black', 
             sc_fill_title = '', 
             sc_fill_value = c("red2","yellow1", "dodgerblue2", "springgreen3", "gray"),
             zone_b = F, col_t = 'm')
  
  print(p)
  
  save_MoI(name = paste0("LISA Map_", var,"_stratified by seismic hazard"))
  
}

