

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
library(votesys)
library(plotly)



w_ <- c(1/6, 1/6,1/6,1/2)
get_copeland <- function(df_, w_ = NULL){
  if (is.null(w_)) w_  <- rep( 1/nrow(df_),nrow(df_))
  sum_m <- initial_m <- matrix(0, nrow = ncol(df_), ncol = ncol(df_))
  initial_m[lower.tri(initial_m)] <- 1
  initial_m[upper.tri(initial_m)] <- -1
  
  for (i in 1:nrow(df_)){
    initial_m_ <- initial_m
    colnames(initial_m_) <- rev(df_[i,])
    rownames(initial_m_) <- rev(df_[i,])
    initial_m_ <- initial_m_[order(rownames(initial_m_)),order(colnames(initial_m_))]
    sum_m <- sum_m + initial_m_*w_[i]
  }
  rm(initial_m_)
  id_l_0 <- sum_m < 0
  id_m_0 <- sum_m > 0
  
  sum_m[id_l_0] <- 1
  sum_m[id_m_0] <- -1
  
  return(rowSums(sum_m))
  gc()
  gc()
}


  pic_path <- paste0("results/pics/", get_filename(), "/")

  create_folder_pics(get_filename())

  
  
  
  
  
  data$Growth_rate <- (data$POP_2018-data$POP_2011)/data$POP_2011
  sc_rot <- 1
  data$Building_PC <-  data %>% select(paste0('E', 8:16)) %>% compositions::clr() %>% princomp(.) %>% .$scores %>% .[,1]*sc_rot
  
  
  sd_ <- data %>% select(paste0('E', 8:16)) %>% compositions::clr() %>% princomp(.) %>%.$sdev %>% .[1]
  
  
  data$IDEM[is.na(data$IDEM)] <- mean(data$IDEM, na.rm = T)
  
  data_temp <- data %>% select(PROCOM, Growth_rate, AGMAX_50, IVSM, Building_PC, DZCOM, POP_2011, IDEM)
  
  rank_df <- 
    rbind(data_temp %>% select(PROCOM, Building_PC) %>% mutate(Building_PC = rank(Building_PC)) %>% arrange(Building_PC) %>% select(PROCOM) %>% pull(),
          data_temp %>% select(PROCOM, Growth_rate, POP_2011) %>% mutate(Building_PC = rank(-Growth_rate)) %>% arrange(Building_PC, -POP_2011) %>% select(PROCOM) %>% pull(),
          data_temp %>% select(PROCOM, IVSM) %>% mutate(Building_PC = rank(IVSM)) %>% arrange(Building_PC) %>% select(PROCOM) %>% pull(),
          data_temp %>% select(PROCOM, IDEM) %>% mutate(Building_PC = rank(IDEM)) %>% arrange(Building_PC) %>% select(PROCOM) %>% pull()
    )
  
  data %>% select(Building_PC, Growth_rate, IVSM, IDEM, PROCOM) %>% filter(PROCOM %in% c('001001','001002','001003','001004'))
  
  cop_list <- list()
  j<- list(c(1:3), c(1:4), c(1:4))
  w__ <- list(NULL, NULL, c(1/6,1/6,1/6,1/2))
  for (i in 1:3){
    print(i)
    cop_list[[i]] <- get_copeland(df_ = rank_df[j[[i]],], w_ = w__[[i]])
  }
  
  
  cop_df <- do.call('rbind',cop_list)%>% t() %>% data.frame() %>% set_colnames(c('x','y','z')) %>% mutate(PROCOM =  sort(rev(rank_df[1,]))) %>%
    left_join(data_temp %>% select(PROCOM, AGMAX_50))
  
  
  
  
  cop_df  %>%
    mutate(ag_cl = (cut(AGMAX_50, c(0, 0.05, 0.15, 0.25, 0.3)))) %>% ggplot(., aes(x = WO_IDEM, col = ag_cl)) +geom_density(size = 2) +
    theme_bw()+
    xlab('Cop. Score') +
    scale_color_manual('AG Class',values = tim.colors(5)[1:4])
  
  quantile(cop_df$WO_IDEM, 0.75)
  save_pic('cop_score_distribution')
  
  
  ggplot(cop_df %>% mutate(AG_cl = cut(AGMAX_50, c(0, 0.05, 0.15, 0.25, 0.35))), aes(x = AG_cl, y = WO_IDEM)) +geom_boxplot() +
    xlab('AG') +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18)
    ) +
    ylab('Copeland score')
  
  save_pic('cop_score_boxplot')
  
  
  
  
  
  #### !
  
  a__1 <- which(data$IVSM >= quantile(data$IVSM, 0.75))
  a__2 <- which(data$Building_PC >= quantile(data$Building_PC, 0.75))
  a__3 <- which(data$VAR_PERC <= quantile(data$VAR_PERC, 0.25))
  
  Reduce(intersect, list(a__1, a__2, a__3))
  View(data %>% select(PROCOM, IVSM, Building_PC, VAR_PERC) %>% .[
    Reduce(intersect, list(a__1, a__2, a__3)),])
  
  
  data %>% select(AGMAX_50, PROCOM) %>% .[  Reduce(intersect, list(a__1, a__2, a__3)),] %>% mutate(ag_cl = cut(AGMAX_50, c(0, 0.05, 0.15, 0.25, 0.35))) %>%
    group_by(ag_cl) %>% summarise(n())
  
  italy_@data$ag_cl
  data$sel_cr <- NA
  data$sel_cr[ Reduce(intersect, list(a__1, a__2, a__3))] <- T
  sum(is.na(data$sel_cr))
  
  data$ag_cl <- data$AGMAX_50 %>% cut(., c(0, 0.05, 0.15, 0.25, 0.35))
  
  italy_@data$ag_cl <- data$AGMAX_50 %>% cut(., c(0, 0.05, 0.15, 0.25, 0.35))
  data %>% filter(sel_cr == T) %>% group_by(ag_cl) %>%summarise(n())
  ggplot(data %>% filter(sel_cr == T), aes(x = ag_cl)) + geom_bar(aes(y = (..count..)/sum(..count..)), width = 0.7) +
    xlab('a(g)') +
    ylab('%')
  
  save_pic("proportions by AG")
  
  
  italy_@data %<>% left_join(data %>% select(PROCOM, ag_cl, sel_cr))
  
  italy_fortified_ <- italy_fortified %>% left_join(italy_@data)%>%
    mutate(ag_cl = as.character(ag_cl))%>% mutate(ag_cl = ifelse(is.na(sel_cr), NA, ag_cl))
  p <- ggplot(italy_fortified_%>% 
                dplyr::select(long, lat, group,ag_cl)
  ) +  
    theme_minimal()+
    geom_polygon( aes(x = long, y = lat, 
                      group = group, 
                      fill =  ag_cl), alpha = 0.9, size = 0.0015, col = 'black'
    )+ 
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18))+
    xlab('longitude') +
    ylab('latitude') +
    scale_x_continuous(breaks = seq(7, 19, by = 2)) +
    scale_y_continuous(breaks = seq(35, 49, by = 2))+
    coord_map(projection = "lambert", parameters = c(lat0 = 35 , lat1 = 49)) +
    scale_fill_manual('ag[max] class',values = tim.colors(20)[c(1, 7, 13, 20)], breaks = levels(data$ag_cl),na.value="gray", na.translate = T)+
    geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
    geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
    geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
    scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                        , "ag[max] > 0.15" = "dashed",
                                        'ag[max] > 0.25'= 'solid')
    ) +
    guides(linetype = guide_legend(order = 1))
  p
  
  
  
  save_MoI(name = paste0("MoI_aggr_cont"))
  
  