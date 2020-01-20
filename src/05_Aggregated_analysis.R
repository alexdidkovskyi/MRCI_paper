

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
  
  
  cop_list <- list()
  j<- list(c(1:3), c(1:4), c(1:4))
  w__ <- list(NULL, NULL, c(1/6,1/6,1/6,1/2))
  for (i in 1:3){
    print(i)
    cop_list[[i]] <- get_copeland(df_ = rank_df[j[[i]],], w_ = w__[[i]])
  }
  
  
  cop_df <- do.call('rbind', cop_list)%>% 
    t() %>% 
    data.frame() %>% 
    set_colnames(c('x','y','z')) %>% 
    mutate(PROCOM =  sort(rev(rank_df[1,]))) %>%
    left_join(data_temp %>% select(PROCOM, AGMAX_50))
  
  
  
  
  ggplot(cop_df %>% 
           mutate(AG_cl = cut(AGMAX_50, c(0, 0.05, 0.15, 0.25, 0.35))), 
         aes(x = AG_cl, y = x)) +
    geom_boxplot() +
    xlab('AG') +
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 18)
    ) +
    ylab('Copeland score')
  
  save_pic('cop_score_boxplot')
  
  
  ###
  p <- MoI_p(data_ = italy_fortified_%>%
               left_join(data %>% 
                           select(PROCOM, DZCOM, AGMAX_50) %>%
                           left_join(
                             cop_df %>% 
                               select(x, PROCOM, AGMAX_50) %>%
                               rename(Cop_sc = x)
                             )
                         ) %>%
               dplyr::select(long, lat, group, Cop_sc),
             fill_ = "Cop_sc",
             m_col = 'black',
             sc_fill_title = 'Cop. Score',
             col_t = 'g_n', 
             sc_fill_value = tim.colors(21),
             zone_b = T)
  p
  
  save_MoI(name = paste0("MoI_aggr_cont"))
  
  ###
  
  
  
  ### 
  
  italy_ <- italy %>% subset(., DZCOM != 'Mappano')
  italy_ <- italy_[order(italy_$PROCOM),]
  
  IVSM_top_25 <- which(data$IVSM >= quantile(data$IVSM, 0.75))
  Building_top_25 <- which(data$Building_PC >= quantile(data$Building_PC, 0.75))
  VAR_PERC_bottom_25 <- which(data$VAR_PERC <= quantile(data$VAR_PERC, 0.25))
  
  Reduce(intersect, list(IVSM_top_25, Building_top_25, VAR_PERC_bottom_25))
  # View(data %>% select(PROCOM, IVSM, Building_PC, VAR_PERC) %>% .[
  #  Reduce(intersect, list(IVSM_top_25, Building_top_25, VAR_PERC_bottom_25)),])
  
  
  ### # of such municipalities by the classes of ag
  data %>% 
    select(AGMAX_50, PROCOM) %>% 
    .[  Reduce(intersect, 
               list(IVSM_top_25, 
                    Building_top_25, 
                    VAR_PERC_bottom_25)),
        ] %>% 
    mutate(ag_cl = cut(AGMAX_50, c(0, 0.05, 0.15, 0.25, 0.35))) %>%
    group_by(ag_cl) %>% summarise(n())
  
  
  
  data$sel_cr <- NA
  data$sel_cr[ Reduce(intersect, list(IVSM_top_25, Building_top_25, VAR_PERC_bottom_25))] <- T
  
  data$ag_cl <- data$AGMAX_50 %>% cut(., c(0, 0.05, 0.15, 0.25, 0.35))
  
  italy_@data$ag_cl <- data$AGMAX_50 %>% cut(., c(0, 0.05, 0.15, 0.25, 0.35))
  
  
  italy_@data %<>% left_join(data %>% select(PROCOM, ag_cl, sel_cr))
  
  italy_fortified_ <- italy_fortified %>% left_join(italy_@data)%>%
    mutate(ag_cl = as.character(ag_cl))%>% mutate(ag_cl = ifelse(is.na(sel_cr), NA, ag_cl))

  
  p <- MoI_p(data_ =italy_fortified_%>% 
               dplyr::select(long, lat, group,ag_cl),
              fill_ = "ag_cl",
             m_col = 'black',
             sc_fill_title = 'ag[max] class',
             col_t = 'm', 
             sc_fill_value = tim.colors(20)[c(1, 7, 13, 20)],
             zone_b = T)
  p
  
  save_MoI(name = paste0("MoI_aggr_cont"))
  
  