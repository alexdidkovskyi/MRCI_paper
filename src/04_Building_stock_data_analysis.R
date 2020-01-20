# This script is the third step of the analysis of MRCI dataset
# Analysis of building stock data

source('src/01_preprocessed_data_loading.R')

library(ggplot2)
library(ggjoy)
library(viridis)
library(tidyverse)
library(compositions)
library(magrittr)
library(cluster)
library(spdep)
library(transport)
library(plotly)
library(fields)
library(ggtext)
library(binsmooth)
library(latex2exp)
library(ggtern)



  pic_path <- paste0("results/pics/", get_filename(), "/")

  create_folder_pics(get_filename())

  
  ### Building stock. Milano
  data %>% 
    filter(DZCOM == 'Milano' ) %>% 
    select(paste0("ERE", 8:16)) %>% gather(x, y) %>% 
    mutate(x = factor(x, 
                      levels = paste0('ERE', 8:16), 
                      labels = c('<1919', '1919-1945','1946-1960', '1961-1970', '1971-1980',
                                 '1981-1990','1991-2000','2001-2005','>2005'))) %>% 
    ggplot(aes(x, y)) + 
    geom_bar(stat ='identity', width = 0.7)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12))+
    xlab('Year') +
    ylab('%') 
  
  save_pic('building stock distribution Milano')
  
  
  
  
  ### Building stock PCA scores plots
  
  cl <- tim.colors(21)
  cl[11] <- 'black'
  
  scores_list <- list()
  l <- 'E'
    
    clr_data  <-  data %>%
      dplyr::select(paste0(l,8:16))  %>%
      clr() %>% princomp()
    print(clr_data$loadings)
      m <- colMeans(data %>%
                    dplyr::select(paste0(l,8:16))  %>%
                    clr())
    
    for (j in 1:3){
      sq <- seq(-5, 5, 0.5)
      i <- j
      res <- lapply(sq, function(x) clrInv(m + x*clr_data$loadings[,i]*clr_data$sdev[i]))
      
      res_df <- do.call(rbind, res) %>% 
        data.frame() %>%
        mutate(alpha = sq)
      
      res_df <- res_df %>% gather(x, y, -alpha)
      res_df <- res_df %>% mutate(x = factor(x, levels = paste0(l, 8:16), 
                                             labels = c('<1919','1919-1945','1946-1960', '1961-1970', '1971-1980',
                                                        '1981-1990','1991-2000','2001-2005','>2005')))
      p <- ggplot(res_df, 
                  aes(x = x, 
                      y = y, 
                      color = factor(alpha), 
                      group = factor(alpha)
                      )
                  )  +
        geom_path(size = 1.1) + geom_point(size = 3) + 
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 12),
              legend.title = element_blank())+
        scale_colour_manual(values = cl)+ 
        xlab('year') +
        ylab('proportion') 
      print(p)
      save_pic(paste0('PC', i,' scores plot'))
      
      
    }
  
    scores_list[[l]] <- cbind(clr_data$scores[,1]/clr_data$sdev[1],
                              clr_data$scores[,2]/clr_data$sdev[2],
                              clr_data$scores[,3]/clr_data$sdev[3])   
  
  
    
  ### DF of PC scores  
  scores_list <- do.call('cbind', scores_list)
  
  scores_list <- data.frame(scores_list)%>% set_colnames(c('PC1', 'PC2', 'PC3'))
  
  scores_list$PROCOM <- data$PROCOM
  scores_list$AGMAX_50 <- data$AGMAX_50
  scores_list$POP_2011 <- data$POP_2011
  
  
  
  
  ### Maps of Italy colored according to scores divided by sd.
  
  cl <- tim.colors(21)
  cl[11] <- 'gray'
    for (i in 1:3){
    
    if (i == 1) {
      cl_ <- cl[2:16]
      br_ <- seq(-4.5, 3, 0.5)
    } else {
      if (i == 2){
        cl_ <- cl[5:20]
        br_ <- seq(-3.5, 4.5, 0.5)
      } else {
        cl_ <- cl[2:21]
        br_ <- seq(-4.5, 5, 0.5)
        
      }
    } 
    fill_ <- paste0('PC', i)  
    
  
    p <- MoI_p(italy_fortified%>% 
                  dplyr::select(long, lat, group,PROCOM ) %>% 
                  left_join(scores_list, by = c('PROCOM'= 'PROCOM')),
                fill_= !!fill_,
                sc_fill_breaks = br_,
                sc_fill_title = fill_,
                sc_fill_value = cl_,
                legend_k_unit = unit(2, 'cm'),
                m_col = 'black',
                col_t ='g_n')
    p
    
    save_MoI(paste0('MoI colored according to PC', i, ' scores'))
  }
  
  
  
  
  data %<>% mutate(ag_class = AGMAX_50 %>% cut(x = ., breaks = c(0, 0.05, 0.15, 0.25, 0.3), right = F))
  
  ### PERMANOVA. Building stock distribution ~ The class of AG
  t__ <- data %>% 
    select(paste0('E',8:16)) %>% 
    compositions::clr() %>% 
    data.frame() %>% 
    mutate(ag_class = data$ag_class)
  t__$id <- 1:length(t__$E8)
  
  vegan_res <- vegan::adonis(t__ %>% select(-ag_class, -id) ~ ag_class,
                             data = data %>% select(ag_class),
                             method = "euclidean", permutations = 1000)
  
  
  
  vegan_res
  
  # Permutation: free
  # Number of permutations: 500
  
  
  #             Df   SumsOfSqs MeanSqs  F.Model      R2   Pr(>F)   
  # ag_class     4      1161    290.248  47.498 0.02326 0.004975 **
  # Residuals 7949     48682    6.12            0.9799            
  # Total     7952     49681                    1.0000          
  # ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  
  
  ### Scatterplot of AG, PC1 scores. Zones of moderate seismic hazard
  ggplot(scores_list, aes(x = AGMAX_50, y = PC1, color = log10(POP_2011) )) + geom_point() +
    xlim(0.15, 0.28) + 
    xlab('ag[max]') +
    ylab('PC1')+
    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 18)
    )+ geom_vline(xintercept = 0.25, lty = 2) +
    scale_color_gradientn(TeX('$\\log_{10}POP'), colors = rainbow(21)[5:21])
  
  save_pic("Building stock PC1" )
  
  
  
  
  ### The same analysis of the reduced dataset. 
  ### The dataset consists of # of buildings built within 3 timeintervals: before 1919, 1919 - 1980, after 1980.
  
  data_3 <- data %>% mutate(`f_g_3` = E8,
                            `s_g_3` = E9 + E10 + E11 + E12,
                            `th_g_3` = E13 + E14 + E15 + E16) %>% select('f_g_3', 's_g_3', 'th_g_3')
  
  
  m <- colMeans(data_3  %>%
                  clr())
  
  cl <- tim.colors(21)
  cl[11] <- 'black'
  
  data_3_princomp <- data_3 %>% clr() %>% princomp()
 
    
  vec_ <- list()
  for (j in 1:2){
    sc_mir <- 1
    sq <- seq(-5, 5, 0.5)
    i <- j
    
    res <- lapply(sq, function(x) clrInv(m + x* sc_mir*data_3_princomp$loadings[,i]*data_3_princomp$sdev[i]))
    
    vec_[[j]]<- do.call(rbind, res) %>% data.frame() %>% mutate(sq = sq)
    
    res_df <- do.call(rbind, res) %>% data.frame()
    res_df$alpha <- sq
    
    res_df <- res_df %>% gather(x, y, -alpha)
    res_df <- res_df %>% mutate(x = factor(x, levels = c('f_g_3', 's_g_3', 'th_g_3'),
                                           labels = c('before 1919','between 1919 & 1980','after 1980')))
    p <- ggplot(res_df, aes(x = x, y = y, color = factor(alpha), group = factor(alpha)))  +
      geom_path(size = 1.1) + geom_point(size = 3) + 
      theme_bw()+
      theme(axis.text.x = element_text(angle = 0, size = 12),
            axis.text.y = element_text(size = 12),
            axis.title = element_text(size = 12),
            legend.title = element_blank())+
      scale_colour_manual(values = cl)+ 
      xlab('year') +
      ylab('proportion') 
    print(p)
    save_pic(paste0('PC', i,' of scores plot_3_classes'))
    
  }
  
  
  scores_list<- cbind(data_3_princomp$scores[,1]/data_3_princomp$sdev[1], 
                      data_3_princomp$scores[,2]/data_3_princomp$sdev[2],
                      data_3_princomp$scores[,3]/data_3_princomp$sdev[3])   
  
  
  
  
    ggtern()  +
      geom_point(data = data_3 %>% 
                   select(f_g_3, s_g_3, th_g_3), 
                 aes(f_g_3, s_g_3, th_g_3), 
                 alpha = 0.1)+  
      geom_path(data = vec_[[1]], 
                aes(x = f_g_3,y = s_g_3, z= th_g_3), 
                col = 'grey15', size = 1.5)+
      geom_path(data = vec_[[2]], 
                aes(x = f_g_3,y = s_g_3, z= th_g_3), 
                col = 'grey35', size = 1.5)+
      geom_point(data = vec_[[1]]%>% rename(sd = sq), 
                 aes(x = f_g_3,y = s_g_3, z= th_g_3, col = sd), size = 3)+
      geom_point(data = vec_[[2]] %>% rename(sd = sq), 
                 aes(x = f_g_3,y = s_g_3, z= th_g_3, col = sd), size = 3)+
      theme_bw() +
    scale_color_gradientn(colors = tim.colors(21))+
    labs(x = '<1919',
         y = '1919-1980',
         z = '>1980')
  
  save_pic('Ternary plot 3 classes')
  
  
  ### Map Of Italy
  
  
  
  
  cl <- tim.colors(21)
  cl[11] <- 'gray'
  
  i <- 1
  
  
  cl_ <- cl[2:17]
  br_ <- seq(-4.5, 3.5, 0.5)
  
  
  fill_ <- paste0('PC', i)  
  
  p <- MoI_p(italy_fortified%>% 
               dplyr::select(long, lat, group,PROCOM ) %>% 
               left_join(scores_list %>% data.frame()%>% set_colnames(paste0('PC', 1:3)) %>%
                           mutate(PROCOM = data$PROCOM), by = c('PROCOM'= 'PROCOM')),
             fill_= !!fill_,
             sc_fill_breaks = br_,
             sc_fill_title = fill_,
             sc_fill_value = cl_,
             legend_k_unit = unit(2, 'cm'),
             m_col = 'black',
             col_t ='g_n')
  p
  save_MoI(paste0('MoI colored according to PC', i, ' scores_3_classes'))
  
  rm(list = ls())
  
  