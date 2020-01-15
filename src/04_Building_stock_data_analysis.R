
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


  pic_path <- paste0("results/pics/", get_filename(), "/")

  create_folder_pics(get_filename())

  
  compositional_data <- data %>%  select(paste0("E", 8:16)) 
  compositional_data <- (compositional_data / rowSums(compositional_data) ) 
  a <- compositional_data %>% clr() %>% princomp()
  summary(a)
  
  data %>% filter(DZCOM == 'Milano' ) %>% select(paste0("ERE", 8:16)) %>% gather(x, y) %>% 
    mutate(x = factor(x, levels = paste0('ERE', 8:16), 
                      labels = c('<1919','1919-1945','1946-1960', '1961-1970', '1971-1980',
                                 '1981-1990','1991-2000','2001-2005','>2005'))) %>% ggplot(aes(x, y)) + geom_bar(stat ='identity', width = 0.7)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12))+
    xlab('Year') +
    ylab('%') 
  save_pic('building stock distribution Milano')
  
  
  
  
  
  cl <- tim.colors(21)
  cl[11] <- 'black'
  
  scores_list <- list()
  for (l in c('E')){
    
    clr_data  <-  data %>%
      dplyr::select(paste0(l,8:16))  %>%
      clr() %>% princomp()
    print(clr_data$loadings)
    
    #Loadings:
    #  Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7 Comp.8 Comp.9
    #E8  -0.762  0.369 -0.344  0.167        -0.105               -0.333
    #E9  -0.397 -0.385  0.321 -0.540  0.244  0.286 -0.207        -0.333
    #E10        -0.505  0.254  0.209 -0.192 -0.487  0.392 -0.311 -0.333
    #E11  0.102 -0.284         0.471 -0.200        -0.340  0.640 -0.333
    #E12  0.188 -0.112 -0.298  0.289  0.129  0.497 -0.177 -0.611 -0.333
    #E13  0.286        -0.427 -0.231  0.359         0.581  0.320 -0.333
    #E14  0.270  0.132 -0.264 -0.373        -0.550 -0.521 -0.115 -0.333
    #E15  0.204  0.335  0.242 -0.259 -0.685  0.316  0.202        -0.333
    #E16  0.140  0.486  0.560  0.268  0.488                      -0.333
    
    
    #Importance of components:
    #                         Comp.1    Comp.2    Comp.3     Comp.4     Comp.5     Comp.6     Comp.7     Comp.8     Comp.9
    #Standard deviation     1.6729607 1.1608988 0.8908210 0.66703153 0.57895553 0.49998339 0.42288672 0.32270197      0
    #Proportion of Variance 0.4475842 0.2155222 0.1269066 0.07115346 0.05360357 0.03997739 0.02859904 0.01665354      0
    #Cumulative Proportion  0.4475842 0.6631064 0.7900130 0.86116646 0.91477003 0.95474743 0.98334646 1.00000000      1
    
    #l <- 'E'
    m <- colMeans(data %>%
                    dplyr::select(paste0(l,8:16))  %>%
                    clr())
    
    #Reverse EME 2nd principle component 
    for (j in 1:3){
      #score mirroring
      #sc_mir <- 1
      sq <- seq(-5, 5, 0.5)
      i <- j
      #if (j == 1) {sc_mir <- -1}
      res <- lapply(sq, function(x) clrInv(m + x* sc_mir*clr_data$loadings[,i]*clr_data$sdev[i]))
      
      res_df <- do.call(rbind, res) %>% data.frame()
      res_df$alpha <- sq
      
      res_df <- res_df %>% gather(x, y, -alpha)
      res_df <- res_df %>% mutate(x = factor(x, levels = paste0(l, 8:16), 
                                             labels = c('<1919','1919-1945','1946-1960', '1961-1970', '1971-1980',
                                                        '1981-1990','1991-2000','2001-2005','>2005')))
      p <- ggplot(res_df, aes(x = x, y = y, color = factor(alpha), group = factor(alpha)))  +
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
    j <- 1
    #if (j == 1) {sc_mir <- -1}
    
    scores_list[[l]] <- cbind(clr_data$scores[,1]/clr_data$sdev[1] * sc_mir, 
                              clr_data$scores[,2]/clr_data$sdev[2],
                              clr_data$scores[,3]/clr_data$sdev[3])   
  }
  
  scores_list <- do.call('cbind', scores_list)
  
  scores_list <- data.frame(scores_list)%>% set_colnames(c('PC1', 'PC2', 'PC3'))
  
  scores_list$PROCOM <- data$PROCOM
  scores_list$AGMAX_50 <- data$AGMAX_50
  scores_list$POP_2011 <- data$POP_2011
  
  ### Custom selection
  
  
  cl <- tim.colors(21)
  cl[11] <- 'gray'
  i <- 1
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
    j <- paste0('PC', i)  
    p <- ggplot(italy_fortified%>% 
                  dplyr::select(long, lat, group,PROCOM ) %>% 
                  left_join(scores_list, by = c('PROCOM'= 'PROCOM'))) +  
      theme_minimal()+
      geom_polygon( aes_string(x = quote(long), y = quote(lat), 
                               group = quote(group), 
                               fill =  j), alpha = 0.9, size = 0.0015, col = 'black'
      )+ 
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 18),
            legend.title = element_text(size = 18),
            legend.key.height = unit(2, 'cm'))+
      xlab('longitude') +
      ylab('latitude') +
      scale_x_continuous(breaks = seq(7, 19, by = 2)) +
      scale_y_continuous(breaks = seq(35, 49, by = 2))+
      coord_map(projection = "lambert", parameters = c(lat0 = 35 , lat1 = 49))+
      scale_fill_gradientn(colors = cl_,breaks=br_)
    p
    save_MoI(paste0('MoI colored according to PC', i, ' scores'))
  }
  
  
  
  
  data %<>% mutate(ag_class = AGMAX_50 %>% cut(x = ., breaks = c(0, 0.05, 0.15, 0.25, 0.3), right = F))
  
  ### Permutation tests
  t__ <- data %>% select(paste0('E',8:16)) %>% compositions::clr() %>% data.frame() %>% mutate(ag_class = data$ag_class)
  t__$id <- 1:length(t__$E8)
  
  vegan_res <- vegan::adonis(t__ %>% select(-ag_class, -id)~ag_class,
                             data = data %>% select(ag_class), method = "euclidean", permutations = 1000)
  
  
  
  vegan_res
  
  # Permutation: free
  # Number of permutations: 500
  
  
  #             Df   SumsOfSqs MeanSqs  F.Model      R2   Pr(>F)   
  # ag_class     4      1161    290.248  47.498 0.02326 0.004975 **
  # Residuals 7949     48682    6.12            0.9799            
  # Total     7952     49681                    1.0000          
  # ---
  #  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  
  
  
  
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
  
  
  
  
  
  data_3 <- data %>% mutate(`f_g_3` = E8,
                            `s_g_3` = E9 + E10 + E11 + E12,
                            `th_g_3` = E13 + E14 + E15 + E16) %>% select('f_g_3', 's_g_3', 'th_g_3')
  
  m <- colMeans(data_3  %>%
                  clr())
  
  cl <- tim.colors(21)
  cl[11] <- 'black'
  
  data_3_princomp <- data_3 %>% clr() %>% princomp()
  #Reverse EME 2nd principle component 
  j <- 1
  vec_ <- list()
  for (j in 1:3){
    #score mirroring
    sc_mir <- 1
    sq <- seq(-5, 5, 0.5)
    i <- j
    #if (j == 1) {sc_mir <- -1}
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
  j <- 1
  #if (j == 1) {sc_mir <- -1}
  
  scores_list<- cbind(data_3_princomp$scores[,1]/data_3_princomp$sdev[1] * sc_mir, 
                      data_3_princomp$scores[,2]/data_3_princomp$sdev[2],
                      data_3_princomp$scores[,3]/data_3_princomp$sdev[3])   
  
  
  
  
  library(ggtern)
  library(viridis)
  ggtern()  +geom_point(data = data_3 %>% select(f_g_3, s_g_3, th_g_3), aes(f_g_3, s_g_3, th_g_3), alpha = 0.1)+  
    geom_path(data = vec_[[1]], aes(x = f_g_3,y = s_g_3, z= th_g_3), col = 'grey15', size = 1.5)+
    geom_path(data = vec_[[2]], aes(x = f_g_3,y = s_g_3, z= th_g_3), col = 'grey35', size = 1.5)+
    geom_point(data = vec_[[1]]%>% rename(sd = sq), aes(x = f_g_3,y = s_g_3, z= th_g_3, col = sd), size = 3)+
    geom_point(data = vec_[[2]] %>% rename(sd = sq), aes(x = f_g_3,y = s_g_3, z= th_g_3, col = sd), size = 3)+theme_bw() +
    scale_color_gradientn(colors = tim.colors(21))+
    labs(x = '<1919',
         y = '1919-1980',
         z = '>1980')
  
  save_pic('Ternary plot 3 classes')
  
  
  ### Map Of Italy
  
  
  
  
  cl <- tim.colors(21)
  cl[11] <- 'gray'
  
  i <- 1
  
  if (i == 1) {
    cl_ <- cl[2:17]
    br_ <- seq(-4.5, 3.5, 0.5)
  }
  
  j <- paste0('PC', i)  
  p <- ggplot(italy_fortified%>% 
                dplyr::select(long, lat, group,PROCOM ) %>% 
                left_join(scores_list %>%data.frame()%>% set_colnames(paste0('PC', 1:3)) %>%
                            mutate(PROCOM = data$PROCOM), by = c('PROCOM'= 'PROCOM'))) +  
    theme_minimal()+
    geom_polygon( aes_string(x = quote(long), y = quote(lat), 
                             group = quote(group), 
                             fill =  j), alpha = 0.9, size = 0.0015, col = 'black'
    )+ 
    theme(axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          legend.text = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.key.height = unit(2, 'cm'))+
    xlab('longitude') +
    ylab('latitude') +
    scale_x_continuous(breaks = seq(7, 19, by = 2)) +
    scale_y_continuous(breaks = seq(35, 49, by = 2))+
    coord_map(projection = "lambert", parameters = c(lat0 = 35 , lat1 = 49))+
    scale_fill_gradientn(colors = cl_,breaks=br_)
  p
  save_MoI(paste0('MoI colored according to PC', i, ' scores_3_classes'))
  
  rm(list = ls())
  
  