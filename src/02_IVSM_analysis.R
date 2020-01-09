source('src/01_preprocessed_data_loading.R')

library(cluster)
library(dendextend)
library(factoextra)
library(transport)
library(viridis)
library(shapes)
library(fdadensity)
library(fdapace)
library(zoo)
library(ggplot2)
library(viridis)
library(magrittr)
library(reshape2)
library(tidyr)
library(dplyr)
library(fields)
library(fda)
library(tidyr)
library(ggjoy)
library(latex2exp)


pic_path <- paste0("results/pics/", get_filename(), "/")
create_folder_pics(get_filename())




### ###############################################  ###
##### Provinces' joyplots based on IVSM distribution #####
### ###############################################  ###
data$ID_REGIONE <- as.numeric(data$ID_REGIONE)
italy_fortified %<>% left_join(data %>% dplyr::select(ID_REGIONE, PROCOM))
ggplot(data %>%
         mutate(ID_REGIONE = if_else(ID_REGIONE > 10, ID_REGIONE + 1, ID_REGIONE),
                ID_REGIONE = if_else(DZREG == "Sardegna", 11, ID_REGIONE))%>%
         arrange(ID_REGIONE)%>%
         mutate(DZREG= gsub('\\/.*','', DZREG))%>% 
         mutate(Region = factor(ID_REGIONE, labels = unique(DZREG) )) , aes(y = as.factor(CODPRO), 
                                                                            x = IVSM, 
                                                                            group =as.factor(CODPRO), # DZREG, 
                                                                            fill = Region )) +
  geom_density_ridges(scale = 10, 
                      color = 'white',
                      alpha = 1) +
  xlab('IVSM') +
  guides(fill = guide_legend(title = 'Region')) +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        #axis.text.y = element_text(size = 18),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_fill_viridis_d(option = 'inferno', begin = 0.25, end = 0.85)

ggsave(filename = paste0(pic_path,'IVSM_pr_joyplot.png'), height = 8, width = 14, dpi = 150)

temp_ <- italy_fortified%>% 
  mutate(ID_REGIONE = if_else(ID_REGIONE > 10, ID_REGIONE + 1, ID_REGIONE),
         ID_REGIONE = if_else(DZREG == "Sardegna", 11, ID_REGIONE)) %>%
  dplyr::select(long, lat, group, DZREG, ID_REGIONE) %>% drop_na%>%
  arrange(ID_REGIONE)


data_joy <- data %>%
  dplyr::select(DZPRO, IVSM, CODPRO,DZREG) %>% 
  filter(DZPRO %in% 
           c('Torino','Milano', 'Brescia','Roma', 'Napoli','Palermo') 
         #c('Milano', "Torino", 'Roma', 'Bologna', 'Venezia','Napoli')
  ) %>%
  left_join(temp_%>%
              dplyr::select(DZREG, ID_REGIONE)  %>% 
              distinct() )%>%
  arrange(ID_REGIONE)


ggplot(data_joy, 
       aes(y = factor(DZPRO, levels = 
                        c('Torino','Milano', 'Brescia','Roma', 'Napoli','Palermo')
                      #c('Torino','Milano','Venezia',  'Bologna' ,'Roma','Napoli') 
       ), 
       x = IVSM, 
       group =as.factor(DZPRO),
       fill= as.factor(ID_REGIONE))) +
  geom_density_ridges(scale = 2, 
                      color = 'white',
                      alpha = 1) +
  ylab('Province') +
  xlab('IVSM') +
  guides(fill = guide_legend(title = 'Region'))+
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = 'none') +
  scale_fill_manual(values = inferno(20, begin= 0.25, end  = 0.85)[unique(data_joy$ID_REGIONE)] )
save_pic('IVSM_top_6_provinces_joyplot_upd')



data %<>% mutate(ag_class = AGMAX_50 %>% cut(x = ., breaks = c(0, 0.05, 0.15, 0.25, 0.3), right = F))



data_for_clust <- data %>% dplyr::select( IVSM, CODPRO, DZPRO)

clust_list <- list()
names <- c()

for (i in unique(data_for_clust$CODPRO)){
  temp_ <- data_for_clust %>% filter(CODPRO == i) %>% dplyr::select( IVSM) %>% pull()
  
  clust_list <- c(clust_list, list(temp_))
  names <- c(names, unique(data_for_clust$DZPRO[data_for_clust$CODPRO == i]))
}
names(clust_list) <- names

res2p <- sapply(clust_list, function(x) {
  lapply(clust_list, 
         function(y) wasserstein1d(x, y, p = 2)
         )
  })

res2p_matr <- matrix(unlist(res2p), nrow = length(names), dimnames = list(names, names))
res2p_hclust_ward <- hclust(as.dist(res2p_matr), method = 'ward.D')

png(paste0(pic_path,'IVSM_wass_pr_dendr.png'), height = 600, width = 1200)
plot(res2p_hclust_ward, xlab=c('distance'), cex = 1, hang = -1, main = '')
dev.off()

cluster <- cutree(res2p_hclust_ward, k = 4)


data_temp <- data %>% 
  select(IVSM, CODPRO, DZPRO, DZREG, ID_REGIONE) %>% 
  left_join(
    data.frame(DZPRO = names(cluster),
                       clt_id = cluster))

grid_IVSM_list <- lapply(sort(unique(data_temp$clt_id)), function(i) {
  id_pro <- data_temp$CODPRO[data_temp$clt_id == i] %>% unique()
  IVSM_t <- data %>% 
    select(IVSM, CODPRO) %>% 
    filter(CODPRO %in% id_pro) %>% 
    select(IVSM) %>% pull() 
  return(sapply(seq(0, 1, by = 0.001), 
                function(x) quantile(IVSM_t, x)))
  
})

data_temp_sb <- 
  data %>% 
  select(IVSM, CODPRO, DZPRO, DZREG, ID_REGIONE) %>% 
  left_join(data.frame(DZPRO = names(cluster),
                       clt_id = cluster))%>% 
  select(IVSM, CODPRO, clt_id) %>%
  group_by(CODPRO, clt_id)%>% 
  summarise(dens = list(density(IVSM, n = 1024)))

id <- 1
res <-
  lapply(1:nrow(data_temp_sb), function(id) {
    dens <- data_temp_sb$dens[[id]]
    clt_id <- data_temp_sb$clt_id[id]
    temp <- approx(dens$x, dens$y, grid_IVSM_list[[clt_id]])
    temp$y[is.na(temp$y)] <- 0
    return(temp$y)
  }) %>% do.call(rbind, .)


res1 <- sapply(1:nrow(data_temp_sb), function(i) {
  1e-14+res[i,]/fdapace::trapzRcpp(grid_IVSM_list[[data_temp_sb$clt_id[i]]], res[i,])
}) %>% t()


distr_list <- c()
i <- 1
for (i in unique(data_temp_sb$clt_id)){
  idx <- which(data_temp_sb$clt_id ==  i)
  distr_list <- c(distr_list, list(getWFmean(res1[idx, ], grid_IVSM_list[[i]], useAlpha = T, alpha = 0.001)))
}
distr_df <- lapply(1:4, function(i) 
  data.frame(density = distr_list[[i]], 
             x = grid_IVSM_list[[i]], 
             Cluster = as.character(i))) %>% 
  do.call('rbind',.)


tbl_ <- list()
for (i in 1:length(distr_list)){
  l_ <- length(distr_list[[i]])
  d_ <- (distr_list[[i]][1:(l_ - 1)] + distr_list[[i]][2:l_])/2 * (grid_IVSM_list[[i]][2:l_] - grid_IVSM_list[[i]][1:(l_ - 1)])
  cs_<- cumsum(d_)
  mn <- sum(d_*(grid_IVSM_list[[i]][2:l_] + grid_IVSM_list[[i]][1:(l_ - 1)])/2)
  mn_2 <- sum(d_*((grid_IVSM_list[[i]][2:l_] + grid_IVSM_list[[i]][1:(l_ - 1)])/2)^2)
  sd_ <- sqrt(mn_2 - mn^2)*sqrt(1001/1000)
  
  id_ <- sapply(seq(0.25, 0.75, 0.25), function(x) min(which(cs_ >= x))) +1
  tbl_[[i]] <- c(mn, sd_, grid_IVSM_list[[i]][id_])
  
}

tbl_ %>% do.call('rbind',.) %>% openxlsx::write.xlsx(., file = 'IVSM_tbl_by_cl')


### Visualizing obtained clusters 

ggplot(distr_df, aes(x = x, y = density, fill = Cluster)) +geom_area(position = "identity", alpha = 0.7) +
  xlab('IVSM')+
  ylab('')+
  theme(title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  scale_fill_manual(values = inferno(4,begin = 0.25, end = 0.85, direction = 1))+ 
  coord_cartesian(ylim = c(0.01, 0.3))



ggplot(distr_df, aes(x = x, y = density, color = Cluster)) +geom_line(size = 2) +
  xlab('IVSM')+
  ylab('')+
  theme(title = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x =  element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.text.y = element_text(size = 16)) +
  scale_color_manual(values = inferno(4,begin = 0.25, end = 0.85, direction = 1))
save_pic('IVSM_wass_pr_4_clt_barycenters')

cluster_1 <- cluster
if ('Cluster' %in% names(italy_fortified)) {
  italy_fortified%<>% select(-Cluster)
}
italy_fortified %<>% left_join(data.frame(DZPRO = names(cluster_1), Cluster = unname(cluster_1),stringsAsFactors = F ))

p <- ggplot(italy_fortified %>%
              dplyr::select(long, lat, group, Cluster)) +
  theme_minimal()+
  geom_polygon( aes(x = long, y = lat, 
                    group = group, 
                    fill = as.character(Cluster)), alpha = 0.9, size = 0.0015, color = 'black'
               
  )+
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18))+
  xlab('longitude') +
  ylab('latitude') +
  scale_x_continuous(breaks = seq(7, 19, by = 2)) +
  scale_y_continuous(breaks = seq(35, 49, by = 2))+
  coord_map(projection = "lambert", parameters = c(lat0 = 35 , lat1 = 49))+
  scale_fill_manual('Cluster', values = inferno(4, begin = 0.25, end = 0.85), na.translate = F)

p

save_MoI('IVSM_wass_pr_4_clt_MoI')

