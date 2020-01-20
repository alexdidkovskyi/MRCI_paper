# This script is the first step of the analysis of MRCI dataset
# It is used to load preprocessed data and necessary functions

lapply(c('feather','dplyr', 'geosphere', 'readr', 'openxlsx', 'sf'), 
       require, character.only = T)


## Useful functions:

get_filename <- function(){
  path <- rstudioapi::getSourceEditorContext()$path
  gsub('.*/|\\.R', '',path)
}


MoI_p <- function(data_, fill_,  m_col = NA,
                  sc_fill_title = "",
                  sc_fill_value, 
                  zone_b = F,
                  col_t = 'm',
                  gradientn_na_val = 'grey50',
                  sc_fill_breaks = waiver(),
                  legend_k_unit = NULL
                  ){
  p <- ggplot(data_) +
    theme_minimal()+
    geom_polygon( aes(x = long, y = lat, 
                      group = group, 
                      fill = !!ensym(fill_)), 
                  alpha = 0.9, size = 0.0015,  
                  color = m_col
                  
    )+
    theme(axis.title   = element_text(size = 18),
          axis.text    = element_text(size = 18),
          legend.text  = element_text(size = 18),
          legend.title = element_text(size = 18),
          legend.key.height = legend_k_unit)+
    xlab('longitude') +
    ylab('latitude') +
    scale_x_continuous(breaks = seq(7, 19, by = 2)) +
    scale_y_continuous(breaks = seq(35, 49, by = 2))+
    coord_map(projection = "lambert", parameters = c(lat0 = 35 , lat1 = 49))
     
    if (col_t == 'm') {
      p <- p + scale_fill_manual(sc_fill_title, values = sc_fill_value, na.translate = F)
    } else {
      if (col_t == 'g_n') {
        p <-p+ scale_fill_gradientn(sc_fill_title, colors = sc_fill_value, na.value = gradientn_na_val,
                                    breaks = sc_fill_breaks
        )
      }
    }
  
    if (zone_b == T) {
      p <- p +
        geom_path(data = italy_ag_high_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.25'), color = "black", size = 0.35)+
        geom_path(data = italy_ag_moderate_borders, aes(x = long, y = lat, group = group, linetype = 'ag[max] > 0.15'),color = "black", size = 0.35)+
        geom_path(data = italy_ag_median_borders, aes(x = long, y = lat, group = group, lty = 'ag[max] > 0.05'),color = "black", size = 0.35) +
        scale_linetype_manual('',values = c("ag[max] > 0.05" = "dotted"
                                            , "ag[max] > 0.15" = "dashed",
                                            'ag[max] > 0.25'= 'solid'))+
        guides(linetype = guide_legend(order = 1))
    }
  
  
  return(p)
}

Scatter_p <- function(data_, y_, ylab_) {
  p <- ggplot(data_, aes(x = AGMAX_50, y = !!(ensym(y_)), color = log10(POP_2011), weight = log10(POP_2011))) + geom_point() +
    xlim(0.15, 0.28)  +
    scale_color_gradientn(TeX('$\\log_{10}POP'), colors = rainbow(21)[5:21])+
    xlab('AG') +
    ylab(ylab_)+
    theme(
      axis.title   = element_text(size = 18),
      axis.text    = element_text(size = 14),
      legend.text  = element_text(size = 14),
      legend.title = element_text(size = 18)
    )+
    geom_vline(xintercept = 0.25, lty = 2)
}

save_pic <- function(name, .p_path = pic_path){
  ggsave(filename = paste0(.p_path, name,".png"), width = 8, height = 6, dpi = 200)
}
save_MoI <- function(name, .p_path = pic_path){
  ggsave(filename = paste0(.p_path, name,".png"), width = 9, height = 6, dpi = 500)
}
save_MoI_small <- function(name, .p_path = pic_path){
  ggsave(filename = paste0(.p_path, name,".png"), width = 20, height = 12, dpi = 100)
}




create_folder_pics <- function(pth = path){
    dir.create(paste0('results/pics/', pth))
  }
  
## Load preprocessed data
temp_data <- readRDS('data/created/all_data_2018.RDS')

### MRCI dataset
data <- temp_data[[1]]

### sp.data.frame Italy (municipalitie`s boudaries)
italy <- temp_data[[2]]

### data.frame Italy (municipalitie`s boudaries)
italy_fortified <- temp_data[[3]]

italy_ag_moderate_borders <- temp_data[[4]]

italy_ag_high_borders <- temp_data[[5]]

italy_ag_median_borders <- temp_data[[6]]

rm(temp_data)
