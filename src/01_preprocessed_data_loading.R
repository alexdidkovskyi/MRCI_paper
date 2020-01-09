lapply(c('feather','dplyr', 'geosphere', 'readr', 'openxlsx', 'sf'), 
       require, character.only = T)


### Useful functions:

get_filename <- function(){
  path <- rstudioapi::getSourceEditorContext()$path
  gsub('.*/|\\.R', '',path)
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
  
### Loading preprocessed data
temp_data <- readRDS('data/created/all_data_2018.RDS')

#### MRCI dataset
data <- temp_data[[1]]

#### sp.data.frame Italy (municipalitie`s boudaries)
italy <- temp_data[[2]]

#### data.frame Italy (municipalitie`s boudaries)
italy_fortified <- temp_data[[3]]

italy_ag_moderate_borders <- temp_data[[4]]

italy_ag_high_borders <- temp_data[[5]]

italy_ag_median_borders <- temp_data[[6]]

rm(temp_data)
