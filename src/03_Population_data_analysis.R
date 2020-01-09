source("src/01_preprocessed_data_loading.R")

library(ggplot2)
library(ggjoy)
library(viridis)
library(tidyverse)
library(compositions)
library(fields)
library(ineq)
library(magrittr)
library(cluster)
library(corrr)
library(spdep)
library(transport)
library(fda)


pic_path <- paste0("results/pics/", get_filename(), "/")
create_folder_pics(get_filename())
