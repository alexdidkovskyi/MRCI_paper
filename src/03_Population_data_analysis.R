# This script is the third step of the analysis of MRCI dataset
# Analyses of growth rate and the third quartile of population age distribution

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



### Scatterplots of Growth rate and Third quartile of population age distribution 
### for the municipalities that belong to the zones of moderate and high seismic hazard
p <- Scatter_p(data_ = data %>% mutate(gr = (POP_2018/POP_2011 - 1)), y_ = 'gr',  ylab_ = 'Growth rate' )
p
save_pic("2012 to 2018 Growth rate" )

p <- Scatter_p(data_ = data, y = 'ETA_Q3', ylab_ = 'ETA_Q3')
p
save_pic("Eta Q3 2018y")


### Map of Italy colored according to Growth rate and Third quartile of population age distribution 

for (i in c('VAR_PERC', "ETA_Q3")) {
  p <- MoI_p(data_ = italy_fortified %>% select(long, lat, group, PROCOM) %>%
               left_join(data %>% select(PROCOM, !!ensym(i), AGMAX_50)),
             fill_ = i,m_col = 'black',
             sc_fill_title = 'VAR_PERC',
             sc_fill_value =  tim.colors(21),
             col_t = 'g_n')
  print(p)
  if( i == 'VAR_PERC') {
    save_MoI('MoI_colored_acc_to_growth_rate')
  } else{
    save_MoI('MoI_colored_acc_to_ETA_Q3')
  }
}



### Old municipalities` population data
### select municipalities that haven't changed since 1991

trim <- function(x) gsub("^\\s+|\\s+$", "", x)
observ_92_01 <- read.xlsx("data/initial/1991-2001.xlsx", sheet = 1) %>%
  mutate(Comune = trim(Comune)) %>%
  mutate(Comune = gsub("/.*", "", Comune)) %>%
  select(-`1991`)
observ_02_11 <- read.xlsx("data/initial/2001-2011.xlsx", sheet = 1) %>%
  mutate(Comune = trim(Comune)) %>%
  mutate(Comune = gsub("/.*", "", Comune)) %>%
  select(-`2001`)



### identify municipalities` names that appear several times
bad_c <- (table(observ_92_01$Comune) > 1) %>%
  which() %>%
  names()

c_names <- Reduce(intersect, list(observ_02_11$Comune, observ_92_01$Comune, data$DZCOM))
c_names <- setdiff(c_names, bad_c)


data %<>% mutate(ag_class = AGMAX_50 %>% cut(x = ., breaks = c(0, 0.05, 0.15, 0.25, 0.35), right = F))
data_for_q <-
  data %>%
  select(DZCOM, POP_2011, AGMAX_50, DZPRO) %>%
  select(-DZPRO) %>%
  filter(DZCOM %in% c_names) %>%
  rename(Comune = DZCOM) %>%
  left_join(observ_92_01 %>%
              filter(Comune %in% c_names)) %>%
  left_join(observ_02_11 %>%
              filter(Comune %in% c_names)) %>%
  rename(`2012` = POP_2011) %>%
  mutate(prop = `2012` / `1992`) %>%
  mutate(ag_m = AGMAX_50 >= 0.2) %>%
  distinct()



data_for_q_scaled <- data_for_q
data_for_q_scaled[, c(2, 4:23)] <- data_for_q_scaled[, c(2, 4:23)] / data_for_q_scaled$`1992`


data_for_q_scaled_l <-
  data_for_q_scaled %>%
  select(-prop, -ag_m) %>%
  gather(year, population, -Comune, -AGMAX_50) %>%
  mutate(year = as.numeric(year))

### Log-growth curves
ggplot(
  data_for_q_scaled_l,
  aes(
    x = year, y = log10(population),
    group = Comune, col = AGMAX_50
  )
) +
  geom_line() +
  scale_color_gradientn('ag[max]', colors = tim.colors(21)) +
  ylim(-1 / 2, 1 / 2) +
  ylab("Log of proportion to 1992") +
  xlab("Year") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )

save_pic(paste0("Log growth curves "))

data_for_q_scaled_m <-
  data_for_q_scaled %>%
  left_join(data %>% select(DZCOM, ag_class) %>%
              distinct(),
            by = c("Comune" = "DZCOM")
  ) %>%
  select(-prop, -ag_m, -AGMAX_50)



data_for_q_scaled_m <- data.frame(data_for_q_scaled_m, stringsAsFactors = F)
data_for_q_scaled_m <- do.call(data.frame, lapply(data_for_q_scaled_m, function(x) replace(x, is.infinite(x), NA)))

### Set all NA proportions to 1
data_for_q_scaled_m[is.na(data_for_q_scaled_m)] <- 1

colnames(data_for_q_scaled_m) <- gsub("X", "", colnames(data_for_q_scaled_m))

### Log growth year-wise medians

median_df <- sapply(unique(data$ag_class), function(ag_cl) {
  data_for_q_scaled_m %>%
    filter(ag_class == ag_cl) %>%
    select(-ag_class) %>%
    set_rownames(.$Comune) %>%
    select(-Comune) %>%
    apply(., 2, median)
}) %>%
  data.frame() %>%
  set_colnames(unique(data$ag_class)) %>%
  mutate(Year = as.numeric(row.names(.))) %>%
  gather(ag_class, prop, -Year) %>%
  mutate(ag_class = factor(ag_class, levels = levels(data$ag_class)))

ggplot(median_df, aes(x = Year, y = log10(prop), col = ag_class)) +
  geom_line(size = 2) +
  scale_color_manual(values = tim.colors(5)) +
  ylim(-0.05, 0.05) +
  ylab("Log of proportion to 1992") +
  xlab("Year") +
  theme_bw()  +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18)
  )+
  labs(color = "ag[max] class") 
save_pic(paste0("Log growth year-wise medians by ag class", lbl))

### ########  ###
##### FPCA: growth rate #####
### ########  ###


data_for_q_scaled_l <-
  data_for_q_scaled %>%
  select(-prop, -ag_m) %>%
  gather(year, population, -Comune, -AGMAX_50) %>%
  mutate(year = as.numeric(year))

data_for_q_scaled_pr <- data_for_q_scaled[, c(4:23,2)] %>% log10(.)

### All NA/inf/NaN set to 0
data_for_q_scaled_pr <-do.call(data.frame,lapply(data_for_q_scaled_pr, function(x) replace(x, is.infinite(x),0)))
data_for_q_scaled_pr <-do.call(data.frame,lapply(data_for_q_scaled_pr, function(x) replace(x, is.nan(x),0)))
data_for_q_scaled_pr <-do.call(data.frame,lapply(data_for_q_scaled_pr, function(x) replace(x, is.na(x),0)))
##### fda

### Functional PCA using 11 bsplines

growth_basis_bspline <- create.bspline.basis(c(1, 21), nbasis = 11, norder = 4)
growth_fd <- smooth.basis(1:21, data_for_q_scaled_pr  %>% 
                            as.matrix() %>% 
                            t(), 
                          growth_basis_bspline)
growth_pca <- pca.fd(growth_fd$fd, nharm = 1)

###  plot harmonics
png(paste0(pic_path,'FPCA_log_growth_rate', '_.png'), width = 600, height = 480)
par(ask=F)
plot.pca.fd(growth_pca, cex = 1.5, expand = c(1:3))
dev.off()


nx <- 21
pcafd <- growth_pca 

argvals <- seq(pcafd[[1]]$basis$rangeval[1], pcafd[[1]]$basis$rangeval[2], length = nx)
meanmat <- eval.fd(argvals, pcafd$meanfd)[,1]

vecharm <- eval.fd(argvals, pcafd[[1]])[,1]

fpca_score <- lapply(sqrt(pcafd$values[1])*c(seq(-3, 3, 0.5)), function(x) meanmat + vecharm*x) %>% do.call('cbind',.) %>% data.frame() %>%
  set_colnames(as.factor(seq(-3, 3, 0.5))) %>%
  mutate(year = seq(1992, 2012, 1)) %>% tidyr::gather('sd_','val', -year) %>%
  mutate(sd_ = factor(sd_, levels = as.character(seq(-3, 3, 0.5))))

cl <- tim.colors(13)
cl[7] <- 'black'

### Fpca first PC score +- 3sd

ggplot(fpca_score, aes(x = year, y = val, color = sd_))+geom_path(size = 1.1) + geom_point(size = 3) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_blank())+
  scale_colour_manual(values = cl) +
  xlab('Year') +
  ylab(latex2exp::TeX('$\\log_{10}PGR'))

save_pic('FPCA_log_gr_score_plot')


data %<>% left_join(data_frame(DZCOM = data_for_q$Comune,
                               growth_PC1 = growth_pca$scores[,1]))


### Map of Italy colored accoding to Functional PC1 scores (in terms of sd). 
### For the purpose of visualization all values higher than 5 set to 5.

p <- MoI_p(data_ =  
             italy_fortified %>% 
             left_join(data %>% 
                         select(PROCOM, growth_PC1) %>% 
                         mutate(growth_PC1 = growth_PC1/sqrt(growth_pca$values[1])
                         )
             )  %>% 
             mutate(growth_PC1 = ifelse(growth_PC1 >5, 5, growth_PC1)),
           fill_ = 'growth_PC1',
           col_t = 'g_n',
           m_col = 'black',
           sc_fill_title = 'GR PC1',
           sc_fill_value = tim.colors(21),
           gradientn_na_val = 'white')

p
