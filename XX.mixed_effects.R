library(lme4)
library(magrittr) 
library(dplyr)
source("00_Functions.R")
library(patchwork)

datasets <- c("data_dbscan",'data_test')

datasets %>%
  sapply(function(i){
    assign(x = i,
           value = read.table(paste0("data/",i,".tsv"), header = T , sep = "\t",
                              na.strings = c(""," ","NA","88","99",
                                             "99.99","777.77","999",
                                             "77.77","77.7","77")),
           envir = .GlobalEnv)
  }) %>% invisible()


lmer( ga ~ I(crl) + I(crl^2) + (1|enrid), data = data_dbscan)

model_lmer <- function(crl){
  6.2679 + crl * 1.53902 - 0.07713 * crl^2
} 

pl_data <- data_test %>% transmute(ga, crl) %>% 
  mutate(garbhini = predict_ga(crl,method = 'Garbhini1') - ga,
         mixed_effects = model_lmer(crl) - ga)
p1 <- pl_data %>% 
  ggplot()+
  geom_density(aes(x = garbhini, color = 'Garbhini')) + 
  geom_density(aes(x = mixed_effects, color = 'lmer')) + 
  theme_bw()
  
p2 <- pl_data %>% 
  ggplot(aes(x = crl))+
  geom_point(aes(y = garbhini, color = 'Garbhini')) + 
  geom_point(aes(y = mixed_effects, color = 'lmer')) + 
  theme_bw()
