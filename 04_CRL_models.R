
datasets <- c("data_main","data_dbscan","data_clinical")
options(digits =4)

methods <- c("Hadlock",
             "Intergrowth")


# Library Imports ---------------------------------------------------------

library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function as . or assigns to first parameter by default 
library(ggplot2) #plotting functions
library(dplyr)
library(patchwork)
library(tidyr)
library(readr)
source("00_Functions.R")

theme_set(theme_bw())
set.seed(666)
seed_list <- sample(1:10^5,10^2)

# loading datasets --------------------------------------------------------

lapply(datasets, function(i){
  tmp <- read.table(paste0("data/",i,".tsv"), sep = "\t", header = TRUE,
                    na.strings = c(""," ","NA","88","99",
                                   "99.99","777.77","999",
                                   "77.77","77.7","77"))
  tmp <- tmp[ !(is.na(tmp$crl)|is.na(tmp$ga)),]
  assign(i,
         value = tmp,
         envir = .GlobalEnv)
}) %>% invisible()

datasets[1] <- "data_training"
data_training <- data_main
# generate fractional polynomial formulae ---------------------------------

formulas_crl <- generate_frac_poly_formulas(outcome = "crl",
                                            predictor = "ga")

formulas_ga <- generate_frac_poly_formulas(outcome = "ga",
                                           predictor = "crl")



# Adding data from hadlock ------------------------------------------------

data_hadlock <- seq(0.1,13,0.1) %>% 
  lapply(function(crl){
    data.frame(ga = rnorm(n = 100 , 
                          mean = predict_ga(crl = crl , method = "Hadlock"),
                          sd = get_sd_published_formulas(crl = crl , method = "Hadlock")),
               crl)
  }) %>% bind_rows() %>% filter(ga > 15, ga < 18) %>% 
  {
    set.seed(666)
    selection <- sample(x = 1:nrow(.),
                        size = 3 * nrow(data_dbscan)/(max(data_dbscan$crl) - min(data_dbscan$crl)),replace = FALSE)
    .[selection,]
  }
data_dbscan_hadlock3 <- rbind(data_dbscan[,1:2], data_hadlock ) %>% 
  mutate(is_from_dbscan = c(rep(TRUE,nrow(data_dbscan)),rep(FALSE, nrow(data_hadlock))))


mccv_results_ga <- lapply(seed_list, function(seed){
  set.seed(seed)
  selection <- sample(x = 1:nrow(data_dbscan_hadlock3),
                      size = as.integer(0.7 * nrow(data_dbscan_hadlock3)))
  train_set <- data_dbscan_hadlock3[selection,]
  test_set <- data_dbscan_hadlock3[-selection,]
  
  lapply(formulas_ga, function(formula){
    
    pred <- lm(formula = as.formula(formula),data = train_set) %>% 
      predict(object = ., test_set)
    
    data.frame(formula,
               seed,
               R2 = R2(test_set$ga, pred),
               RMSE = RMSE(test_set$ga, pred))
  }) %>% bind_rows()
}) %>% bind_rows()

# p1 <- ggplot(mccv_results_ga, aes(x = R2 , y = formula)) +
#   geom_violin()
# p2 <- ggplot(mccv_results_ga, aes(x = RMSE , y = formula)) +
#   geom_violin()
# p1 | p2

#selecting formula with best median R2 in mccv
best_formula_hadlock3 <- mccv_results_ga %>% group_by(formula) %>% 
  summarise(R2 = median(R2),
            RMSE = median(RMSE), .groups = "drop") %>% 
  arrange(-R2) %>% 
  {.$formula[1]}

mccv_results_ga %>% group_by(formula) %>% 
  summarise( R2 = median(R2),
            RMSE = median(RMSE), .groups = "drop") %>% 
  arrange(-R2) %>% write_tsv(file = "tables/crl_models_+hadlock3.tsv")

model_dbscan_hadlock3 <-lm(data = data_dbscan_hadlock3 ,
                           formula = best_formula_hadlock3)

ggplot(data_dbscan_hadlock3, aes( x = crl , y = ga)) +
  geom_point(aes(shape = is_from_dbscan),alpha = 0.2) +
  geom_line( data = data.frame( x = data_dbscan_hadlock3$crl,
                                y = model_dbscan_hadlock3$fitted.values),
             aes(x =x , y = y ), size = 1.5)

# Creating Supplementary fig 1 ------------------------------------------------

plot_list <- lapply(datasets, function(dataset){
  
  data <- get(dataset) %>% select(ga , crl)
  data_hadlock <- seq(0.1,13,0.1) %>% 
    lapply(function(crl){
      data.frame(ga = rnorm(n = 100 , 
                            mean = predict_ga(crl = crl , method = "Hadlock"),
                            sd = get_sd_published_formulas(crl = crl , method = "Hadlock")),
                 crl)
    }) %>% bind_rows() %>% filter(ga > 15, ga < 18) %>% 
    {
      set.seed(666)
      selection <- sample(x = 1:nrow(.),
                          size = 3 * nrow(data)/(max(data$crl) - min(data$crl)),replace = FALSE)
      .[selection,]
    }
  
  data_hadlock3 <- rbind(data, data_hadlock)
  data_hadlock3$data_source <- "Hadlock"
  data_hadlock3$data_source[1:nrow(data)] <- "GARBHINI"
  
  ggplot(data_hadlock3, aes(x = crl, y = ga , color = data_source )) +
    geom_point(alpha =0.1)+
    labs(title = paste0(dataset, " n = ", sum(data_hadlock3$data_source == "GARBHINI"),
                        "; hadlock n = ", sum(!data_hadlock3$data_source == "Hadlock")),
         x = "Crown Rump Length(cm)",
         y = "Gestational Age(weeks)") + 
    ylim(0,24) + xlim(0,10)
  
  })

wrap_plots(plot_list, ncol = 1, guides = "collect") %>% 
  ggsave(filename = "plots/Suppl_Fig1.png", plot = . ,
         height = 7.5 , width = 6 , dpi = 600)

# Creating Supplementary fig 2 ------------------------------------------------

plot_list <- lapply(datasets, function(dataset){
  
  data <- get(dataset) %>% select(ga , crl)
  data_hadlock <- seq(0.1,13,0.1) %>% 
    lapply(function(crl){
      data.frame(ga = rnorm(n = 100 , 
                            mean = predict_ga(crl = crl , method = "Hadlock"),
                            sd = get_sd_published_formulas(crl = crl , method = "Hadlock")),
                 crl)
    }) %>% bind_rows() %>% filter(ga > 15, ga < 18) %>% 
    {
      set.seed(666)
      selection <- sample(x = 1:nrow(.),
                          size = 3 * nrow(data)/(max(data$crl) - min(data$crl)),replace = FALSE)
      .[selection,]
    }
  
  data_hadlock3 <- rbind(data, data_hadlock)
  
  mccv_results <- lapply(seed_list, function(seed){
    set.seed(seed)
    selection <- sample(x = 1:nrow(data_hadlock3),
                        size = as.integer(0.7 * nrow(data_hadlock3)))
    train_set <- data_hadlock3[selection,]
    test_set <- data_hadlock3[-selection,]
    
    lapply(formulas_ga, function(formula){
      
      pred <- lm(formula = as.formula(formula),data = train_set) %>% 
        predict(object = ., test_set)
      
      data.frame(formula,
                 seed,
                 R2 = R2(test_set$ga, pred),
                 RMSE = RMSE(test_set$ga, pred))
    }) %>% bind_rows()
  }) %>% bind_rows()
  
  best_formulas <- mccv_results %>% group_by(formula) %>% 
    summarise(R2 = median(R2),
              RMSE = median(RMSE), .groups = "drop") %>% 
    arrange(-R2) %>% 
    {.$formula[1:3]}
  
  plot_data <- sapply(best_formulas, function(formula){
    lm(formula, data = data_hadlock3)$fitted.values
  }) %>% 
    cbind(data_hadlock3,.) %>% 
    pivot_longer(cols = -c(ga,crl) ,
                 names_to = "formula",
                 values_to = "ga_pred")
  ggplot(plot_data, aes(x = crl )) +
    geom_point(aes(y = ga), alpha = 0.01) +
    geom_line(aes(y = ga_pred, color = formula), size = 1.1)+
    labs(title = paste0(dataset, ", n = ", nrow(get(dataset))),
         x = "Crown Rump Length(cm)",
         y = "Gestational Age(weeks)") + 
    ylim(0,24) + xlim(0,10)
})

wrap_plots(plot_list, ncol = 1)  %>% 
  ggsave(filename = "plots/Suppl_Fig2.png", plot = . ,
         height = 7.5 , width = 6 , dpi = 600)
