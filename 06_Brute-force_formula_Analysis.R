seed_range <- 100:1099
options(stringsAsFactors = FALSE)
# Library Imports ---------------------------------------------------------
library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function as . or assigns to first parameter by default 
library(reshape2) #used for melt() to transform the data from wide to long
library(ggplot2)
library(readr)
library(caret)
# reading union of all features picked during feature selection -----------

features <- list.files("outputs/feature_selection/",full.names = TRUE) %>%
  grep(pattern = "data_main",x = ., value = T) %>%
  sapply(readLines) %>% unlist() %>% unique()
features <- features[features != "derived_parity"] # preg_num is already included

# reading data_dbscan -----------------------------------------------------

data_dbscan <- read.table("data/data_dbscan.tsv", header = T , sep = "\t",
                          na.strings = c(""," ","NA","88","99",
                                         "99.99","777.77","999",
                                         "77.77","77.7","77"))
data_dbscan <- data_dbscan[,colnames(data_dbscan) %in% c("ga", features)]

# reading results of bootsrap model fitting -------------------------------

bootstrap <- read.table("outputs/bootstrap_output.tsv", sep = "\t",
                        header = TRUE)

#-----------------------------------------------------------------
# Evaluating performance of 1000 iterations ------------------------------------------------------------
#-----------------------------------------------------------------

bootstrap_stats <- sapply(unique(bootstrap$formula), function(i){
  c(formula = i ,
    seed_count = sum(bootstrap$formula == i),
    var_count = sum(sapply(features,function(j) grepl(j,i))),
    R2 = mean(bootstrap$R2[bootstrap$formula ==i]),
    RMSE = mean(bootstrap$RMSE[bootstrap$formula ==i]))
}, USE.NAMES = FALSE) %>% t() %>% as.data.frame()

bootstrap_stats$seed_count <- as.integer(bootstrap_stats$seed_count)
bootstrap_stats$var_count <- as.integer(bootstrap_stats$var_count)
bootstrap_stats$R2 <- as.numeric(bootstrap_stats$R2) %>% round(digits = 3)
bootstrap_stats$RMSE <- as.numeric(bootstrap_stats$RMSE)

bootstrap_stats <- bootstrap_stats[bootstrap_stats$seed_count > 0.8*length(seed_range),]

n <- nrow(data_dbscan)
p <- sapply(1:nrow(bootstrap_stats),function(i){
  if(grepl("polym",bootstrap_stats$formula[i])){
    if(grepl("degree = 2", bootstrap_stats$formula[i])){
      (bootstrap_stats$var_count[i] + 1)*bootstrap_stats$var_count[i]/2 #nC2 terms = n(n-1)/2 for polym 2
    }else if(grepl("degree = 3", bootstrap_stats$formula[i])){
      (bootstrap_stats$var_count[i] + 1)*bootstrap_stats$var_count[i]*(bootstrap_stats$var_count[i]-1)/6 #nC3 terms = n(n-1)(n-2)/6 for polym 3
    }else if(grepl("degree = 4", bootstrap_stats$formula[i])){
      (bootstrap_stats$var_count[i] + 1)*bootstrap_stats$var_count[i]*(bootstrap_stats$var_count[i]-1)*(bootstrap_stats$var_count[i]-2)/24 
    }
  }else{
    bootstrap_stats$var_count[i] + 1
  }
})
bootstrap_stats$AdjR2 <- round(1 - (1-bootstrap_stats$R2)*(n-1)/(n-p), digits = 4)
rm(n,p)

bootstrap_stats <- bootstrap_stats[order(bootstrap_stats$AdjR2, decreasing = T),!(colnames(bootstrap_stats) %in% c("seed_count"))]


# plotting on test dataset----------------------------------------------------------------

data_test <- read.table("data/data_test.tsv", header = T , sep = "\t",
                        na.strings = c(""," ","NA","88","99",
                                       "99.99","777.77","999",
                                       "77.77","77.7","77"))

test_formulas <- data.frame(formula = bootstrap_stats$formula[] ,
                            error = t(sapply(bootstrap_stats$formula[], function(formula){
                              model <- lm(as.formula(formula), data = data_dbscan)
                              predict(model, data_test) - data_test$ga
                            }))) 
test_formulas <- melt(test_formulas, id.vars = "formula")
colnames(test_formulas) <- c("formula","sample","error")
plot <- ggplot(test_formulas, aes(y=error, x= formula))+
  geom_violin(fill = "black")+ggtitle("Distribution across formulas")+
  scale_y_continuous(breaks = -8:8, limits = c(-8,8))+ theme_bw()+
  theme(legend.position = "none") + coord_flip()
ggsave(filename = "plots/Suppl_Fig2.svg", 
       plot = plot,
       height = 12,
       width = 7)

write_tsv(bootstrap_stats, "tables/Suppl_table_6.tsv")

formula_stats <- data.frame(formula = bootstrap_stats$formula,
                            R2 = sapply(bootstrap_stats$formula,function(j){
                              predict(lm(as.formula(j),data_dbscan),data_test)%>%
                                R2(data_test$ga, na.rm = T)
                            }),
                            RMSE = sapply(bootstrap_stats$formula,function(j){
                              predict(lm(as.formula(j),data_dbscan),data_test)%>%
                                RMSE(data_test$ga, na.rm = T)
                            }),
                            var_count = sapply(bootstrap_stats$formula,function(i){
                              sum(sapply(features,function(j) grepl(j,i)) %>% unlist())
                            }))
n <- nrow(data_test)
p <- sapply(1:nrow(formula_stats),function(i){
  if(grepl("polym",formula_stats$formula[i])){
    if(grepl("degree = 2", formula_stats$formula[i])){
      (formula_stats$var_count[i] + 1)*formula_stats$var_count[i]/2 #nC2 terms = n(n-1)/2 for polym 2
    }else if(grepl("degree = 3", formula_stats$formula[i])){
      (formula_stats$var_count[i] + 1)*formula_stats$var_count[i]*(formula_stats$var_count[i]-1)/6 #nC3 terms = n(n-1)(n-2)/6 for polym 3
    }else if(grepl("degree = 4", bootstrap_stats$formula[i])){
      (formula_stats$var_count[i] + 1)*formula_stats$var_count[i]*(formula_stats$var_count[i]-1)*(formula_stats$var_count[i]-2)/24 
    }
  }else{
    formula_stats$var_count[i] + 1
  }
})
formula_stats$AdjR2 <- round(1 - (1-formula_stats$R2)*(n-1)/(n-p), digits = 4)

write_tsv(formula_stats, "tables/Suppl_table_6_test.tsv")
