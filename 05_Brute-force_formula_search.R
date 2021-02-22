seed_range <- 100:1099
options(stringsAsFactors = FALSE)
# Library Imports ---------------------------------------------------------

library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function as . or assigns to first parameter by default 
library(reshape2) #used for melt() to transform the data from wide to long
library(parallel) #used to run formula fitting in parallel
library(ggplot2)
library(readr)
source("00_Functions.R")
 

# reading union of all features picked during feature selection -----------

features <- list.files("outputs/feature_selection/",full.names = TRUE) %>%
  grep(pattern = "data_dbscan",x = ., value = T) %>%
  sapply(readLines) %>% unlist() %>% unique()
for(i in c("part_edu","drnk_wtr") ){
  features[grepl(i,features)] <- i
}
rm(i)

features[] # preg_num is already included
# reading data_dbscan -----------------------------------------------------

data_dbscan <- read.table("data/data_dbscan.tsv", header = T , sep = "\t",
                          na.strings = c(""," ","NA","88","99",
                                         "99.99","777.77","999",
                                         "77.77","77.7","77"))
data_dbscan <- data_dbscan[,colnames(data_dbscan) %in% c("ga", features)]

#-----------------------------------------------------------------
# formula generation ------------------------------------------------------------
#-----------------------------------------------------------------

# function to generate transformation of features for formula -----------
transform_feature <- function(x){
  return(c("",paste0("+ ",x),
           paste0("+ ",c("log(","sqrt("),x,")")))
}

formulas <- expand.grid("ga ~ " ,
                        lapply(features, transform_feature) %>% 
                          expand.grid(stringsAsFactors = FALSE) %>% apply(1,function(i)paste0(i,collapse="")),
                         c("", ", degree = 2)", ", degree = 3)"),
                        stringsAsFactors = FALSE) %>%
  apply(1,function(i)paste0(i,collapse=""))

formulas <- sub("~ \\+","~",formulas)
formulas[grep("degree", formulas)] <- paste0(gsub("\\+ ",
                                                  ", ",
                                                  formulas[grep("degree", formulas)]))
formulas[grep("degree", formulas)] <- sub("~|~,", "~ polym(",formulas[grep("degree", formulas)] )

#-----------------------------------------------------------------
# Bootstrapping formulas performance over 1000 iterations ------------------------------------------------------------
#-----------------------------------------------------------------
t <- proc.time()
bootstrap <- mclapply(X = seed_range,mc.cores = 40,FUN =  function(seed_value){
  temp <- return_split_data(data = data_dbscan, ratio = 0.70, seed_value = seed_value)
  training_set <- temp[[1]]
  test_set <- temp[[2]]
  ga_pred <- lapply(formulas, function(f){
    tryCatch(expr = predict(lm(formula = as.formula(f), data = training_set), test_set),
             error = function(e) NA,
             warning = function(w) "")
  })
  data.frame(seed_value,
             formula = formulas,
             R2 = sapply(ga_pred, function(i) ifelse(length(i)==1,
                                                     NA,
                                                     R2(unlist(i),               test_set$ga))),
             RMSE = sapply(ga_pred, function(i) ifelse(length(i)==1,
                                                       NA,
                                                       RMSE(unlist(i),test_set$ga))))
  
})
bootstrap <- Reduce(rbind,bootstrap) %>% na.omit
t <- proc.time()-t
print(t)

write.table(bootstrap, file = "outputs/bootstrap_output.tsv", sep = "\t",
            row.names = FALSE)
