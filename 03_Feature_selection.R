seed_value <- 666
seed_range <- 100:1099
# Library Imports ---------------------------------------------------------
library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function as . or assigns to first parameter by default 
library(Boruta) #used for the boruta random forests implementation of feature selection
library(caret) #used for createDataPartition() and preProcess()

if(!dir.exists("outputs/feature_selection/")) dir.create("outputs/feature_selection/")

t <- proc.time()
for(i in c("data_main.tsv","data_dbscan.tsv","data_clinical.tsv")){
  
# loading datasets --------------------------------------------------------
 data <- read.table(paste0("data/",i), sep = "\t", header = TRUE,
                     na.strings = c(""," ","NA","88","99",
                                    "99.99","777.77","999",
                                    "77.77","77.7","77"))# setting the various encoded missing values as NA

#-----------------------------------------------------------------
# Pre-Feature Selection Cleaning ------------------------------------------------------------
#-----------------------------------------------------------------
  
  data <- data[, colnames(data) %in% c("ga","crl",
                                       readLines("features_list.txt"))]
  data <- data[,colnames(data) != "enrid"]
  #removing features with a high number of missing values 
  #and then omitting rows with missing values
  data <- data[,apply(data,2,function(i) sum(is.na(i)) < 200)] %>% na.omit()
  #converting categorical variables into factors
  for(j in c("state","rlgn","fmly_typ",
             "part_edu","part_occ","fmly_mem",
             "drnk_wtr","derived_ses_mks_2019",
             "derived_state")){
    data[,j] <- as.factor(data[,j])
  }
  
#-----------------------------------------------------------------
# Feature Selection  ------------------------------------------------------------
#-----------------------------------------------------------------
  
# Random Forests implementation ------------------------------------------------------------------
  
  set.seed(seed_value)
  boruta_output <- Boruta(ga ~ .,
                          data= data,
                          doTrace=0,
                          maxRuns = length(seed_range))
  names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) %>%
    writeLines(paste0("outputs/feature_selection/",sub(".tsv","",i),"_boruta_features.txt"))
  
# linear model implementation ------------------------------------------------------------
  
  features_glm <- lapply(seed_range, function(j){
    set.seed(j)
    # randomly sampling 50% of the data
    trainIndex <- createDataPartition(data$ga, p = 0.5, 
                                      list = FALSE, 
                                      times = 1)
    #normalizing the different features 
    preProcValues <- preProcess(data, method = c("range"))
    training_set <- predict(preProcValues, data[trainIndex,])
    #building the model
    fit_glm <- glm(formula= ga ~. ,family = "gaussian",
                   data =training_set)
    tmp <- summary(fit_glm)$coefficients
    rownames(tmp)[tmp[,4] < 0.05]
  })
  features_glm %>% unlist %>% table %>%
    sort(decreasing = T) %>% {.[. > length(seed_range)/2]} %>% names() %>%
    writeLines(paste0("outputs/feature_selection/",sub(".tsv","",i),"_glm_features.txt"))
}
 t <- proc.time() - t
 print(t)
 

write("Random Forests:",file = "outputs/Suppl_Table_4.txt")
write(readLines("outputs/feature_selection/data_dbscan_boruta_features.txt") %>% paste0("\t",.),
      file = "outputs/Suppl_Table_4.txt", append = TRUE)
write("GLM:",file = "outputs/Suppl_Table_4.txt", append = TRUE)
write(readLines("outputs/feature_selection/data_dbscan_glm_features.txt")%>% paste0("\t",.),
      file = "outputs/Suppl_Table_4.txt", append = TRUE)