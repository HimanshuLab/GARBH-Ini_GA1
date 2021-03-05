#-----------------------------------------------------------------
# Function to return predicted GA values from CRL------------------------------------------------------------
#-----------------------------------------------------------------

predict_ga <- function(crl,
                       method = c("Hadlock",
                                  "McLennan-Schluter",
                                  "Robinson-Fleming",
                                  "Sahota",
                                  "Verburg",
                                  "Intergrowth",
                                  "Garbhini1")){
  ga <- numeric()
  if("Hadlock" %in% method){
    ga <- c(ga,Hadlock = exp(1.684969+0.315646*crl-0.049306*(crl^2)+0.004057*(crl^3)-0.0001204568*(crl^4)))
  }
  if("McLennan-Schluter" %in% method){
    ga <- c(ga,McLennan = (32.61967 + (2.62975*crl*10)-(0.42399*log(crl*10)*(crl*10)))/7)
  }
  if("Robinson-Fleming" %in% method){
    ga <- c(ga,Robinson = (8.052*sqrt(crl*10)+23.73)/7)
  }
  if("Sahota" %in% method){
    ga <- c(ga,Sahota = (26.643+7.822*sqrt(crl*10))/7)
  }
  if("Verburg" %in% method){
    ga <- c(ga,Verburg = exp(1.4653 + 0.001737*crl*10 + 0.2313*log(crl*10)))
  }
  if("Intergrowth" %in% method){
    ga <- c(ga,Intergrowth = (40.9041 + 3.2 * (crl*10)^0.5 + 0.348956*crl*10 )/7)
  }
  if("Garbhini1" %in% method){
    ga <- c(ga,Garbhini1 = 6.73526 + 1.15018*(crl) - 0.02294*(crl^2))
  }
  if("GA_mixed-effects" %in% method){
    ga <- c(ga,GA_mixed_effects = 6.2679 + (crl * 1.53902) - (0.07713 * (crl^2)))
  }
  return(ga)
}

#-----------------------------------------------------------------
# Function to get standard deviation of published formula------------------------------------------------------------
#-----------------------------------------------------------------

get_sd_published_formulas <- function(crl,
                                      method = c("Hadlock",
                                                 "McLennan-Schluter",
                                                 "Robinson-Fleming",
                                                 "Sahota",
                                                 "Verburg",
                                                 "Intergrowth")){
  sd <- NA
  if("Hadlock" %in% method){
    sd <- exp(0.04413)
  } else if("McLennan-Schluter" %in% method){
    sd <- (1.24668+0.00000329*(crl*10)^3)/7
  } else if("Robinson-Fleming" %in% method){
    sd <- 4.7/(2*7)
  } else if("Sahota" %in% method){
    sd <- 3.58/7
  } else if("Verburg" %in% method){
    sd <- exp(0.0459)
  } else if("Intergrowth" %in% method){
    sd <- (2.39102 + 0.0193474*crl*10)/7
  }
  return(sd)
}

#-----------------------------------------------------------------
# Function to return training and validation datasets------------------------------------------------------------
#-----------------------------------------------------------------
`%nin%` <- function(x,y) !(`%in%`(x,y)) # a function for not %in%
library(caret)
return_split_data <- function(data, ratio , seed_value,
                              exclude ="ga_birth",
                              var = "ga"){
  set.seed(seed_value)
  trainIndex <- createDataPartition(data[,var], p = ratio, 
                                    list = FALSE, 
                                    times = 1)
  exclude_training <- data[trainIndex,colnames(data) %in% exclude]
  exclude_validation <- data[-trainIndex,colnames(data) %in% exclude]
  training_set <- data[trainIndex,colnames(data) %nin% exclude]
  validation_set <- data[-trainIndex,colnames(data) %nin% exclude]
  
  rm(trainIndex)
  preProcValues <- preProcess(training_set, method = c("range"))
  
  training_set <- predict(preProcValues, training_set)
  validation_set <- predict(preProcValues, validation_set)
  training_set <- cbind(training_set,exclude_training)
  validation_set <- cbind(validation_set,exclude_validation)
  if(length(exclude) ==1 ){
    colnames(training_set)[colnames(training_set) == "exclude_training"] <- exclude
    colnames(validation_set)[colnames(validation_set) == "exclude_validation"] <- exclude
  }
  return(list(training_set, validation_set))
}

return_split_data <- function(data, ratio , seed_value,
                              exclude ="ga_birth",
                              var = "ga"){
  set.seed(seed_value)
  trainIndex <- createDataPartition(data[,var], p = ratio, 
                                    list = FALSE, 
                                    times = 1)
  exclude_training <- data[trainIndex,colnames(data) %in% exclude]
  exclude_validation <- data[-trainIndex,colnames(data) %in% exclude]
  training_set <- data[trainIndex,colnames(data) %nin% exclude]
  validation_set <- data[-trainIndex,colnames(data) %nin% exclude]
  
  rm(trainIndex)
  preProcValues <- preProcess(training_set, method = c("range"))
  
  training_set <- predict(preProcValues, training_set)
  validation_set <- predict(preProcValues, validation_set)
  training_set <- cbind(training_set,exclude_training)
  validation_set <- cbind(validation_set,exclude_validation)
  if(length(exclude) ==1 ){
    colnames(training_set)[colnames(training_set) == "exclude_training"] <- exclude
    colnames(validation_set)[colnames(validation_set) == "exclude_validation"] <- exclude
  }
  return(list(training=training_set, test=validation_set))
}


# fractional poly formulas ------------------------------------------------


generate_frac_poly_formulas <- function(outcome = "ga" , predictor = "crl"){
  
  c(-2,-1,-0.5,0,0.5,1,2,3) %>% { #possible powers in Fractional Poly models
    c(combn(x = . , m = 2 , simplify = FALSE) , # taking all 2 term combinations
      as.list(.), # adding single terms
      lapply(.[. != 0], function(i) c(i,i))) # adding repeated terms except for 0,0
  }  %>% {
    x <- sapply(. ,function(i){
      if(length(i) == 1){
        paste0("I(",predictor,"^",i,")")
      }else {
        if(i[1] == i[2]){
          paste0("I(",predictor,"^",i[1],") + I(log(",predictor,"^",i[2],"))")
        }else{
          paste0("I(",predictor,"^",i[1],") + I(",predictor,"^",i[2],")")
        }
      }
    }, USE.NAMES = FALSE) %>% 
      paste0(outcome," ~ ",.)
    names(x) <- sapply(. , function(i) {
      paste0("FP(",paste0(i,collapse = ", "),")")
    })
    x
  } %>% 
    gsub(pattern = paste0(predictor,"\\^0)"),
         replacement = paste0("log(",predictor,"))"),
         x = .)  #replacing all predictor^0 with log(predictor)
}
