# Library Imports ---------------------------------------------------------
library(haven) # used for read_dta()
library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function 
library(dplyr) #used for mutate()
library(dbscan)
library(ggplot2) #used for plotting
library(cowplot) #used for plot_grid()
source("00_Functions.R") # for predict_ga() containing Garbhini1 formulas
#-----------------------------------------------------------------

options(stringsAsFactors = FALSE)
# Description of variables ------------------------------------------------------------
#-----------------------------------------------------------------

# data_raw: Raw Dataset
# data_preprocessed: Preprocessed Dataset
# data_main: Training Dataset
# data_clinical: Dataset filtered using clinical gold standard criteria
# data_dbscan: Dataset filtered using dbscan clustering to remove noise

#-----------------------------------------------------------------
# Data Processing ---------------------------------------------------------
#-----------------------------------------------------------------

# Raw Dataset read
if(!("data_raw" %in% ls())){
  data_raw <- read_dta("Dataset for 1st tri dating analysis.dta")
  data_raw <- data_raw[!duplicated(data_raw$enrid),]
}

features <- readLines("features_list.txt") # List of features used in this study
filter_parameters <- c("smok_his","tob_chew","bmi","cncv","contr_bcp","brfeed","smok_prs","alch")#additional clinical features

# Preprocessing of the raw dataset
df <- data_raw
# individuals with at least one CRL measurement
crl_present <- df[grepl("^crl_qc_",colnames(df))] %>%
  apply(1,function(i) !all(is.na(i))) 
# individuals for whom the LMP was known
lmp_known <- df$lmp_by == 1 
# individuals with singleton pregnancy
is_singleton <- (df$gest_ds1 == 1 | df$fetal_num_em1 == 1) 
#outcome of the pregnancy was not abortion
not_abortion <- df$op_abort!=1 | is.na(df$op_abort)
# excluding individuals that had a missing value for LMP known or had singleton pregnancy
missing_value <- is.na(lmp_known & is_singleton) 

# selecting individuals with atleast one CRL value present, LMP known, singleton pregnancies, did not results in an abortion
df <- df[which(crl_present & lmp_known & is_singleton & not_abortion & !missing_value),]
#selecting only the columns required for this study
df <- df[, (colnames(df) %in% c(features,filter_parameters)) | grepl("^vdt_scr$|^vdt_(ds|em)\\d$|^crl_qc_(ds11|em)[12]$|^crl_(ga|qa)_(days|wks)_qc_(ds11|em)[12]$|^pog_lmp_|^op_dt_event", colnames(df)) ]

# converting the columns with dates to the appropritate date type
for( i in grep("^vdt_ds\\d$|^vdt_em\\d$|^vdt_scr$",colnames(df),value=T)){
  df[,i] <- as.Date(df[,i] %>% unlist, format = "%Y-%m-%d" , origin = "1970-01-01") 
}
for( i in grep("^op_dt_",colnames(df),value=T)){
  df[,i] <- as.Date(df[,i] %>% unlist, format = "%Y-%m-%d" , origin = "1960-01-01")
}
#converting kuppuswamy scale feature to character from the stata factor type
df$derived_ses_mks_2019 <- df$derived_ses_mks_2019 %>% as_factor() %>% as.character()

rm(i,crl_present,lmp_known,is_singleton,not_abortion,missing_value) #cleaning up variables no longer needed

# data_preprocessed <- df

# calculating the gestational ages and period of gestation from the appropriate dates
df  <- df %>% mutate(date_diff_ds1 =vdt_ds1 - vdt_scr,
                     date_diff_ds2 =vdt_ds2 - vdt_scr,
                     date_diff_em1 =vdt_em1 - vdt_scr,
                     date_diff_em2 =vdt_em2 - vdt_scr,
                     pog_birth=pog_lmp_w+(pog_lmp_d+op_dt_event - vdt_scr)/7,
                     total_pog_lmp_ds1 = pog_lmp_w + (pog_lmp_d + date_diff_ds1)/7,
                     total_pog_lmp_em1 = pog_lmp_w + (pog_lmp_d + date_diff_em1)/7,
                     total_pog_lmp_ds2 = pog_lmp_w + (pog_lmp_d + date_diff_ds2)/7,
                     total_pog_lmp_em2 = pog_lmp_w + (pog_lmp_d + date_diff_em2)/7)
# spliting up each visit into its own row
df <- rbind(data.frame(ga= df$total_pog_lmp_ds1,crl= df$crl_qc_ds111, ga_birth=df$pog_birth, type= "ds1", df[,colnames(df) %in% c(features, filter_parameters)]),
            data.frame(ga= df$total_pog_lmp_ds2,crl= df$crl_qc_ds112, ga_birth=df$pog_birth, type= "ds2" , df[,colnames(df) %in% c(features, filter_parameters)]),
            data.frame(ga= df$total_pog_lmp_em1,crl= df$crl_qc_em1, ga_birth=df$pog_birth, type= "em1" , df[,colnames(df) %in% c(features, filter_parameters)]),
            data.frame(ga= df$total_pog_lmp_em2,crl= df$crl_qc_em2, ga_birth=df$pog_birth, type= "em2" , df[,colnames(df) %in% c(features, filter_parameters)]))
# removing visit that did not occur
df <- df[!is.na(df$ga) & !is.na(df$crl),]
# converting gestational age at visit and at pregnancy to numeric values from date
df$ga <- df$ga %>% as.numeric()
df$ga_birth <- as.numeric(df$ga_birth)

# excluding pateints who had abnormally high GA values at visit or delivery occured extremely preterm
df <- df[df$ga > 0 & df$ga < 30,]
df <- df[df$crl > 0 & df$crl <= 10,]

# Processed Dataset
data_main <- df

# filters used for Dataset1 and Dataset2
no_smoking_history <- df$smok_his ==11
no_chewing_tobacco <-  df$tob_chew == 11
normal_bmi <-  df$bmi >= 18.5 & df$bmi <25
spontaneous_conception <- df$cncv ==11
no_birth_control_use <- df$contr_bcp ==2
breastfeed_prior_conc. <- df$brfeed != 1
no_smoker_presence <- df$smok_prs ==2
no_alcohol_use <- df$alch == 11

dataset1_filters <- no_birth_control_use & spontaneous_conception & normal_bmi & breastfeed_prior_conc.
dataset2_filters <- no_smoker_presence & no_smoking_history & no_chewing_tobacco & no_alcohol_use & normal_bmi

#Dataset 1 excludes individuals prone to error in LMP, Dataset 2 excludes individuals prone to error in CRL
dataset1 <- df[dataset1_filters,]
dataset2 <- df[dataset2_filters,]
rm(dataset1_filters, dataset2_filters)
rm(df)


# Importing the test dataset and processing it ----------------------------

data_test <- read_dta("New dataset for testing dating models (3500- 4500).dta")
data_test <- data_test[!duplicated(data_test$enrid),]
data_test <- data_test[ data_test[grepl("^crl_qc_",colnames(data_test))] %>% apply(1,function(i) !all(is.na(i))),]
#filtering by LMP known , singleton pregnancy ,and no abortion
data_test <- data_test[which(data_test$lmp_by == 1 & (data_test$gest_ds1 == 1 | data_test$fetal_num_em1 == 1) & (data_test$op_abort!=1 | is.na(data_test$op_abort))),]  
#converting required dates columns to date
for( i in grep("^vdt_ds\\d$|^vdt_em\\d$|^vdt_scr$",colnames(data_test),value=T)){
  data_test[,i] <- as.Date(data_test[,i] %>% unlist, format = "%Y-%m-%d" , origin = "1970-01-01")
}
for( i in grep("^op_dt_",colnames(data_test),value=T)){
  data_test[,i] <- as.Date(data_test[,i] %>% unlist, format = "%Y-%m-%d" , origin = "1960-01-01")
}
rm(i)
# calculating the gestational ages and period of gestation from the appropriate dates
data_test <- data_test %>% mutate(date_diff_ds1 =vdt_ds1 - vdt_scr,
                                  date_diff_ds2 =vdt_ds2 - vdt_scr,
                                  date_diff_em1 =vdt_em1 - vdt_scr,
                                  date_diff_em2 =vdt_em2 - vdt_scr,
                                  pog_birth=pog_lmp_w+(pog_lmp_d+op_dt_event - vdt_scr)/7,
                                  total_pog_lmp_ds1 = pog_lmp_w + (pog_lmp_d + date_diff_ds1)/7,
                                  total_pog_lmp_em1 = pog_lmp_w + (pog_lmp_d + date_diff_em1)/7,
                                  total_pog_lmp_ds2 = pog_lmp_w + (pog_lmp_d + date_diff_ds2)/7,
                                  total_pog_lmp_em2 = pog_lmp_w + (pog_lmp_d + date_diff_em2)/7)
# spliting up each visit into its own row
data_test <- rbind(data.frame(ga= data_test$total_pog_lmp_ds1,crl= data_test$crl_qc_DS111, ga_birth=data_test$pog_birth, type= "ds1", data_test[,colnames(data_test) %in% c(features, filter_parameters)]),
                   data.frame(ga= data_test$total_pog_lmp_ds2,crl= data_test$crl_qc_DS112, ga_birth=data_test$pog_birth, type= "ds2" , data_test[,colnames(data_test) %in% c(features, filter_parameters)]),
                   data.frame(ga= data_test$total_pog_lmp_em1,crl= data_test$crl_qc_EM1, ga_birth=data_test$pog_birth, type= "em1" , data_test[,colnames(data_test) %in% c(features, filter_parameters)]),
                   data.frame(ga= data_test$total_pog_lmp_em2,crl= data_test$crl_qc_EM2, ga_birth=data_test$pog_birth, type= "em2" , data_test[,colnames(data_test) %in% c(features, filter_parameters)]))

data_test$ga <- data_test$ga %>% as.numeric()
data_test <- data_test[!is.na(data_test$ga) & !is.na(data_test$crl),]
data_test <- data_test[data_test$ga > 0 & data_test$ga < 30,]
data_test <- data_test[data_test$crl > 0 & data_test$crl <= 10,]

#-----------------------------------------------------------------
# Filtering of Data for Model building -----------------------------------------------------------------
#-----------------------------------------------------------------

# clinical filtering

selected_clinical_filter <- no_smoker_presence & 
  no_chewing_tobacco &
  normal_bmi &
  spontaneous_conception & 
  no_birth_control_use & 
  breastfeed_prior_conc. &
  no_smoker_presence &
  no_alcohol_use
# data points with missing values in the clinical filtering criteria are not selected
selected_clinical_filter[is.na(selected_clinical_filter)] <- FALSE 
data_clinical <- data_main[selected_clinical_filter,]
rm(no_smoker_presence , no_chewing_tobacco , normal_bmi ,spontaneous_conception ,
   no_birth_control_use , breastfeed_prior_conc. ,no_smoking_history, no_alcohol_use)



# Effect of hyperparameters on DBSCAN filtering ---------------------------


suppl_table3 <- sapply(seq(0.2,1,0.1), function(eps){
  sapply(seq(10,40,1), function(minPts){
    clusters <- dbscan(x = data_main[,colnames(data_main) %in% c("crl","ga")],
                       minPts = minPts,
                       eps = eps)$cluster
    paste0(sum(clusters != 0),
           ", ",
           length(unique(clusters))-1)
  })
}) %>% as.data.frame(row.names = paste0("minPts = ",seq(10,40,1)))
colnames(suppl_table3) <- paste0("eps = ",seq(0.2,1,0.1))
write.table(suppl_table3,"tables/Suppl_table_3.tsv")

# DBSCAN filtering

eps <- 0.5 # distance cutoff
minPts <- 20 # minimum number of neighbors
selected_dbscan_filter <- dbscan(x = data_main[,colnames(data_main) %in% c("ga","crl")],
                                 minPts = 20,
                                 eps = 0.5)$cluster %>% {
                                   . != 0 #taking all points not classified as noise
                                 }
rm(eps,minPts)
data_dbscan <- data_main[selected_dbscan_filter,]

#-----------------------------------------------------------------
# Filtering of Test data -----------------------------------------------------------------
#-----------------------------------------------------------------

# clinical filtering ------------------------------------------------------
no_smoking_history <- data_test$smok_his ==11
no_chewing_tobacco <-  data_test$tob_chew == 11
normal_bmi <-  data_test$bmi >= 18.5 & data_test$bmi <25
spontaneous_conception <- data_test$cncv ==11
no_birth_control_use <- data_test$contr_bcp ==2
breastfeed_prior_conc. <- data_test$brfeed != 1
no_smoker_presence <- data_test$smok_prs ==2
no_alcohol_use <- data_test$alch == 11

selected_clinical_filter_test <- no_smoker_presence & 
  no_chewing_tobacco &
  normal_bmi &
  spontaneous_conception & 
  no_birth_control_use & 
  breastfeed_prior_conc. &
  no_smoker_presence &
  no_alcohol_use
# data points with missing values in the clinical filtering criteria are not selected
selected_clinical_filter_test[is.na(selected_clinical_filter_test)] <- FALSE 
data_test_clinical <- data_test[selected_clinical_filter_test,]
rm(no_smoker_presence , no_chewing_tobacco , normal_bmi ,spontaneous_conception ,
   no_birth_control_use , breastfeed_prior_conc. ,no_smoking_history, no_alcohol_use)

# dbscan filtering --------------------------------------------------------

selected_dbscan_filter_test <- dbscan(x = data_test[,colnames(data_test) %in% c("ga","crl")],
                                 minPts = 20,
                                 eps = 0.6)$cluster %>% {
                                   . != 0 #taking all points not classified as noise
                                 }

data_test_dbscan <- data_test[selected_dbscan_filter_test,]


#-----------------------------------------------------------------
# Writing results to files -----------------------------------------------------------------
#-----------------------------------------------------------------

if(!dir.exists("data")) dir.create("data")

grep("^data",ls(), value = TRUE) %>% 
  sapply(function(i){ 
    write.table(x = get(i),
                file = paste0("data/",i,".tsv"),
                sep = "\t",
                row.names = FALSE)
  }) %>% invisible()

#-----------------------------------------------------------------
# Plotting output of filtering  -----------------------------------------------------------------
#-----------------------------------------------------------------
selected_dbscan_filter <- factor(selected_dbscan_filter, levels = c(TRUE,FALSE))
selected_clinical_filter <- factor(selected_clinical_filter, levels = c(TRUE,FALSE))

dbscan_plot <- ggplot(data_main, aes(x=crl ))+ 
  geom_point(aes(shape =selected_dbscan_filter, y= ga)) +
  geom_line(aes(y= predict_ga(crl,method = "Garbhini1")), size = 1.2, color = "gray")+
  scale_shape_manual(values = c(19,21),
                     name = "Selected using DBSCAN filtering",
                     labels = c(paste0("TRUE \n",sum(selected_dbscan_filter == TRUE)),
                                paste0("FALSE \n",sum(selected_dbscan_filter == FALSE)))) +
  theme_bw() +
  theme(legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="Crown Rump Length (cm)", y =" Gestational Age (weeks)")

clinical_plot <- ggplot(data_main, aes(x=crl ))+ 
  geom_point(aes(shape =selected_clinical_filter, y= ga)) +
  geom_line(aes(y= predict_ga(crl,method = "Garbhini1")), size = 1.2, color = "gray") +
  scale_shape_manual(values = c(19,21),
                     name = "Selected using clinical filtering",
                     labels = c(paste0("TRUE \n",sum(selected_clinical_filter == TRUE)),
                                paste0("FALSE \n",sum(selected_clinical_filter == FALSE)))) +
                       theme_bw() +
                       theme(legend.position = "bottom",
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank())+
  labs(x="Crown Rump Length (cm)", y =" Gestational Age (weeks)")

plot_grid(clinical_plot, dbscan_plot, ncol = 2 ,labels = 
            c("A","B")) %>% 
  ggsave(filename = "plots/Figure3.pdf", plot = . ,
         height = 7,
         width = 11)
