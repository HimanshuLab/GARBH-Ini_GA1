
# Library Imports ---------------------------------------------------------
library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function as . or assigns to first parameter by default 
library(ggplot2)
library(dplyr)
library(BlandAltmanLeh) #used for bland.altman.stats()
library(cowplot)
library(ready)
source("00_Functions.R")

# custom function definitions ---------------------------------------------

return_median_IQR <- function(x){
  paste0("\t",median(x, na.rm = T) %>% round(digits = 2),
         "\t(", quantile(x,na.rm = T)[2] %>% round(digits = 2),
         ",",quantile(x,na.rm = T)[4] %>% round(digits = 2),")")
}
return_mean_sd <- function(x){
  paste0("\t",mean(x, na.rm = T) %>% round(digits = 3),"±", sd(x, na.rm = T) %>% round(digits = 2))
}
return_percentages <- function(x){
  (table(x)*100/length(x)) %>% round(3) %>% {
    paste0("\t",names(.),"\t",.,"%", collapse = "\n")
  } 
}
# reading main data, dataset1 and dataset2-------------------------------------------------------
c("data_raw","dataset1","dataset2","data_main") %>%
  sapply(function(i){
    assign(x = i,
           value = read.table(paste0("data/",i,".tsv"), header = T , sep = "\t",
                              na.strings = c(""," ","NA","88","99",
                                             "99.99","777.77","999",
                                             "77.77","77.7","77")),
           envir = .GlobalEnv)
  }) %>% invisible()

#-----------------------------------------------------------------
# Description of participants ------------------------------------------------------------
#-----------------------------------------------------------------

sink("outputs/3.1_descriptions.txt",type = "output",split = TRUE)
cat("sample size",nrow(data_main),"\n")
cat("Age",return_median_IQR(data_raw$derived_age_approx_age),"\n")
cat("Parity (primigravida = 0)\n",return_percentages(data_raw$derived_parity),"\n")
cat("Height",return_median_IQR(data_raw$hght),"\n")
cat("Weight",return_median_IQR(data_raw$anc_cur_wt),"\n")
cat("BMI",return_median_IQR(data_raw$bmi),"\n")
cat("Family type\n",return_percentages(data_raw$fmly_typ),"\n")
cat("Gestation Age\n",return_median_IQR(data_main$ga),"\n")
sink()

#-----------------------------------------------------------------
# Comparison of USG-Hadlock and LMP-based dating methods ------------------------------------------------------------
#-----------------------------------------------------------------

sink("outputs/3.2_USGvsLMP.txt",type = "output",split = TRUE)
#Bland-Altman Analysis output
BA_stats <- bland.altman.stats(predict_ga(data_main$crl,method = "Hadlock"), data_main$ga, conf.int = 0.95)
BA_stats_dataset1 <- bland.altman.stats(predict_ga(dataset1$crl,method = "Hadlock"), dataset1$ga, conf.int = 0.95)
BA_stats_dataset2 <- bland.altman.stats(predict_ga(dataset2$crl,method = "Hadlock"), dataset2$ga, conf.int = 0.95)

#The fit of regression line between means and differences
BA_fit <- data.frame( means = BA_stats$means , diffs = BA_stats$diffs) %>% lm(diffs ~ means, data = .)
cat("Bland-Altman Statistics\n")
cat("USG-LMP ",return_mean_sd(BA_stats$diffs),"\n")
cat("Limits of Agreement",BA_stats$lower.limit,BA_stats$upper.limit,"\n")
cat("Precentage outside:",
    sum(BA_stats$diffs < BA_stats$lower.limit |
          BA_stats$diffs > BA_stats$upper.limit)*100/length(BA_stats$diffs),"\n")
cat("USG-LMP dataset1",return_mean_sd(BA_stats_dataset1$diffs),"\n")
cat("Limits of Agreement dataset1",BA_stats_dataset1$lower.limit,BA_stats_dataset1$upper.limit,"\n")
cat("Precentage outside dataset1:",
    sum(BA_stats_dataset1$diffs < BA_stats_dataset1$lower.limit |
          BA_stats_dataset1$diffs > BA_stats_dataset1$upper.limit)*100/length(BA_stats_dataset1$diffs),"\n")
cat("USG-LMP dataset1",return_mean_sd(BA_stats$diffs[data_main$ga > 11 & data_main$ga < 14]),"\n")
cat("Limits of Agreement dataset1",BA_stats_dataset2$lower.limit,BA_stats_dataset2$upper.limit,"\n")
cat("Precentage outside  dataset1:",
    sum(BA_stats_dataset2$diffs < BA_stats_dataset2$lower.limit |
          BA_stats_dataset2$diffs > BA_stats_dataset2$upper.limit)*100/length(BA_stats_dataset2$diffs),"\n")
sink()

# analysis on dataset1 ----------------------------------------------------
bland.altman.stats(predict_ga(dataset1$crl,method = "Hadlock"), dataset1$ga, conf.int = 0.95) %>% 
{
  cat("USG-LMP datasetl ",return_mean_sd(.$diffs),"\n")
  cat("Limits of Agreement dataset1",.$lower.limit,.$upper.limit,"\n")
} 

# analysis on dataset2 ----------------------------------------------------
bland.altman.stats(predict_ga(dataset2$crl,method = "Hadlock"), dataset2$ga, conf.int = 0.95) %>% 
{
  cat("USG-LMP dataset2 ",return_mean_sd(.$diffs),"\n")
  cat("Limits of Agreement dataset2",.$lower.limit,.$upper.limit,"\n")
}

# plotting figure 2ab -----------------------------------------------------

fig2a <- data.frame(diff = predict_ga(data_main$crl, method = "Hadlock") - data_main$ga) %>% 
  ggplot(aes(x=diff))+ 
  geom_histogram(binwidth = 0.5, color = "black", fill = "white")+
  annotate(geom="text",x=6,y=350,label = paste0("n = ",nrow(data_main) ))+
  theme_bw()+
  labs(title = "Distribution of difference between USG- & LMP-based GA",
       x= "Difference in GA", y = "Count")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

fig2b <- data.frame( means = BA_stats$means , diffs = BA_stats$diffs , fit = BA_fit$fitted.values) %>%
  ggplot(aes(x= means)) + geom_point(aes(y= diffs), alpha = 0.7) + 
  geom_line(aes(y = fit), linetype = "dashed", size = 1.2)+
  annotate(geom="text",x=18,y=8,label = paste0("n = ",length(BA_stats$means)))+
  geom_ribbon(aes(ymin = fit - sd(BA_stats$diffs)* qnorm(1-0.05/2),
                  ymax = fit + sd(BA_stats$diffs)* qnorm(1-0.05/2)),
              alpha = 0.3) + 
  labs(x= "Mean of Hadlock and LMP-based GA",
       y = "Hadlock-based GA - LMP-based GA") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
plot_grid(fig2a,fig2b, ncol = 1 , labels = c("A","B")) %>%
ggsave(filename = "plots/Figure2.pdf",plot = .,
       height = 12, width = 9)

rm(BA_stats, BA_fit)

#-----------------------------------------------------------------
# Writing the results of Table 3 ------------------------------------------------------------
#-----------------------------------------------------------------

file.create("tables/Table1_data.txt")
write(x = "#Age", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_median_IQR(data_raw$derived_age_approx_age),file = "tables/Table1_data.txt" ,
      append = T )

write(x = "#GA estational at enrolment by LMP ", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_mean_sd(data_main$ga),file = "tables/Table1_data.txt" ,
      append = T )

write(x = "#GA at enrolment by USG-Hadlock ", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_mean_sd(predict_ga(data_main$crl, method = "Hadlock")),file = "tables/Table1_data.txt" ,
      append = T )

write(x = "#BMI ", file = "tables/Table1_data.txt" ,
      append = T)
data_raw$bmi %>% sapply(function(i){
  if(is.na(i)){
    NA
  }else if(i < 18.5){
    "Underweight"
  }else if(i >= 18.5 & i < 25){
    "Normal"
  }else if( i >= 25 & i < 30){
    "Obese"
  } else "Overweight"
})  %>% return_percentages() %>%
  write(x= . ,file = "tables/Table1_data.txt" ,
        append = T )

write(x = "#Haemoglobin  ", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_median_IQR(data_raw$hemglo),file = "tables/Table1_data.txt" ,
      append = T )

write(x = "#Height   ", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_median_IQR(data_raw$hght),file = "tables/Table1_data.txt" ,
      append = T )

write(x = "#Socioeconomic Status ", file = "tables/Table1_data.txt" ,
      append = T)
data_raw$derived_ses_mks_2019 %>%  as.character %>% 
{ 
  .[is.na(.)] <- "Undetermined"
  .
} %>% return_percentages() %>%
  write(x= .,file = "tables/Table1_data.txt" ,
        append = T )

write(x = "# Parity", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_percentages(data_raw$derived_parity),file = "tables/Table1_data.txt" ,
      append = T )

write(x = "# Education", file = "tables/Table1_data.txt" ,
      append = T)
data_raw$part_edu %>% sapply(function(i){
  if(is.na(i)){
    NA
  }else if ( i == 11){
    "Illiterate"
  }else if ( i == 12){
    "Literate or primary school  "
  }else if ( i == 13){
    "Middle school "
  }else if ( i == 14){
    "High school "
  }else if ( i == 15){
    "Post high school diploma "
  }else if ( i == 16){
    "Graduate "
  }else if ( i == 17){
    "Post-graduate "
  }
}) %>% return_percentages() %>% 
  write(x= .,file = "tables/Table1_data.txt" ,
        append = T )

write(x = "# Occupation", file = "tables/Table1_data.txt" ,
      append = T)
data_raw$part_occ %>% sapply(function(i){
  if(is.na(i)){
    NA
  }else if ( i == 11){
    "Unemployed"
  }else if ( i == 12){
    "Unskilled worker "
  }else if ( i == 13){
    "Semi-skilled worker "
  }else if ( i == 14){
    "Skilled worker "
  }else if ( i == 15){
    "Clerk, shop, farm owner "
  }else if ( i == 16){
    "Semi-professional "
  }else if ( i == 17){
    "Professional "
  }
}) %>% return_percentages() %>% 
  write(x= .,file = "tables/Table1_data.txt" ,
        append = T )

write(x = "# Religion", file = "tables/Table1_data.txt" ,
      append = T)
data_raw$rlgn %>% sapply(function(i){
  if(is.na(i)){
    NA
  }else if ( i == 11){
    "Hindu"
  }else if ( i == 12){
    "Muslim"
  }else if ( i == 13){
    "Sikh"
  }else if ( i == 14){
    "Christian"
  }else if ( i == 15){
    "Buddhist"
  }else if ( i == 16){
    "Other"
  }else if ( i == 17){
    "More than one"
  }
}) %>% return_percentages() %>% 
  write(x= .,file = "tables/Table1_data.txt" ,
        append = T )

write(x = "#Secondhand tobacco smoke", file = "tables/Table1_data.txt" ,
      append = T)
data_raw$smok_prs %>% sapply(function(i){
  if(is.na(i)){
    "Undertermined"
  }else if(i == 1){
    "Exposed"
  }else if(i ==2){
    "Unexposed"
  }
}) %>% return_percentages() %>%
  write(x= .,file = "tables/Table1_data.txt" ,
        append = T )

data_raw %>% 
  mutate(safe_fuel = ifelse(fuel == 11, 1,0),
         safe_water = ifelse(drnk_wtr == 11 | drnk_wtr == 19, 1, 0),
         chro_ill = ifelse(htn==1 | diab==1 |hypothy==1|hyprthy==1,1,0),
         past_htp = if_else(prclm==1|(swld>30 & high_bp==1),1,0,missing = 0)) %>% 
  select(safe_fuel, safe_water, chro_ill, past_htp) %>% colMeans(na.rm = TRUE) %>%
  round(4) %>% {data.frame(col = names(.), val = ., stringsAsFactors = FALSE)} %>%
  apply(1,function(i){
    write(x = paste0("# ",i[1],"\n\t","TRUE = ", 100*as.numeric(i[2]) ,"\n\tFALSE = ",100*(1-as.numeric(i[2]))),
          file = "tables/Table1_data.txt" ,
          append = T )
  }) %>% invisible()

write(x = "# Contraceptive History", file = "tables/Table1_data.txt" ,
      append = T)
write(x= return_percentages(data_main$contr_bcp),file = "tables/Table1_data.txt" ,
      append = T )

plot <- data_main %>% 
  transmute(crl = crl,
            GA_LMP = ga_birth,
            GA_hd = (ga_birth - ga) +predict_ga(crl, method = "Hadlock"),
            color = paste0(GA_LMP < 37, GA_hd <37)) %>% 
  filter(!(is.na(GA_hd)|is.na(GA_LMP))) %>% 
  ggplot(aes(x = GA_hd, y = GA_LMP))+
  geom_point(aes(color = color))+
  geom_hline(yintercept = 37,linetype = 2)+
  geom_vline(xintercept = 37,linetype = 2)+
  xlim(20,50)+ ylim(20,50)+
  scale_color_manual(values = c("darkgreen","blue","#8B008B","darkred"))+
  theme_bw() +
  labs( x = "GA by Hadlock at birth(weeks)",
        y = "GA by LMP at birth(weeks)")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(filename = "plots/Figure2c.pdf",
       plot = plot,
       height = 7, width = 7)


# Supplementary table S10 -------------------------------------------------

x <- data_main %>% mutate( ga = as.integer(ga)) %>% 
  group_by(ga) %>% 
  summarise(n = n())

colnames(x) <- c("Gestational Age", 'count')
write_tsv(x, file= 'tables/TableS11.tsv')
