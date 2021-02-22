
# Library Imports ---------------------------------------------------------
library(magrittr) #used for the pipe operator %>% which feeds the output of the previous function to the next function as . or assigns to first parameter by default 
library(caret)
library(BlandAltmanLeh)
library(DescTools)
library(reshape2)
library(dplyr)
library(readr)
library(svglite)
source("00_Functions.R")

datasets <- c("data_main","data_dbscan","data_clinical","data_test","data_test_dbscan","data_test_clinical")

datasets %>%
  sapply(function(i){
    assign(x = i,
           value = read.table(paste0("data/",i,".tsv"), header = T , sep = "\t",
                              na.strings = c(""," ","NA","88","99",
                                             "99.99","777.77","999",
                                             "77.77","77.7","77")),
           envir = .GlobalEnv)
  }) %>% invisible()

methods =  c("Hadlock",
             "McLennan-Schluter",
             "Robinson-Fleming",
             "Sahota",
             "Verburg",
             "Intergrowth",
             "Garbhini1")
#-----------------------------------------------------------------
# Evaluation of formulae for gest age prediction ------------------------------------------------------------
#-----------------------------------------------------------------

predicted_test <- sapply(data_test$crl, predict_ga) %>% t %>% as.data.frame
predicted_main <- sapply(data_main$crl, predict_ga) %>% t %>% as.data.frame

performance_formulae <- data.frame(R2 = apply(predicted_test,2,R2,data_test$ga),
                                   RMSE = apply(predicted_test,2,RMSE,data_test$ga)) 

predicted_test_dbscan <- sapply(data_test_dbscan$crl, predict_ga) %>% t %>% as.data.frame
predicted_test_clinical <- sapply(data_test_clinical$crl, predict_ga) %>% t %>% as.data.frame

data.frame(R2_test_dbscan = apply(predicted_test_dbscan,2,R2,data_test_dbscan$ga),
           R2_test_clinical = apply(predicted_test_clinical,2,R2,data_test_clinical$ga)) 

performance_formulae %>%
  write_tsv("outputs/Suppl_table_4.tsv")

sapply(colnames(predicted_test), function(i){
  sapply(colnames(predicted_test), function(j){
    mean(predicted_test[,i]-predicted_test[,j])
  })
}) %>% {
  .[lower.tri(.)] <- NA
  .
} %>% 
  unlist %>%  range(na.rm=T)

sapply(methods, function(method){
  pred_test <- predict_ga(crl = data_test$crl, method = method)
  pred_test_dbscan <- predict_ga(crl = data_test_dbscan$crl, method = method)
  pred_test_clinical <- predict_ga(crl = data_test_clinical$crl, method = method)
  
  c(R2(pred_test, data_test$ga) , RMSE(pred_test, data_test$ga),
    R2(pred_test_clinical, data_test_clinical$ga) , RMSE(pred_test_clinical, data_test_clinical$ga),
    R2(pred_test_dbscan, data_test_dbscan$ga) , RMSE(pred_test_dbscan, data_test_dbscan$ga))
}) %>% t()
# generating table 4 ------------------------------------------------------

pairwise_BA_main <- sapply(predicted_main, function(i){
  sapply(predicted_main,function(j){
    BA_stats <- bland.altman.stats(i,j, conf.int = 0.95)
    paste0(BA_stats$mean.diffs %>% round(3),
           "( ",BA_stats$lower.limit %>% round(3),
           ", ",BA_stats$upper.limit %>% round(3),")")
    
  })
})

pairwise_BA_test <- sapply(predicted_test, function(i){
  sapply(predicted_test,function(j){
    BA_stats <- bland.altman.stats(i,j, conf.int = 0.95)
    paste0(BA_stats$mean.diffs %>% round(3),
           "( ",BA_stats$lower.limit %>% round(3),
           ", ",BA_stats$upper.limit %>% round(3),")")
    
  })
})

table4 <- matrix(nrow = ncol(predicted_main), ncol = ncol(predicted_main),
                 dimnames = list(colnames(predicted_main),colnames(predicted_main)))
table4[upper.tri(table4)] <- pairwise_BA_main[upper.tri(pairwise_BA_main)]
table4[lower.tri(table4)] <- pairwise_BA_test[lower.tri(pairwise_BA_test)]
table4 %>% as.data.frame() %>% 
  write_tsv("outputs/Table2.tsv")

#-----------------------------------------------------------------
# Preterm Analysis ------------------------------------------------------------
#-----------------------------------------------------------------

data_main_preterm<- data_main
data_main_preterm$type<-as.factor(data_main_preterm$type)
data_main_preterm[data_main_preterm$type=='ds1','level']<-4
data_main_preterm[data_main_preterm$type=='ds2','level']<-3
data_main_preterm[data_main_preterm$type=='em1','level']<-2
data_main_preterm[data_main_preterm$type=='em2','level']<-1

data_main_preterm <- data_main_preterm[with(data_main_preterm, ave(level, enrid, FUN=max)==level),]

preterm <- data_main_preterm$crl %>% sapply(predict_ga) %>% t() %>% as.data.frame()
preterm <- preterm + (data_main_preterm$ga_birth - data_main_preterm$ga)
preterm$LMP <- data_main_preterm$ga_birth
preterm <- preterm %>% na.omit()
preterm<-preterm[!rowSums(preterm<=20),]
preterm_plot<-preterm[c("Hadlock", "LMP")]
preterm_plot[preterm_plot$Hadlock>=37 & preterm_plot$LMP>=37,'level']<-'Term'
preterm_plot[preterm_plot$Hadlock>=37 & preterm_plot$LMP<37,'level']<-'Preterm, LMP'
preterm_plot[preterm_plot$Hadlock<37 & preterm_plot$LMP>=37,'level']<-'Preterm, Hadlock'
preterm_plot[preterm_plot$Hadlock<37 & preterm_plot$LMP<37,'level']<-'Preterm'


plot_preterm<-ggplot(preterm_plot) +
  geom_point(aes(x=Hadlock, y= LMP, shape=level), size=1, alpha=0.5) +
  geom_hline(yintercept=37, linetype="dashed", size=0.5) +
  geom_vline(xintercept=37, linetype="dashed", size=0.5) +
  scale_colour_manual(values=c("Term" = "darkgreen", "Preterm" = "red", "Preterm, LMP" = "blue", "Preterm, Hadlock" = "orange"))+
  labs(x="GA by Hadlock at birth (weeks)", y="GA by LMP at birth (weeks)", 
       shape="PTB by formula")+
  ggtitle('Classification of preterm birth: LMP-based vs USG (Hadlock)')+
  theme_bw()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(filename = "plots/Figure2c.svg", plot_preterm)
#agreement between formulae -> intersection by union
preterm_class<-preterm
preterm_class[preterm_class<37]<-1
preterm_class[!preterm_class<37]<-0
preterm_vals<-sapply(colnames(preterm), function(i){
  sapply(colnames(preterm), function(j){
    round((nrow(preterm_class[ which(preterm_class[,i] == 1 & preterm_class[,j] == 1) , ])*100)/nrow(preterm_class[ which(preterm_class[,i] == 1 | preterm_class[,j] == 1) , ]),3)
  })
}) %>% data.frame() %>% write_tsv("outputs/Tablee9.tsv")

#sensitivuty, specificity
preterm_new<- as.factor(preterm_class)
preterm_acc<- sapply(c('Hadlock','McLennan','Robinson','Sahota','Verburg','Intergrowth','Garbhini1'), function(i){
  confmat<-confusionMatrix(factor(preterm_class$LMP), factor(preterm_class[,i]),positive='1')
  round(confmat$byClass,3)
} 
                     )
data.frame(Indicators = row.names(preterm_acc),preterm_acc) %>%write_tsv("outputs/Tablee8.tsv")


preterm %>% sapply(function(i) BinomCI(sum(i < 37),length(i), conf.level = 0.95, method = "clopper-pearson")) %>% 
  t() %>% `*`(100) %>% round(2) %>%
  {
    data.frame(formula = row.names(.),
               PTB_rate = .[,1],
               CI = paste0("(", .[,2] , ", " ,
                           .[,3]," )"),
               row.names = NULL)
  } %>% 
  write_tsv("outputs/Table5.tsv")

apply(preterm,2,function(i) sum(i < 37)*100/length(i)) %>% range()

sapply(colnames(preterm), function(i){
  sapply(colnames(preterm), function(j){
    (c((preterm[,i] < 37) %>% table , (preterm[,j] < 37) %>% table) %>%
       matrix(nrow=2) %>% 
       fisher.test())$p.value
  })
}) %>% (function(x) {
  x[lower.tri(x)] <- NA
  x
}) %>%  melt() %>% na.omit() %>% mutate( value = p.adjust(value, method = "bonferroni")) %>%
  filter(Var1 != Var2) %>% transmute( Method1 = Var1,
                                      Method2 = Var2,
                                      `Fisher exact test p-value (Bonferroni corrected)` = value) %>% 
  write_tsv("tables/eTable8.tsv")

(c((preterm$Hadlock < 37) %>% table , (preterm$LMP < 37) %>% table) %>%
    matrix(nrow=2) %>% 
    fisher.test())$p.value

# Supplementary figure 3a -------------------------------------------------

cbind(data_test %>% select(ga, crl) , predicted_test) %>% 
  melt(id.vars = c("ga","crl")) %>% 
  transmute( ga , crl ,
             Formula = variable,
             ga_pred = value) %>% 
  {
    i <-.
    i$Formula <- as.character(i$Formula)
    i$Formula[i$Formula == "Garbhini1"] <- "Garbhini-GA1"
    i
  } %>% 
  ggplot(aes(x = crl )) +
  geom_point(aes(y = ga), alpha = 0.01) +
  geom_line(aes(y = ga_pred , color = Formula ), size = 1.1) +
  theme_bw() + 
  xlim(0,8) + ylim(5,14) +
  labs(caption = paste0("n =", nrow(data_test)),
       x = "Crown Rump Length(cm)",
       y = "Gestational Age(weeks)") +
  theme(legend.title = element_text("Formula"))

cbind(data_test %>% select(ga, crl) , predicted_test) %>% 
  melt(id.vars = c("ga","crl"))  %>% 
  transmute( ga , crl ,
             Formula = variable,
             Error = value - ga) %>% 
  {
    i <-.
    i$Formula <- as.character(i$Formula)
    i$Formula[i$Formula == "Garbhini1"] <- "Garbhini-GA1"
    i
  } %>% 
  ggplot(aes(x = Error, y = Formula)) + 
  geom_violin(fill = "black") +
  theme_bw()
