# Biblio
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(readr)
library(reshape2)

# Data imports
DE_GOUT_vs_HC <- read_delim("Data for Report Assessment-20221017/DE_GOUT_vs_HC.csv", 
                            delim = "\t", escape_double = FALSE, 
                            trim_ws = TRUE)

DE_SA_vs_HC <- read_delim("Data for Report Assessment-20221017/DE_SA_vs_HC.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

Expression_Table <- read_delim("Data for Report Assessment-20221017/Expression_Table.csv", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)

Sample_Information <- read_delim("Data for Report Assessment-20221017/Sample_Information.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

Gene_BG <- read_delim("Data for Report Assessment-20221017/Annotations.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

# Merging data sets with Gene_BG to get names
## Expression table merged with Gene_BG
ET_annotated <-  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID', )
DE_GOUT_vs_HC <- merge(x= Gene_BG,y= DE_GOUT_vs_HC,by.x= 'ID',by.y = 'ID', )
DE_SA_vs_HC <- merge(x= Gene_BG,y= DE_SA_vs_HC,by.x= 'ID',by.y = 'ID', )

## Removing cols : ID, CHROMOSOME, START, STOP, and BIOTYPE
ET_annotated <-  ET_annotated[-c(3:6)]
DE_GOUT_vs_HC <- DE_GOUT_vs_HC[-c(3:6)]
DE_SA_vs_HC <- DE_SA_vs_HC[-c(3:6)]

# Filtering for Significant genes

SIG_GOUT_v_HC <- DE_GOUT_vs_HC %>%
  filter(p.adj < 0.05) 

SIG_SA_v_HC <- DE_SA_vs_HC %>% 
  filter(p.adj < 0.01) 

# SIG_ET_annotated

SIG_ET_GOUT_v_HC <- Expression_Table %>% 
  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID',) 

SIG_ET_GOUT_v_HC <- SIG_ET_GOUT_v_HC[-c(1,3:6)] %>% 
  filter(SIG_ET_GOUT_v_HC$ID %in% SIG_GOUT_v_HC$ID)

#SIG_ET_GOUT_v_HC[2:28] <- log2(SIG_ET_GOUT_v_HC[2:28] + .25)

SIG_ET_GOUT_v_HC <- setNames(data.frame(t(SIG_ET_GOUT_v_HC[ , - 1])), SIG_ET_GOUT_v_HC[ , 1])

#SIG_ET_GOUT_v_HC <-  merge(SIG_ET_GOUT_v_HC,Sample_Information, by.x = 0, by.y = 'SAMPLE') %>% 
  #relocate('SAMPLE_GROUP':'NEUTROPHILS', .before = 'TF')

# Removing SEPSIS group
SIG_ET_GOUT_v_HC <- SIG_ET_GOUT_v_HC[-c(19:27),]

#names(SIG_ET_GOUT_v_HC)[1] <- "Sample"

#sum_stats of gout

sum_stats_gout <- data.frame(t(sapply(SIG_ET_GOUT_v_HC, function(SIG_ET_GOUT_v_HC) c( "Stand dev" = sd(SIG_ET_GOUT_v_HC), 
                         "Mean"= mean(SIG_ET_GOUT_v_HC,na.rm=TRUE),
                         "n" = length(SIG_ET_GOUT_v_HC),
                         "Median" = median(SIG_ET_GOUT_v_HC),
                         "CoeffofVariation" = sd(SIG_ET_GOUT_v_HC)/mean(SIG_ET_GOUT_v_HC,na.rm=TRUE),
                         "Minimum" = min(SIG_ET_GOUT_v_HC),
                         "Maximun" = max(SIG_ET_GOUT_v_HC),
                         "Upper Quantile" = quantile(SIG_ET_GOUT_v_HC,.95)
)
)))

# Merge with Diff express table and cleaning df
sum_stats_gout <- merge(sum_stats_gout,DE_GOUT_vs_HC, by.x = 0, by.y = 'symbol', )
names(sum_stats_gout)[1] = "Symbol" 
sum_stats_gout <-  sum_stats_gout[-10]

# Filtering for low SD, Coeff of Vara, 

sum_stats_gout <- sum_stats_gout %>% 
  filter(p.adj < 0.01)




