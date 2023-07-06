                            ###//// IMPORTS\\\\###
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(readr)
library(reshape2)
library(logistf)
library(ggpubr)
library(gt)

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
                    
                        ###//// DE_GOUT_vs_SA table creation \\\\###

DE_GOUT_vs_SA <-  as.data.frame(matrix(0,ncol = 2,nrow = nrow(Expression_Table)))
names(DE_GOUT_vs_SA) = c('log2Fold','P')
DE_GOUT_vs_SA <- DE_GOUT_vs_SA %>% 
  add_column(Expression_Table[1],.before = 1)

# get the SA and Gout expression tables
# split age into 2 groups

Sample_Information <- Sample_Information %>% 
  remove_rownames() %>% 
  column_to_rownames(var = 'SAMPLE') 

samples_sa = row.names(subset(Sample_Information, SAMPLE_GROUP == "SEPSIS"))
samples_gout = row.names(subset(Sample_Information, SAMPLE_GROUP == "GOUT"))
em_sa = Expression_Table[,samples_sa]
em_gout = Expression_Table[,samples_gout]

# start the loop - nore this is simply the code from "test a gene" with minor modifications.
for (row in 1:nrow(Expression_Table))
{
  gene_data_sa = as.numeric(em_sa[row,])
  gene_data_gout = as.numeric(em_gout[row,])
  mean_sa = mean(gene_data_sa)
  mean_gout = mean(gene_data_gout)
  log2fold = log2(mean_gout) - log2(mean_sa)
  p = t.test(gene_data_gout,gene_data_sa)
  p = p$p.value
  
  DE_GOUT_vs_SA[row,"log2Fold"] = log2fold
  DE_GOUT_vs_SA[row,"P"] = p
  DE_GOUT_vs_SA[row,"p.adj"] = p.adjust(p, n = nrow(Expression_Table))
}

                      ###//// CLEANINNG and ANNOTATING DATA FRAMES \\\\###

# Removing duplicate symbol values
Gene_BG <- Gene_BG[!duplicated(Gene_BG$symbol),]

# Merging data sets with Gene_BG to get names
## Expression table merged with Gene_BG
ET_annotated <-  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID', )
DE_GOUT_vs_HC <- merge(x= Gene_BG,y= DE_GOUT_vs_HC,by.x= 'ID',by.y = 'ID', )
DE_SA_vs_HC <- merge(x= Gene_BG,y= DE_SA_vs_HC,by.x= 'ID',by.y = 'ID', )
DE_GOUT_vs_SA<- merge(x= Gene_BG,y= DE_GOUT_vs_SA,by.x= 'ID',by.y = 'ID', )

## Removing cols : ID, CHROMOSOME, START, STOP, and BIOTYPE
ET_annotated <-  ET_annotated[-c(3:6)]
DE_GOUT_vs_HC <- DE_GOUT_vs_HC[-c(3:6)]
DE_SA_vs_HC <- DE_SA_vs_HC[-c(3:6)]
DE_GOUT_vs_SA <- DE_GOUT_vs_SA[-c(3:6)]

#Filter for significance and absolute fold change
SIG_GOUT_v_HC <- DE_GOUT_vs_HC %>%
  filter(p.adj < 0.05) %>% 
  filter(abs(log2Fold) > 1)

SIG_SA_v_HC <- DE_SA_vs_HC %>% 
  filter(p.adj < 0.05) %>% 
  filter(abs(log2Fold) > 2)

SIG_GOUT_v_SA <- DE_GOUT_vs_SA %>% 
  filter(p.adj < 0.05) %>% 
  filter(abs(log2Fold) > 2)
# SA is different because it had a lot of significant gene for greater filtering

# SIG_ET_GOUT_annotated

SIG_ET_GOUT_v_HC <- Expression_Table %>% 
  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID',) 

SIG_ET_GOUT_v_HC <- SIG_ET_GOUT_v_HC[-c(1,3:6)] %>% 
  filter(SIG_ET_GOUT_v_HC$ID %in% SIG_GOUT_v_HC$ID)

#SIG_ET_GOUT_v_HC[2:28] <- log2(SIG_ET_GOUT_v_HC[2:28] + .25)

SIG_ET_GOUT_v_HC <- setNames(data.frame(t(SIG_ET_GOUT_v_HC[ , - 1])), SIG_ET_GOUT_v_HC[ , 1])

#Removing SEPSIS group
SIG_ET_GOUT_v_HC <- SIG_ET_GOUT_v_HC[-c(19:27),]


# SIG_ET_SA_annotated

SIG_ET_SA_v_HC <- Expression_Table %>% 
  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID',) 

SIG_ET_SA_v_HC <- SIG_ET_SA_v_HC[-c(1,3:6)] %>% 
  filter(SIG_ET_SA_v_HC$ID %in% SIG_SA_v_HC$ID)

#SIG_ET_SA_v_HC[2:28] <- log2(SIG_ET_SA_v_HC[2:28] + .25)

SIG_ET_SA_v_HC <- setNames(data.frame(t(SIG_ET_SA_v_HC[ , - 1])), SIG_ET_SA_v_HC[ , 1])

#Removing GOUT group
SIG_ET_SA_v_HC <- SIG_ET_SA_v_HC[-c(10:18),]

# SIG_ET_GOUT_SA_annotated

SIG_ET_SA_v_GOUT <- Expression_Table %>% 
  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID',) 

SIG_ET_SA_v_GOUT <- SIG_ET_SA_v_GOUT[-c(1,3:6)] %>% 
  filter(SIG_ET_SA_v_GOUT$ID %in% SIG_GOUT_v_SA$ID)

#SIG_ET_SA_v_HC[2:28] <- log2(SIG_ET_SA_v_HC[2:28] + .25)

SIG_ET_SA_v_GOUT <- setNames(data.frame(t(SIG_ET_SA_v_GOUT[ , - 1])), SIG_ET_SA_v_GOUT[ , 1])

#Removing HC group
SIG_ET_SA_v_GOUT <- SIG_ET_SA_v_GOUT[-c(1:9),]




                        ###//// SUMMARY STATISTICS TABLES \\\\###

# Separating groups
SIG_ET_GOUT <- SIG_ET_GOUT_v_HC[c(10:18),]
SIG_ET_SA   <- SIG_ET_SA_v_HC[c(10:18),]
SIG_ET_GOUT2 <- SIG_ET_SA_v_GOUT[c(10:18),]

# Column range for numeric conversion
numeric_cols <- c(2:9)

#sum_stats of Gout
sum_stats_gout<- data.frame(t(sapply(SIG_ET_GOUT, function(SIG_ET_GOUT) c(
                                                                           "Minimum" = min(SIG_ET_GOUT),
                                                                           "Mean"= mean(SIG_ET_GOUT,na.rm=TRUE),
                                                                           "Median" = median(SIG_ET_GOUT),
                                                                           "Maximun" = max(SIG_ET_GOUT),
                                                                           "Standard Deviation" = sd(SIG_ET_GOUT), 
                                                                           "Coefficient of Variation" = sd(SIG_ET_GOUT)/mean(SIG_ET_GOUT,
                                                                                                                             na.rm=TRUE))

)))

# Merge with Diff express table and cleaning df
sum_stats_gout <- merge(sum_stats_gout,DE_GOUT_vs_HC, by.x = 0, by.y = 'symbol', )
names(sum_stats_gout)[1] = "Symbol" 
sum_stats_gout <-  sum_stats_gout[-8]
sum_stats_gout <-  sum_stats_gout %>% mutate_at(vars(numeric_cols), as.numeric)

#sum_stats of SA
sum_stats_sa<- data.frame(t(sapply(SIG_ET_SA, function(SIG_ET_SA) c( 
                                                                     "Minimum" = min(SIG_ET_SA),
                                                                     "Mean"= mean(SIG_ET_SA,na.rm=TRUE),
                                                                     "Median" = median(SIG_ET_SA),
                                                                     "Maximun" = max(SIG_ET_SA),
                                                                     "Standard Deviation" = sd(SIG_ET_SA), 
                                                                     "Coefficient of Variation" = sd(SIG_ET_SA)/mean(SIG_ET_SA,
                                                                                                                           na.rm=TRUE))

)))

# Merge with Diff express table and cleaning df
sum_stats_sa <- merge(sum_stats_sa,DE_SA_vs_HC, by.x = 0, by.y = 'symbol', )
names(sum_stats_sa)[1] = "Symbol" 
sum_stats_sa <-  sum_stats_sa[-8]
sum_stats_sa <- sum_stats_sa %>% 
                mutate_at(vars(numeric_cols), as.numeric)%>% 
                filter(Coefficient.of.Variation < 0.2)

#sum_stats of SA v GOUT
sum_stats_sa_v_gout <- data.frame(t(sapply(SIG_ET_GOUT2, function(SIG_ET_GOUT2) c( 
  "Minimum" = min(SIG_ET_GOUT2),
  "Mean"= mean(SIG_ET_GOUT2,na.rm=TRUE),
  "Median" = median(SIG_ET_GOUT2),
  "Maximun" = max(SIG_ET_GOUT2),
  "Standard Deviation" = sd(SIG_ET_GOUT2), 
  "Coefficient of Variation" = sd(SIG_ET_GOUT2)/mean(SIG_ET_GOUT2,
                                                  na.rm=TRUE))
  
)))

# Merge with Diff express table and cleaning df
sum_stats_sa_v_gout <- merge(sum_stats_sa_v_gout,DE_SA_vs_HC, by.x = 0, by.y = 'symbol', )
names(sum_stats_sa_v_gout)[1] = "Symbol" 
sum_stats_sa_v_gout <-  sum_stats_sa_v_gout[-8]
sum_stats_sa_v_gout <- sum_stats_sa_v_gout %>% 
  mutate_at(vars(numeric_cols), as.numeric)%>% 
  filter(Coefficient.of.Variation < 0.2)


# -------------------------------------------------------- #


# Adding Sample information to SIG gene GOUT v HC
SIG_ET_GOUT_v_HC <-  merge(SIG_ET_GOUT_v_HC,Sample_Information, by.x = 0, by.y = 'SAMPLE') %>% 
  relocate('SAMPLE_GROUP':'NEUTROPHILS', .after = 1)
names(SIG_ET_GOUT_v_HC)[1] <- "Sample"

SIG_ET_GOUT_v_HC <- SIG_ET_GOUT_v_HC %>% 
  mutate(SAMPLE_GROUP = recode(SAMPLE_GROUP, 'HC' = 0, 'GOUT' = 1))


# -------------------------------------------------------- #
# Resetting SIG_ET_SA_annotated using names from sum_stats_sa
SIG_ET_SA_v_HC <- Expression_Table %>% 
  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID',) 

SIG_ET_SA_v_HC <- SIG_ET_SA_v_HC[-c(1,3:6)] %>% 
  filter(SIG_ET_SA_v_HC$symbol %in% sum_stats_sa$Symbol)

SIG_ET_SA_v_HC <- setNames(data.frame(t(SIG_ET_SA_v_HC[ , - 1])), SIG_ET_SA_v_HC[ , 1])

#Removing GOUT group
SIG_ET_SA_v_HC <- SIG_ET_SA_v_HC[-c(10:18),]


# Resetting Sample information
Sample_Information <- read_delim("Data for Report Assessment-20221017/Sample_Information.csv", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

# Adding Sample information to SIG gene SA v HC
SIG_ET_SA_v_HC <-  merge(SIG_ET_SA_v_HC,Sample_Information, by.x = 0, by.y = 'SAMPLE') %>% 
  relocate('SAMPLE_GROUP':'NEUTROPHILS', .after = 1)
names(SIG_ET_SA_v_HC)[1] <- "Sample"

SIG_ET_SA_v_HC <- SIG_ET_SA_v_HC %>% 
  mutate(SAMPLE_GROUP = recode(SAMPLE_GROUP, 'HC' = 0, 'SEPSIS' = 1))

# -------------------------------------------------------- #

# Resetting SIG_ET_SA_V_GOUT annotated using names from sum_stats_sa_v_gout

SIG_ET_SA_v_GOUT <- Expression_Table %>% 
  merge(x= Gene_BG,y= Expression_Table,by.x= 'ID',by.y = 'ID',) 

SIG_ET_SA_v_GOUT <- SIG_ET_SA_v_GOUT[-c(1,3:6)] %>% 
  filter(SIG_ET_SA_v_GOUT$symbol %in% sum_stats_sa_v_gout$Symbol)

SIG_ET_SA_v_GOUT <- setNames(data.frame(t(SIG_ET_SA_v_GOUT[ , - 1])), SIG_ET_SA_v_GOUT[ , 1])

#Removing HC group
SIG_ET_SA_v_GOUT <- SIG_ET_SA_v_GOUT[-c(1:9),]

# Adding Sample information to SIG gene SA v GOUT
SIG_ET_SA_v_GOUT <-  merge(SIG_ET_SA_v_GOUT,Sample_Information, by.x = 0, by.y = 'SAMPLE') %>% 
  relocate('SAMPLE_GROUP':'NEUTROPHILS', .after = 1)
names(SIG_ET_SA_v_GOUT)[1] <- "Sample"








                            ###//// ANALYSIS OF NEUTROPHILS \\\\###

# Distribution of Neutrophil counts among sample group in clinical information
neutrophil_distribution <- ggplot(Sample_Information, aes(fill = SAMPLE_GROUP, x = NEUTROPHILS)) + 
                            geom_density(alpha = 0.5) +
                          labs(x = "Neutrophil Count", y = "Density",
                          title = "Distribution of neutrophil amongst groups") +
                          theme_pubclean() +
                          labs_pubr() 
neutrophil_distribution

# lm of sample group and neutrophil count
Sample_Information$SAMPLE_GROUP <- as.factor(Sample_Information$SAMPLE_GROUP)
neutro_model <- lm(Sample_Information$NEUTROPHILS ~ 0 + Sample_Information$SAMPLE_GROUP)
summary(neutro_model)
anova(neutro_model)




                          ###//// ANALYSIS OF SIG GENES IN GOUT vs HC \\\\###

# Melted SIG ET GOUT table for dot plots

melted_SIG_ET_GOUT_v_HC <- melt(SIG_ET_GOUT_v_HC,id.vars='SAMPLE_GROUP', measure.vars=c(4:19))
melted_SIG_ET_GOUT_v_HC <- melted_SIG_ET_GOUT_v_HC %>% 
  mutate(SAMPLE_GROUP = recode(SAMPLE_GROUP, '0' = 'HC', '1' ='GOUT'))

# Box plot exploration of data

Gout_boxplot <- ggplot(melted_SIG_ET_GOUT_v_HC, aes(x= as.factor(SAMPLE_GROUP), y=value)) +
  facet_wrap(~variable, scales = "free") +
  stat_compare_means(method = "t.test", label.x.npc = 'left', label.y.npc = "top" ,
                     label = "p.format", vjust = 1.5) +
  labs(x = NULL, y = "Gene Expression",
       title = "Significant differential gene expression  and  neutrophil count in Gout patients",
       subtitle = "P values derived from unpaired T-test") +
  theme_pubclean() +
  labs_pubr() +
  geom_boxplot(width = 0.3)

Gout_boxplot

# Box plot exploration of data

Gout_dotplot <- ggplot(melted_SIG_ET_GOUT_v_HC,aes(y= SAMPLE_GROUP, x = value)) + 
                geom_point() +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE) +
                facet_wrap(~variable, scales = 'free')
Gout_dotplot


# Binomial GLM
binomial_gout_model <- glm(SAMPLE_GROUP ~ PCP2,
                           data = SIG_ET_GOUT_v_HC, family = binomial)

# This model can't be fitted due to low sample size causing a Hauck-Donner phenomenon 
# due to low sample size

#"MYO3B","TF","SPP1","SULT4A1","GATD3A","PCP2","KLHDC7A","EGFL6" ,"SNORA73B","AC116407.1","RPL35P5","SMKR1","SNHG25","Z82217.1", AC010300.1


summary(binomial_gout_model)
anova(binomial_gout_model)
hist(rstandard(binomial_gout_model))
autoplot(binomial_gout_model)

# Volcano (not used in report)

volcano_plot_GOUT_HC <- ggplot(DE_GOUT_vs_HC, aes(log2Fold, -log(p.adj,10))) + # -log10 conversion  
  geom_point(size = 0.3,) +
  xlab(expression("log"[2]*"log2Fold")) + 
  ylab(expression("-log"[10]*"p.adj")) + 
  geom_label(data = SIG_GOUT_vs_HC, 
                    mapping = aes(log2Fold, -log(p.adj,10)), label=SIG_GOUT_vs_HC$symbol)


volcano_plot_GOUT_HC

                           ###//// ANALYSIS OF SIG GENES IN SA vs HC \\\\###

# Melted SIG ET GOUT table for dot plots

melted_SIG_ET_SA_v_HC <- melt(SIG_ET_SA_v_HC,id.vars='SAMPLE_GROUP', measure.vars=c(4:19))
melted_SIG_ET_SA_v_HC <- melted_SIG_ET_SA_v_HC%>% 
  mutate(SAMPLE_GROUP = recode(SAMPLE_GROUP, '0' = 'HC', '1' ='SEPSIS'))

# Box plot exploration of data

SA_boxplot <- ggplot(melted_SIG_ET_SA_v_HC, aes(x= as.factor(SAMPLE_GROUP), y=value)) +
  facet_wrap(~variable, scales = "free") +
  stat_compare_means(method = "t.test", label.x.npc = 'left', label.y.npc = 'top' , label = "p.format", vjust = 1.5) +
  labs(x = NULL, y = "Gene Expression",
       title = "Significant differential gene expression  and  neutrophil count in Septic Arthritis patients",
       subtitle = "P values derived from unpaired T-test") +
  theme_pubclean() +
  labs_pubr() +
  geom_boxplot(width = 0.3)

SA_boxplot

SA_dotplot <- ggplot(melted_scaled_sa,aes(y= SAMPLE_GROUP, x = value)) + 
  geom_point() +
  stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE) +
  facet_wrap(~variable, scales = 'free')
SA_dotplot


# Binomial GLM
binomial_SA_model <- logistf(SAMPLE_GROUP ~ NOD2,
                           data = SIG_ET_SA_v_HC,  family = binomial)
flac(binomial_SA_model, data =SIG_ET_SA_v_HC, model = TRUE)
#FAR2 + TTC39A + PECR + EHF + GALNT6 + S100A8 + TPBG + MPZL2 + TMEM45B + PDZK1IP1 + NOD2 + HPSE + GCNT4 + FAM110C + C10orf99

# This model can't be fitted due to low sample size causing a Hauck-Donner phenomenon 
# due to low sample size, this module logistf was used to over come the effect but couldn't 
# understant the new output and decided to make box plots with T test to calculate significance

plot(profile(binomial_SA_model,variable="NOD2"))
summary(binomial_SA_model)
anova(binomial_SA_model)
hist(rstandard(binomial_SA_model))
autoplot(binomial_SA_model)

# Volcano plots (not used in report)

volcano_plot_SA_HC <- ggplot(DE_SA_vs_HC, aes(log2Fold, -log(p.adj,10))) + # -log10 conversion  
  geom_point(size = 0.3,) +
  xlab(expression("log"[2]*"log2Fold")) + 
  ylab(expression("-log"[10]*"p.adj")) + 
  geom_label(data = SIG_SA_vs_HC, 
             mapping = aes(log2Fold, -log(p.adj,10)), label=SIG_SA_vs_HC$symbol)


volcano_plot_SA_HC




                      ###//// ANALYSIS OF SIG GENES IN SA vs GOUT \\\\###

# Melted SIG ET GOUT table for dot plots

melted_SIG_ET_SA_v_GOUT <- melt(SIG_ET_SA_v_GOUT,id.vars='SAMPLE_GROUP', measure.vars=c(4:13))

# Box plot exploration of data

SA_v_GOUT_boxplot <- ggplot(melted_SIG_ET_SA_v_GOUT, aes(x= as.factor(SAMPLE_GROUP), y=value)) +
  facet_wrap(~variable, scales = "free") +
  stat_compare_means(method = "t.test", label.x.npc = 'left', label.y.npc = 'top' , label = "p.format", vjust = 1.5) +
  labs(x = NULL, y = "Gene Expression",
       title = "Significant differential gene expression  and  neutrophil count between Septic Arthritis and Gout patients",
       subtitle = "P values derived from unpaired T-test") +
  theme_pubclean() +
  labs_pubr() +
  geom_boxplot(width = 0.3)

SA_v_GOUT_boxplot

# No models were tried for this, unlike the others due to low sample size that makes the comparison of 
# the data look like near complete separation and cannot approximate for 0 or 1



                      ###//// PRINCIPAL COMPONENT ANALYSIS \\\\###

PCA = prcomp(t(ET_annotated[,3:29]))
pca_coordinates = data.frame(PCA$x)

# PCA plot
PCA_plot = ggplot(pca_coordinates, aes(x=PC1, y= PC2,
                                       colour = Sample_Information$SAMPLE_GROUP,
                                       size = 3)) +
  labs(title = "Principal Component Analysis of the Expression table",
       subtitle = "PC1 and PC2 account for a majority of the tables' variance") +
  scale_size_continuous(guide = "none") +
  guides(color=guide_legend("Sample Groups")) +
  theme_pubclean() +
  geom_point()

PCA_plot 

ggpar(p = PCA_plot,
      main = "Principal Component Analysis of the Expression table",
      submain = "PC1 and PC2 account for a majority of the tables' variance",
      legend.title = "Sample Groups",
      ggtheme = theme_pubclean())



                        ###//// Appendix tables \\\\###

appendix1 <- sum_stats_gout%>%              
  gt()%>% 
  # Formatting table
  tab_stubhead(label = "Gene Symbols") %>% 
  cols_align(align = "center", columns = everything()) %>% 
  tab_header(title = "Summary Statistics of Significant Genes",
             subtitle = "Sample group: Gout") %>% 
  tab_spanner(label = "Gene Expression Statistics",
              columns = c(2:7)) %>% 
  tab_spanner(label = "Differential Expression Statistics",
              columns = c(8:10)) %>% 
  tab_style(style = cell_text(align = "left"),
            locations = cells_stub(rows = TRUE))

 
appendix2 <- sum_stats_sa%>%              
  gt()%>% 
  # Formatting table
  tab_stubhead(label = "Gene Symbols") %>% 
  cols_align(align = "center", columns = everything()) %>% 
  tab_header(title = "Summary Statistics of Significant Genes",
             subtitle = "Sample group: Spetic arthritis") %>% 
  tab_spanner(label = "Gene Expression Statistics",
              columns = c(2:7)) %>% 
  tab_spanner(label = "Differential Expression Statistics",
              columns = c(8:10)) %>% 
  tab_style(style = cell_text(align = "left"),
            locations = cells_stub(rows = TRUE))

appendix3 <- sum_stats_sa_v_gout%>%              
  gt()%>% 
  # Formatting table
  tab_stubhead(label = "Gene Symbols") %>% 
  cols_align(align = "center", columns = everything()) %>% 
  tab_header(title = "Summary Statistics of Significant Genes",
             subtitle = "Sample group: Gout") %>% 
  tab_spanner(label = "Gene Expression Statistics",
              columns = c(2:7)) %>% 
  tab_spanner(label = "Differential Expression Statistics",
              columns = c(8:10)) %>% 
  tab_style(style = cell_text(align = "left"),
            locations = cells_stub(rows = TRUE))


# Exporting tables
gtsave(appendix1, filename = "appendix1.pdf",)
gtsave(appendix2, filename = "appendix2.pdf",)
gtsave(appendix3, filename = "appendix3.pdf",)
