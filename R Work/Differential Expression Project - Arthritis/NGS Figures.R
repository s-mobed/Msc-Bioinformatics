library(tidyverse)
library(ggplot2)
library(ggfortify)
library(readr)
library(ggpubr)
library(gt)

# Import data 
NGS_table <- read_csv("~/Desktop/NGS table.csv")
View(NGS_table)



table <- NGS_table%>%              
  gt()%>% 
  # Formatting table
  fmt_number(decimals=0, columns = c(2:3),rows = c(3:10)) %>% 
  fmt_number(decimals=3, columns = c(2:3),rows = c(1:2)) %>% 
  cols_align(align = "center", columns = everything()) %>% 
  tab_header(title = "Quast assessment metrics") %>% 
  tab_style(style = cell_text(align = "left"),
            locations = cells_stub(rows = TRUE))

gtsave(table, filename = "NGS table.pdf",path="~/Desktop/")

# Import phred
Phred_scores <- read_csv("~/Desktop/Phred scores.csv")

ggplot(Phred_scores, aes(x=`Phred score`)) +
  geom_line(aes(y=Trimmed),colour='#b2182b',size=1) +
  geom_line(aes(y=Raw),colour="#2166ac",size=1) +
  theme_pubclean() +
  labs(x="Phred scores",
       y="Read count",
       title="Per Sequence quality scores")

NGx_data <- read_csv("~/Desktop/NGx data.csv", 
                     col_types = cols(...4 = col_skip(), ...5 = col_skip(), 
                                      ...6 = col_skip(), ...7 = col_skip(), 
                                      ...8 = col_skip()))

NGx_data[3] <- NGx_data[3]/1000

ggplot(NGx_data, aes(x=Percentage, y=Size, group=Assembly, colour=Assembly)) +
  geom_line() +
  xlim(0,100) +
  theme_pubclean() +
  labs(x="Percentage reference genome covered",
       y="Size (Kbp)",
       title="Combined NGx assembly graphs",
       subtitle = "The effect of long read data on de-novo assembly")

