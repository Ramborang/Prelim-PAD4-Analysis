library("ggplot2")
library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)
library(ggpubr)
library(rstatix)

# Prelim analysis for aged mice. Import data and cut columns
Aged_Data <- read_excel("J:Park/Projects/Park Mouse Lab LLC/Aged Mouse Data 9-11-2026.xlsx")
Aged_Data$Groups <- paste(Aged_Data$Group, Aged_Data$Drug, sep = "-")
Aged_Data <- subset(Aged_Data, select = -c(LT, PH, ttPeak, ETP, `Vel Index`, `Start Tail`))
Aged_Data <- Aged_Data[ , !(names(Aged_Data) %in% cols_to_drop)]

# Start with Group Count
table(Aged_Data$Groups)
table(Aged_Data$Drug)
table(Aged_Data$Group)

################################################################################
# Boxplot of CitH3 Data for each group.
ggplot(Aged_Data, aes(x = Groups, y = `CitH3 Counts`, fill = `Groups`))+
  labs(title = " Histology CitH3 Counts per Group in Aged Mice",
       x = "Drug and Trauma Groups",
       y = "CitH3  Counts")+
  geom_boxplot(alpha = 0.7)+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.format', hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Get Counts of MPO categories per group
print("MPO counts of spleen samples")
data <- Aged_Data %>% count(Groups, MPO)
matrix_data <- Aged_Data %>% 
  count(Groups, MPO) %>% 
  pivot_wider(names_from = MPO, values_from = n, values_fill = 0) %>% 
  tibble::column_to_rownames("Groups") %>% 
  as.matrix()
contingency_table <- as.data.frame(matrix_data)
print(contingency_table)
fisher.test(contingency_table)

ggplot(data, aes(x = `Groups`, y = `n`, fill = `MPO`))+
  geom_bar(position = "dodge", stat = "identity")+
  labs(title = "MPO Counts in Aged Mice",
       y = "MPO Counts")+
  geom_text(aes(label = `n`),
            position = position_dodge(width = 0.9),
            vjust = -0.25)

################################################################################
# subset the CAT data
Aged_CAT <- subset(Aged_Data, select = c(Groups, `Mean LT`, `Mean PH`, `Mean ttPeak`, `Mean ETP`))
# Create long data
Aged_CAT_long <- Aged_CAT %>%
  pivot_longer(-`Groups`, names_to = "variables", values_to = "value")

# Create boxplots of CAT Data per group
plot(ggplot(data = Aged_CAT_long, aes(x = variables, y = value, fill = factor(`Groups`))) +
       geom_boxplot()+
       labs(title = "Aged Mice CAT Data by Group")+
       theme_classic() +
       geom_pwc(method = "dunn_test", label = 'p.format', hide.ns = TRUE)+
       scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
       facet_wrap(~variables, scale = "free"))

# Kruskal Test
stat.test <- Aged_CAT_long %>%
  group_by(variables) %>%
  kruskal_test(value ~ `Groups`)
print(stat.test)

################################################################################
# Subset Volition data
Aged_Vol <- subset(Aged_Data, select = c(Groups,`[Nucleosome H3.1]-mean (ng/mL)`))
#pivot long
Aged_Vol_long <- Aged_Vol %>%
  pivot_longer(-`Groups`, names_to = "variables", values_to = "value")

# Boxplot of nucleosome data from Volition
ggplot(Aged_Data, aes(x = `Groups`, y = `[Nucleosome H3.1]-mean (ng/mL)`, fill = `Groups`))+
  geom_boxplot(alpha = 0.7)+
  labs(title = "Aged Mice Volition Data by Group")+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.format', hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Kruskal Wallis test
stat.test <- Aged_Vol_long %>%
  group_by(variables) %>%
  kruskal_test(value ~ `Groups`)
stat.test
