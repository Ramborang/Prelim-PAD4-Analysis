library("ggplot2")
library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)
library(ggpubr)
library(rstatix)

# We only have data for 3 hour defibrotide (we don't have any defibrotide without trauma)
# Defibrotide Analysis including histology, and CAT data. The spreadsheet has 3 different tabs.
Defib_Vol <- read_excel("J:Park/Projects/Park Mouse Lab LLC/Defibrotide Initial Analysis.xlsx", sheet = 'Volition')
Defib_CAT <- read_excel("J:Park/Projects/Park Mouse Lab LLC/Defibrotide Initial Analysis.xlsx", sheet = 'CAT')
Defib_Hist <- read_excel("J:Park/Projects/Park Mouse Lab LLC/Defibrotide Initial Analysis.xlsx", sheet = 'Histology')
Defib_Lung <- read_excel("J:Park/Projects/Park Mouse Lab LLC/Defibrotide Initial Analysis.xlsx", sheet = 'Lung')

# Begin with the Volition analysis
Defib_Vol$Groups <- paste(Defib_Vol$Group, Defib_Vol$Drug, sep = "-")
Defib_Vol <- Defib_Vol %>%
  add_count(Groups, name = 'Count')
table(Defib_Vol$Groups)
table(Defib_Vol$Group)
table(Defib_Vol$Drug)
Defib_Vol$Groups <- paste(Defib_Vol$Groups, " (n = ", Defib_Vol$Count, ")", sep = "")

Defib_Vol <- subset(Defib_Vol, Groups %in% c('Sham-None (n = 32)', '72h-None (n = 15)', '72h-Defibrotide (n = 10)'))
Defib_Vol <- Defib_Vol %>%
  mutate(Groups = fct_relevel(Groups, c('Sham-None (n = 32)', '72h-None (n = 15)','72h-Defibrotide (n = 10)')))

# Create boxplots of Volition Data per group
ggplot(Defib_Vol, aes(x = `Groups`, y = `[Nucleosome H3.1]-mean (ng/mL)`, fill = `Groups`))+
  geom_boxplot(alpha = 0.7)+
  labs(title = "Defibrotide Mice Volition Data by Group")+
  # Add a custom text box inside the plot coordinates
  annotate(
    "text", 
    x = Inf, y = -Inf, # Places it in the bottom right corner
    hjust = 0, vjust = -6,
    label = "p-value <= 0.05: *\np-value <= 0.01: **\np-value <= 0.001: ***\np-value <= 0.0001: ****",
    fontface = "plain", size = 3.5, color = "black", bg = "white")+
  coord_cartesian(clip = "off")+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

kruskal.test(PAD4_Vol$`[Nucleosome H3.1]-mean (ng/mL)`~PAD4_Vol$Groups)
################################################################################
# CAT analysis for Defibrotide mice
Defib_CAT$Groups <- paste(Defib_CAT$Group, Defib_CAT$Drug, sep = '-')
table(Defib_CAT$Groups)
table(Defib_CAT$Drug)

# subset the CAT data
Defib_CAT <- subset(Defib_CAT, select = c(Groups, `Mean LT`, `Mean PH`, `Mean ttPeak`, `Mean ETP`))
# Create long data
Defib_CAT_long <- Defib_CAT %>%
  pivot_longer(-`Groups`, names_to = "variables", values_to = "value")

# Create boxplots of CAT Data per group
plot(ggplot(data = Defib_CAT_long, aes(x = variables, y = value, fill = factor(`Groups`))) +
       geom_boxplot()+
       labs(title = "Defibrotide Mice CAT Data by Group")+
       theme_classic() +
       geom_pwc(method = "dunn_test", label = 'p.format', hide.ns = TRUE)+
       scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
       facet_wrap(~variables, scale = "free"))

# Kruskal Test
stat.test <- Defib_CAT_long %>%
  group_by(variables) %>%
  kruskal_test(value ~ `Groups`)
print(stat.test)
################################################################################
# Finally we have the histology data
# Boxplot of CitH3 Data for each group.
Defib_Hist$Groups <- paste(Defib_Hist$Group, Defib_Hist$Drug, sep = '-')
table(Defib_Hist$Groups)
table(Defib_Hist$Drug)

ggplot(Defib_Hist, aes(x = Groups, y = `CitH3 Counts`, fill = `Groups`, 
                       levels()))+
  labs(title = " Histology Defibrotide Counts in Spleen per Group in Mice",
       x = "Drug and Trauma Groups",
       y = "CitH3  Counts")+
  geom_boxplot(alpha = 0.7)+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Get Counts of MPO categories per group
# First replace -, +, ++, +++ with 0, 1, 2, 3
Defib_Hist$MPO[Defib_Hist$MPO == "-"] <- 0.0
Defib_Hist$MPO[Defib_Hist$MPO == "+"] <- 1.0
Defib_Hist$MPO[Defib_Hist$MPO == "++"] <- 2.0
Defib_Hist$MPO[Defib_Hist$MPO == "+++"] <- 3.0
Defib_Hist$MPO <- as.integer(Defib_Hist$MPO)
print("MPO counts of spleen samples")
data <- Defib_Hist %>% count(Groups, MPO)
matrix_data <- Defib_Hist %>% 
  count(Groups, MPO) %>% 
  pivot_wider(names_from = MPO, values_from = n, values_fill = 0) %>% 
  tibble::column_to_rownames("Groups") %>% 
  as.matrix()
contingency_table <- as.data.frame(matrix_data)
print(contingency_table)
#fisher.test(contingency_table)

# ggplot(data, aes(x = `Groups`, y = `n`, fill = `MPO`))+
#   geom_bar(position = "dodge", stat = "identity")+
#   labs(title = "MPO Counts in Mice Given Defibrotide",
#        y = "MPO Counts")+
#   geom_text(aes(label = `n`),
#             position = position_dodge(width = 0.9),
#             vjust = -0.25)

ggplot(Defib_Hist, aes(x = Groups, y = `MPO`, fill = `Groups`))+
  labs(title = " Histology Defibrotide Counts in Spleen per Group in Mice",
       x = "Drug and Trauma Groups",
       y = "MPO  Counts")+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  geom_boxplot(alpha = 0.7)+
  theme_classic()

################################################################################
# Histology Lung Data
# Boxplot of CitH3 Data for each group.
Defib_Lung$Groups <- paste(Defib_Lung$Group, Defib_Lung$Drug, sep = '-')
table(Defib_Lung$Groups)
table(Defib_Lung$Drug)

ggplot(Defib_Lung, aes(x = Groups, y = `CitH3 Counts`, fill = `Groups`, 
                       levels()))+
  labs(title = " Histology Defibrotide Counts in Lung per Group in Mice",
       x = "Drug and Trauma Groups",
       y = "CitH3  Counts")+
  geom_boxplot(alpha = 0.7)+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Get Counts of MPO categories per group
# First replace -, +, ++, +++ with 0, 1, 2, 3
Defib_Lung$MPO[Defib_Lung$MPO == "-"] <- 0.0
Defib_Lung$MPO[Defib_Lung$MPO == "+"] <- 1.0
Defib_Lung$MPO[Defib_Lung$MPO == "++"] <- 2.0
Defib_Lung$MPO[Defib_Lung$MPO == "+++"] <- 3.0
Defib_Lung$MPO <- as.integer(Defib_Lung$MPO)
print("MPO counts of spleen samples")
data <- Defib_Lung %>% count(Groups, MPO)
matrix_data <- Defib_Lung %>% 
  count(Groups, MPO) %>% 
  pivot_wider(names_from = MPO, values_from = n, values_fill = 0) %>% 
  tibble::column_to_rownames("Groups") %>% 
  as.matrix()
contingency_table <- as.data.frame(matrix_data)
print(contingency_table)
#fisher.test(contingency_table)

# ggplot(data, aes(x = `Groups`, y = `n`, fill = `MPO`))+
#   geom_bar(position = "dodge", stat = "identity")+
#   labs(title = "MPO Counts in Mice Given Defibrotide",
#        y = "MPO Counts")+
#   geom_text(aes(label = `n`),
#             position = position_dodge(width = 0.9),
#             vjust = -0.25)

ggplot(Defib_Lung, aes(x = Groups, y = `MPO`, fill = `Groups`))+
  labs(title = " Histology Defibrotide Counts in Lung per Group in Mice",
       x = "Drug and Trauma Groups",
       y = "MPO  Counts")+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  geom_boxplot(alpha = 0.7)+
  theme_classic()