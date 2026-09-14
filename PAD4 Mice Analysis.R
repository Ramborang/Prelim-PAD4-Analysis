library("ggplot2")
library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)
library(ggpubr)
library(rstatix)

# PAD4 Analysis including Volition, histology, and CAT data. The spreadsheet has 3 different tabs.
PAD4_Vol <- read_excel("J:Park/Projects/Park Mouse Lab LLC/PAD4 Initial Analysis.xlsx", sheet = 'Volition')
PAD4_CAT <- read_excel("J:Park/Projects/Park Mouse Lab LLC/PAD4 Initial Analysis.xlsx", sheet = 'CAT')
PAD4_Hist <- read_excel("J:Park/Projects/Park Mouse Lab LLC/PAD4 Initial Analysis.xlsx", sheet = 'Histology')

# Begin with the Volition analysis
#PAD4_Vol$Groups <- paste(PAD4_Vol$Trauma, PAD4_Vol$Drug, sep = "-") # Grouped by trauma/no trauma. Not enough in the individual time points.
PAD4_Vol$Groups <- paste(PAD4_Vol$Group, PAD4_Vol$Drug, sep = "-")
table(PAD4_Vol$Groups)
table(PAD4_Vol$Group)
table(PAD4_Vol$Drug)

# Create boxplots of Volition Data per group
ggplot(PAD4_Vol, aes(x = `Groups`, y = `[Nucleosome H3.1]-mean (ng/mL)`, fill = `Groups`))+
  geom_boxplot(alpha = 0.7)+
  labs(title = "PAD4 Mice Volition Data by Group")+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.format', hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

kruskal.test(PAD4_Vol$`[Nucleosome H3.1]-mean (ng/mL)`~PAD4_Vol$Groups)

################################################################################
# CAT Analysis on PAD4 Mice (may have to groups 3 hour and 72 hour trauma again for PAD4)
# 7x 3 hour samples and 9x 72 hour samples in the PAD4 sample
PAD4_CAT$Groups <- paste(PAD4_CAT$Group, PAD4_CAT$Drug, sep = '-')
table(PAD4_CAT$Groups)
table(PAD4_CAT$Drug)

# subset the CAT data
PAD4_CAT <- subset(PAD4_CAT, select = c(Groups, `Mean LT`, `Mean PH`, `Mean ttPeak`, `Mean ETP`))
# Create long data
PAD4_CAT_long <- PAD4_CAT %>%
  pivot_longer(-`Groups`, names_to = "variables", values_to = "value")

# Create boxplots of CAT Data per group
plot(ggplot(data = PAD4_CAT_long, aes(x = variables, y = value, fill = factor(`Groups`))) +
       geom_boxplot()+
       labs(title = "PAD4 Mice CAT Data by Group")+
       theme_classic() +
       geom_pwc(method = "dunn_test", label = 'p.format', hide.ns = TRUE)+
       scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
       facet_wrap(~variables, scale = "free"))

# Kruskal Test
stat.test <- PAD4_CAT_long %>%
  group_by(variables) %>%
  kruskal_test(value ~ `Groups`)
print(stat.test)

################################################################################
# Still need to collect data on CitH3 and MPO.