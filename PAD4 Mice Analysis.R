library("ggplot2")
library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)
library(ggpubr)
library(rstatix)
theme_set(theme_classic(base_size = 20))

# PAD4 Analysis including Volition, histology, and CAT data. The spreadsheet has 3 different tabs.
PAD4_Vol <- read_excel("J:Park/Projects/Park Mouse Lab LLC/PAD4 Initial Analysis.xlsx", sheet = 'Volition')
PAD4_CAT <- read_excel("J:Park/Projects/Park Mouse Lab LLC/PAD4 Initial Analysis.xlsx", sheet = 'CAT')
PAD4_Hist <- read_excel("J:Park/Projects/Park Mouse Lab LLC/PAD4 Initial Analysis.xlsx", sheet = 'Histology')

# Begin with the Volition analysis
PAD4_Vol$Groups <- paste(PAD4_Vol$Group, PAD4_Vol$Drug, sep = "-") # Grouped by trauma/no trauma. Not enough in the individual time points.
table(PAD4_Vol$Groups)
table(PAD4_Vol$Group)
table(PAD4_Vol$Drug)
PAD4_Vol <- PAD4_Vol %>%
  add_count(Groups, name = 'Count')
PAD4_Vol$Groups <- paste(PAD4_Vol$Groups, " (n = ", PAD4_Vol$Count, ")", sep = "")
print(unique(PAD4_Vol$Groups))

PAD4_Vol <- PAD4_Vol %>%
  mutate(Groups = fct_relevel(Groups, c('Sham-None (n = 22)', 'Sham-PAD4 (n = 10)',
                                        '3h-None (n = 13)', '3h-PAD4 (n = 4)',
                                        '72h-None (n = 13)', '72h-PAD4 (n = 7)')))

# Create boxplots of Volition Data per group
p<-ggplot(PAD4_Vol, aes(x = `Groups`, y = `[Nucleosome H3.1]-mean (ng/mL)`, fill = `Groups`))+
  geom_boxplot(alpha = 0.7, linewidth = 0.5)+
  #labs(title = "PAD4 Mice Volition Data by Group")+
  # Add a custom text box inside the plot coordinates
  annotate(
    "text", 
    x = Inf, y = -Inf, # Places it in the bottom right corner
    hjust = -0.1, vjust = -3,
    label = "p-value <= 0.05: *\np-value <= 0.01: **\np-value <= 0.001: ***\np-value <= 0.0001: ****",
    fontface = "plain", size = 5, color = "black", bg = "white")+
  coord_cartesian(clip = "off")+
  #theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
p+
  scale_x_discrete(labels = c("Sham\n-None\n(n = 22)", "Sham-\nPAD4 Inh\n(n = 10)",
                              '3h-\nNone\n(n = 13)', '3h-\nPAD4 Inh\n(n = 4)',
                              '72h-\nNone\n(n = 13)', '72h-\nPAD4 Inh\n(n = 7)'))+
  scale_fill_discrete(labels = c("Sham-None (n = 22)", "Sham-PAD4 Inhibitor (n = 10)",
                                 '3h-None (n = 13)', '3h-PAD4 Inhibitor (n = 4)',
                                 '72h-None (n = 13)', '72h-PAD4 Inhibitor (n = 7)'))+
  theme(legend.key.spacing.y = unit(0.25, "cm"), legend.justification = c("right", "bottom"))

# p-value cutoffs for *'s
# p > 0.05: ns
# p <= 0.05: *
# p <= 0.01: **
# p <= 0.001: ***
# p <= 0.0001: ****

kruskal.test(PAD4_Vol$`[Nucleosome H3.1]-mean (ng/mL)`~PAD4_Vol$Groups)
pairwise.t.test(PAD4_Vol$`[Nucleosome H3.1]-mean (ng/mL)`,PAD4_Vol$Groups)

################################################################################
p <-ggplot(data = PAD4_CAT_long, aes(x = variables, y = value, fill = factor(`Groups`))) +
      geom_boxplot(linewidth = 0.5, outlier.size = 1)+
      #labs(title = "PAD4 Inhibitor Mice CAT Data by Group")+
      theme_classic(base_size = 24) +
      geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
      facet_wrap(~variables, scale = "free")+
      coord_cartesian(clip = "off")+
      theme_classic(base_size = 24)+
      geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
      
p + scale_fill_discrete(labels = c("Sham-None (n = 38)", "Sham-PAD4\nInhibitor (n = 12)",
                               '3h-None (n = 16)', '3h-PAD4\nInhibitor (n = 7)',
                               '72h-None (n = 17)', '72h-PAD4\nInhibitor (n = 9)'))+
  theme(legend.key.spacing.y = unit(0.5, "cm"))

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
       #labs(title = "PAD4 Inhibitor Mice CAT Data by Group")+
       theme_classic() +
       geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
       scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
       facet_wrap(~variables, scale = "free"))

# Kruskal Test
stat.test <- PAD4_CAT_long %>%
  group_by(variables) %>%
  kruskal_test(value ~ `Groups`)
print(stat.test)

################################################################################
# Still need to collect data on CitH3 and MPO.
PAD4_Hist$Groups <- paste(PAD4_Hist$Group, PAD4_Hist$Drug, sep = "-")
table(PAD4_Hist$Groups)
table(PAD4_Hist$Group)
table(PAD4_Hist$Drug)

ggplot(PAD4_Hist, aes(x = Groups, y = `CitH3 Counts`, fill = `Groups`))+
  labs(title = " Histology CitH3 Counts per Group in PAD4 Inhibitor Arm",
       x = "Drug and Trauma Groups",
       y = "CitH3  Counts")+
  geom_boxplot(alpha = 0.7)+
  theme_classic()+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# Get Counts of MPO categories per group
PAD4_Hist$MPO[PAD4_Hist$MPO == "-"] <- 0
PAD4_Hist$MPO[PAD4_Hist$MPO == "+"] <- 1
PAD4_Hist$MPO[PAD4_Hist$MPO == "++"] <- 2
PAD4_Hist$MPO[PAD4_Hist$MPO == "+++"] <- 3
PAD4_Hist$MPO <- as.integer(PAD4_Hist$MPO)
print("MPO counts of spleen samples")
data <- PAD4_Hist %>% count(Groups, MPO)
matrix_data <- PAD4_Hist %>% 
  count(Groups, MPO) %>% 
  pivot_wider(names_from = MPO, values_from = n, values_fill = 0) %>% 
  tibble::column_to_rownames("Groups") %>% 
  as.matrix()
contingency_table <- as.data.frame(matrix_data)
print(contingency_table)

ggplot(PAD4_Hist, aes(x = Groups, y = `MPO`, fill = `Groups`))+
  labs(title = " Histology MPO Counts per Group in Mice given PAD4 Inhibitor",
       x = "Drug and Trauma Groups",
       y = "MPO  Counts")+
  geom_pwc(method = "dunn_test", label = 'p.signif', hide.ns = TRUE, use.four.stars = TRUE)+
  geom_boxplot(alpha = 0.7)+
  theme_classic()
