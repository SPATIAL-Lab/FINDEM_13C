

#SETUP

#Packages

library(readxl)
library(tidyverse)
library(corrplot)
library(ggrepel)
library(usmap)       
library(zipcodeR)    
library(patchwork)   
library(scales)
library(sf)
library(ggforce)
library(GGally)




#Read data

path <- "data/comp.xlsx"

iso_raw <- read_excel(path, sheet = "iso")
ind_raw <- read_excel(path, sheet = "ind")
res_raw <- read_excel(path, sheet = "res")
#------------------------------------------------------------------------------------------------------------------------

#PREP


#prep iso
iso <- iso_raw %>%
  mutate(
    d18O = as.numeric(d18O),                           #ensure d180, d13c,d180.sd, and d13c.sd numeric
    d13C = as.numeric(d13C),
    d18O.sd = as.numeric(na_if(d18O.sd, "NA")),
    d13C.sd = as.numeric(na_if(d13C.sd, "NA")),                               
    n = as.integer(n),                                            #ensure n reads as integer
    arch = factor(arch),
    cohort = case_when(str_starts(participant_id, "LAR") ~ "LAR",        #assign cohorts
                       str_starts(participant_id, "CMU") ~ "CMU",
                       str_starts(participant_id, "TSU") ~ "TSU",
                       str_starts(participant_id, "FINDEM") ~ "FINDEM",
                       str_starts(participant_id, "UTK") ~ "UTK",
                       TRUE ~ "Other")
  )




#iso by participant id, ensure one per participant
iso_pid <- iso %>%
  group_by(participant_id) %>%
  summarise(
    mean_d13C = mean(d13C, na.rm = TRUE),       #calculate mean (d13C, d18O)
    mean_d18O = mean(d18O, na.rm = TRUE),
    n_samples = n(),                            #number of samples for each participant
    cohort = first(cohort),                     #create cohort
    .groups = "drop"                            #remove grouping structure
  )

#prep ind 
ind <- ind_raw %>%
  mutate(ancestry = factor(ancestry,                                #assign factors
                           levels = c(1,2,3,4,5,6,7,77,9),
                           labels = c("Am. Indian/AK Nat.",
                                      "Asian/Asian Am.",
                                      "Nat. HI/PI", 
                                      "Black/African Am.",
                                      "White/Caucasian",
                                      "Hispanic/Latinx",
                                      "Mixed",
                                      "Other",
                                      "No answer")),
         child_ecostatus = factor(child_ecostatus, 
                                  levels = 1:5,
                                  labels = c("Lower",
                                             "Lower/Middle",
                                             "Middle",
                                             "Middle/Upper",
                                             "Upper")),
         across(c(vegan, 
                  vegetarian, 
                  pescatarian, 
                  local_foods, 
                  trad_foods, 
                  lactose_tolerance, 
                  food_allergies),
                ~ factor(.x, 
                         levels = c(0,1), 
                         labels = c("No","Yes")))
  )



#ind numeric version for statistical operations 
ind_num <- ind_raw %>%
  mutate(across(c(ancestry,                  
                  child_ecostatus, 
                  vegan, 
                  vegetarian, 
                  pescatarian,
                  local_foods, 
                  trad_foods, 
                  lactose_tolerance, 
                  food_allergies),
                ~ as.numeric(na_if(.x, "NA")))
  )

#------------------------------------------------------------------------------------------------------------------------



#CHILDHOOD ECONOMIC STATUS

#define labels
eco_labels <- c("1" = "Lower",
                "2" = "Lower/\nMiddle",
                "3" = "Middle",
                "4" = "Middle/\nUpper",
                "5" = "Upper"
)

#prep data
joined_eco <- iso_pid %>%
  inner_join(ind_num %>%
               select(participant_id,
                      child_ecostatus),
             by = "participant_id"
  ) %>%
  filter(!is.na(mean_d13C), 
         !is.na(child_ecostatus),
         child_ecostatus %in% 1:5
  ) %>%
  mutate(eco_label = unname(eco_labels[as.character(child_ecostatus)]),
         eco_label = factor(eco_label, levels = unname(eco_labels), ordered = TRUE)
  )


#summary stats by childhood eco status
summary_eco <- joined_eco %>%
  group_by(eco_label
  ) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C, na.rm = TRUE),
            median_d13C = median(mean_d13C, na.rm = TRUE),
            sd_d13C = sd(mean_d13C, na.rm = TRUE),
            se_d13C = sd_d13C /sqrt(n_participants),
            ci_low = if_else(n_participants > 1, mean_d13C_group - qt(0.975, df = n_participants - 1) * se_d13C,
                             NA_real_),
            ci_high = if_else(n_participants > 1, mean_d13C_group + qt(0.975, df = n_participants - 1) * se_d13C,
                              NA_real_),
            .groups = "drop"
  )


#kruskal-Wallis test
kw_eco <- kruskal.test(mean_d13C ~ eco_label,
                       data = joined_eco
)

#pairwise Wilcoxon post-hoc tests
pw_eco <- pairwise.wilcox.test(x = joined_eco$mean_d13C,
                                     g = joined_eco$eco_label,
                                     p.adjust.method = "BH",
                                     exact = FALSE
)

#convert pairwise p-values to tidy data frame
pw_eco_tidy <- as.data.frame(as.table(pw_eco$p.value)
) %>%
  rename(group1 = Var1,
         group2 = Var2,
         p_adj = Freq
  ) %>%
  filter(!is.na(p_adj)
  ) %>%
  mutate(significance = case_when(p_adj < 0.001 ~ "***",
                                  p_adj < 0.01  ~ "**",
                                  p_adj < 0.05  ~ "*",
                                  TRUE          ~ "ns")
  ) %>%
  arrange(p_adj
  )



#create Kruskal-Wallis label for figure
kw_eco_p <- kw_eco$p.value
kw_eco_label <- if (kw_eco_p < 0.001) {"Kruskal-Wallis*','~~italic(p)<0.001"
} else {paste0("Kruskal-Wallis*','~~italic(p)==",
               formatC(kw_eco_p,
                       format = "f",
                       digits = 3))
}


#determine annotation positions
y_min_eco <- min(joined_eco$mean_d13C, na.rm = TRUE
)
y_max_eco <- max(joined_eco$mean_d13C, na.rm = TRUE
)
y_n_eco <- y_min_eco - 1.7
y_test_eco <- y_max_eco + 0.35


#p4 childhood eco status figure
p4 <- ggplot(joined_eco,
             aes(x = eco_label,
                 y = mean_d13C)
) +
  geom_violin(aes(fill = eco_label,                         #violin distribution
                  color = eco_label),
              linewidth = 0.8,
              alpha = 0.30,
              trim = FALSE,
              scale = "width"
  ) +
  ggforce::geom_sina(aes(color = eco_label),              #individual participants
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42
  ) +
  geom_boxplot(width = 0.13,                         #median and interquartile range
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.55,
               alpha = 0.75
  ) +
  geom_errorbar(data = summary_eco,                       #95% confidence interval around group mean
                aes(x = eco_label,
                    ymin = ci_low,
                    ymax = ci_high),
                inherit.aes = FALSE,
                width = 0.07,
                linewidth = 0.7,
                color = "black"
  ) +
  geom_point(data = summary_eco,                       #group mean shown as a diamond
             aes(x = eco_label,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
  ) +
  geom_text(data = summary_eco,                          #number of participants in each group
            aes(x = eco_label,
                y = y_n_eco,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "grey35"
  ) +
  annotate(geom = "text",                                #overall Kruskal-Wallis result
           x = Inf,
           y = y_test_eco,
           label = kw_eco_label,
           parse = TRUE,
           hjust = 1.1,
           size = 4
  ) +
  scale_fill_brewer(palette = "Dark2"                 #make room for sample-size labels and test result
  ) +
  scale_color_brewer(palette = "Dark2"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))
  ) +
  labs(title = "δ¹³C by Childhood Economic Status",
       subtitle = paste("Violin = density | Box = median and IQR |",
                        "Diamond = mean | Error bar = 95% CI"),
       x = "Childhood economic status",
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,
                                   lineheight = 1.05,
                                   color = "grey20"),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.title = element_text(face = "bold",
                                  size = 15),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40"),
        plot.margin = margin(t = 10,
                             r = 15,
                             b = 10,
                             l = 10)
  )



#save p4
ggsave(filename = file.path("outputs", "d13C_by_ecostat.png"),
       plot = p4,
       width = 9,
       height = 6,
       dpi = 300
)


