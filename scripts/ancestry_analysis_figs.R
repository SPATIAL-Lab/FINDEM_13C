

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
#-----------------------------------------------------------------------------------------------------------------------

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


#----------------------------------------------------------------------------------------------------------------------


#ANCESTRY

#assign ancestry labels
anc_labels <- c(
  "1"  = "Am. Indian /\nAK Native",
  "2"  = "Asian /\nAsian Am.",
  "3"  = "Nat. HI/Pacific\nIslander",
  "4"  = "Black /\nAfrican Am.",
  "5"  = "White /\nCaucasian",
  "6"  = "Hispanic /\nLatinx",
  "7"  = "Mixed",
  "77" = "Other"
)


#join isotopic data from iso_pid and ancestry data from ind_num
joined_anc <- iso_pid %>%
  inner_join(
    ind_num %>%
      select(participant_id, ancestry),
    by = "participant_id"
  ) %>%
  filter(!is.na(mean_d13C), !is.na(ancestry), as.character(ancestry) %in% names(anc_labels) #filter and convert numeric ancestry to character
  ) %>%
  mutate(ancestry_label = unname(anc_labels[as.character(ancestry)]   #create new ancestry_label column, unname removes "#"
  ),
  ancestry_label = forcats::fct_reorder(ancestry_label,  #forcats package to reorder groups, arrange by median δ13C
                                        mean_d13C,
                                        .fun = median,
                                        na.rm = TRUE)
  )


#summary stats by ancestry group
summary_anc <- joined_anc %>%
  group_by(ancestry_label) %>%                                          #by ancestry calculate:
  summarise(n_participants = n(),                                       #number of participants
            median_d13C = median(mean_d13C, na.rm = TRUE),              #median
            mean_d13C_group = mean(mean_d13C, na.rm = TRUE),            #mean
            sd_d13C = sd(mean_d13C, na.rm = TRUE),                      #standard deviation
            se_d13C = sd_d13C / sqrt(n_participants),                   #standard error
            ci_low = if_else(n_participants > 1,                        #low confidence interval (lower bound of a 95% ci); if_else(condition, value_if_true, value_if_false)
                             mean_d13C_group - qt(0.975, df = n_participants - 1) * se_d13C,
                             NA_real_),
            ci_high = if_else(n_participants > 1,                       #high confidence interval (upper bound of the 95% ci)
                              mean_d13C_group + qt(0.975, df = n_participants - 1) * se_d13C,
                              NA_real_),
            .groups = "drop"
  )


#Is there any significant difference in δ13C among ancestry groups?
# Kruskal-Wallis test
kw_anc <- kruskal.test(
  mean_d13C ~ ancestry_label,             #mean_d13C as a function of ancestry_label
  data = joined_anc                       #data pulled from joined_anc
)

# Pairwise Wilcoxon post-hoc comparisons; chosen due to group size and possibility of non-normal distribution
pw_anc <- pairwise.wilcox.test(
  x = joined_anc$mean_d13C,                    #compare mean_d13c from joined_anc
  g = joined_anc$ancestry_label,               #grouping variable
  p.adjust.method = "BH",                      #Benjamini-Hochberg correction
  exact = FALSE                                #use normal approximation rather than exact p-values 
)

# Convert pairwise results to a tidy data frame
pw_anc_tidy <- as.data.frame(
  as.table(pw_anc$p.value)) %>%
  rename(group1 = Var1,
         group2 = Var2,
         p_adj = Freq) %>%
  filter(!is.na(p_adj)) %>%
  mutate(significance = case_when(p_adj < 0.001 ~ "***",
                                  p_adj < 0.01  ~ "**",
                                  p_adj < 0.05  ~ "*",
                                  TRUE          ~ "ns")) %>%
  arrange(p_adj)





#median d13C for each ancestry group
anc_medians <- joined_anc %>%
  group_by(ancestry_label) %>%
  summarise(median_d13C = median(mean_d13C,
                                 na.rm = TRUE),
            n_participants = n(),
            .groups = "drop"
  )

#add group medians and direction to pairwise comparisons
pw_anc_interpret <- pw_anc_tidy %>%
  left_join(anc_medians %>%
              rename(group1 = ancestry_label,
                     median_group1 = median_d13C,
                     n_group1 = n_participants),
            by = "group1") %>%
  left_join(anc_medians %>%
              rename(group2 = ancestry_label,
                     median_group2 = median_d13C,
                     n_group2 = n_participants),
            by = "group2") %>%
  
  mutate(median_difference =median_group1 - median_group2,
         direction = case_when(median_difference > 0 ~paste(group1,"has higher δ13C"),
                               median_difference < 0 ~paste(group2,"has higher δ13C"),
                               TRUE ~ "Equal medians")) %>%
  arrange(p_adj)




#create label for Kruskal-Wallis result
kw_p <- kw_anc$p.value

kw_label <- if (kw_p < 0.001) {
  "Kruskal-Wallis*','~~italic(p)<0.001"
} else {
  paste0("Kruskal-Wallis*','~~italic(p)==",
         formatC(
           kw_p,
           format = "f",
           digits = 3))
}


#y-axis positions for annotations
y_min_anc <- min(joined_anc$mean_d13C,
                 na.rm = TRUE
)
y_max_anc <- max(joined_anc$mean_d13C,
                 na.rm = TRUE
)
y_n_anc <- y_min_anc - 2.5
y_test_anc <- y_max_anc + 1.5




#plot 1: ancestry plot
p1 <- ggplot(joined_anc,
             aes(x = ancestry_label,
                 y = mean_d13C)
) +
  geom_violin(aes(fill = ancestry_label,          #violin distribution
                  color = ancestry_label),
              linewidth = 0.8,
              alpha = 0.30,
              trim = FALSE,
              scale = "width"
  ) +
  ggforce::geom_sina(aes(color = ancestry_label),    #individual participant values
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42
  ) +
  geom_boxplot(width = 0.13,                       #median and interquartile range
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.55,
               alpha = 0.75
  ) +
  geom_errorbar(data = summary_anc,                  #95% confidence interval around the mean
                aes(x = ancestry_label,
                    ymin = ci_low,
                    ymax = ci_high),
                inherit.aes = FALSE,
                width = 0.07,
                linewidth = 0.7,
                color = "black"
  ) +
  geom_point(data = summary_anc,                           #mean shown as a diamond
             aes(x = ancestry_label,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black",
             alpha = 0.5
  ) +
  geom_text(data = summary_anc,                                   #sample sizes
            aes(x = ancestry_label,
                y = y_n_anc,
                label = paste0("n = ",n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "grey35"
  ) +
  annotate(geom = "text",                                          #overall Kruskal-Wallis result
           x = Inf,
           y = y_test_anc,
           label = kw_label,
           parse = TRUE,
           hjust = 1.1,
           size = 4
  ) +
  scale_fill_brewer(palette = "Dark2"
  ) +
  scale_color_brewer(palette = "Dark2"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))           #make room for n labels and statistical annotation
  ) +
  labs(title = "δ¹³C by Ancestry",
       subtitle = paste("Violin = density | Box = median and IQR |",
                        "Diamond = mean | Error bar = 95% CI"),
       x = NULL,
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9.5,
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



#save p1
ggsave(filename = file.path("outputs",
                            "vio_d13C_anc.png"),
       plot = p1,
       width = 10,
       height = 6,
       dpi = 300
)


