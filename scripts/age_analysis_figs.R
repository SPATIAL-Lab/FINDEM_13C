

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

#Prep age and sex data
demographic_data <- iso_pid %>%
  inner_join(                                       #merge iso and ind: only select rows kept
    ind_raw %>%
      select(participant_id, age_range, sex),
    by = "participant_id"
  ) %>%
  mutate(age_range = as.numeric(age_range),        #age_range and sex to numeric
         sex = as.numeric(sex),
         
         # Apply descriptive labels from the database key
         age_group = factor(age_range,                     #assign factors age_range 
                            levels = 1:11,
                            labels = c("18–24",
                                       "25–29",
                                       "30–34",
                                       "35–39",
                                       "40–44",
                                       "45–49",
                                       "50–54",
                                       "55–59",
                                       "60–64",
                                       "65–70",
                                       "70+"),
                            ordered = TRUE),
         sex_label = factor(sex,                              #assign factors sex
                            levels = c(1, 2, 77, 3),
                            labels = c("Female",
                                       "Male",
                                       "Other",
                                       "Prefer not to answer"))
  )

#------------------------------------------------------------------------------------------------------------------------




#AGE 

#prep data
age_data <- demographic_data %>%
  filter(!is.na(mean_d13C), 
         !is.na(age_group)
         )

#summary stats by age group
age_summary <- age_data %>%
  group_by(age_group) %>%
  summarise(
    n_participants = n(),
    mean_d13C_group = mean(mean_d13C, na.rm = TRUE),
    median_d13C = median(mean_d13C, na.rm = TRUE),
    sd_d13C = sd(mean_d13C, na.rm = TRUE),
    se_d13C = sd_d13C / sqrt(n_participants),
    ci_low = if_else(n_participants > 1, 
                     mean_d13C_group - qt(0.975, df = n_participants - 1) * se_d13C,
                     NA_real_),
    ci_high = if_else(n_participants > 1,
                      mean_d13C_group +qt(0.975, df = n_participants - 1) * se_d13C,
                      NA_real_),
    .groups = "drop"
  )


#spearman correlation
#convert ordered age groups to numeric ranks
age_test_data <- age_data %>%
  mutate(age_rank = as.numeric(age_group)
  )

#test whether d13C changes monotonically with age group
age_spearman <- cor.test(age_test_data$mean_d13C,
                         age_test_data$age_rank,
                         method = "spearman",
                         exact = FALSE
)


#kruskal-Wallis test
#test for differences in d13C among age groups
age_kw <- kruskal.test(mean_d13C ~ age_group,
                            data = age_data
)




#create statistical annotation

age_kw_p <- age_kw$p.value

age_kw_label <- if (age_kw_p < 0.001) {"Kruskal-Wallis*','~~italic(p)<0.001"
} else {
  paste0("Kruskal-Wallis*','~~italic(p)==", formatC(age_kw_p,
                                                    format = "f",
                                                    digits = 3))
}

#determine annotation positions
y_min_age <- min(age_data$mean_d13C, na.rm = TRUE
)
y_max_age <- max(age_data$mean_d13C, na.rm = TRUE
)

y_n_age <- y_min_age - 1.3
y_test_age <- y_max_age + 1



# P5 d13C by age range
p5 <- ggplot(age_data, aes(x = age_group,
                          y = mean_d13C)
) +
  geom_violin(aes(fill = age_group,                    #violin distribution
                  color = age_group),
              linewidth = 0.8,
              alpha = 0.30,
              trim = FALSE,
              scale = "width"
  ) +
  ggforce::geom_sina(aes(color = age_group),         #individual participants
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
  geom_errorbar(data = age_summary,            #95% confidence interval around group mean
                aes(x = age_group,
                    ymin = ci_low,
                    ymax = ci_high),
                inherit.aes = FALSE,
                width = 0.07,
                linewidth = 0.7,
                color = "black"
  ) +
  geom_point(data = age_summary,                         #group mean shown as a diamond
             aes(x = age_group,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
  ) +
  geom_text(data = age_summary,                                    #number of participants in each age group
            aes(x = age_group,
                y = y_n_age,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "grey35"
  ) +
  annotate(geom = "text",                                 #kruskal-Wallis result
           x = Inf, 
           y = y_test_age,
           label = age_kw_label,
           parse = TRUE,
           hjust = 1.1,
           size = 4
  ) +
  
  scale_fill_viridis_d(option = "D"             #changed from brewer to viridis due to increase in # of categories
  ) +
  scale_color_viridis_d(option = "D"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))                #make room for annotations
  ) +
  labs(title = "δ¹³C by Age Group",
       subtitle = paste("Violin = density | Box = median and IQR |", "Diamond = mean | Error bar = 95% CI"),
       x = "Age group",
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 10,
                                   color = "grey20"),
        axis.title.y = element_text(margin = margin(r = 8)),
        plot.title = element_text(face = "bold", 
                                  size = 15),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40"),
        plot.margin = margin(
          t = 10,
          r = 15,
          b = 10,
          l = 10)
  )


#save p5
ggsave(filename = file.path("outputs", "d13C_by_age_range.png"),
       plot = p5,
       width = 9,
       height = 6,
       dpi = 300
)

