

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
#-------------------------------------------------------------------------------------------------------------------------


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
#-------------------------------------------------------------------------------------------------------------------------



#SEX

#prep sex data
sex_data <- demographic_data %>%
  filter(
    !is.na(mean_d13C),
    !is.na(sex_label),
    sex_label != "Prefer not to answer"
  ) %>%
  droplevels()


#summary stats by sex
sex_summary <- sex_data %>%
  group_by(sex_label) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C, na.rm = TRUE),
            median_d13C = median(mean_d13C, na.rm = TRUE),
            sd_d13C = sd(mean_d13C, na.rm = TRUE),
            se_d13C = sd_d13C / sqrt(n_participants),
            ci_lower = mean_d13C_group - qt(0.975, df = n_participants - 1) * se_d13C,
            ci_upper = mean_d13C_group + qt(0.975, df = n_participants - 1) * se_d13C,
            .groups = "drop"
  )

#female vs Male Wilcoxon test
sex_two_group <- sex_data %>%
  filter(sex_label %in% c("Female", "Male")) %>%
  droplevels()


sex_wilcox <- wilcox.test(mean_d13C ~ sex_label,
                          data = sex_two_group,
                          exact = FALSE,
                          conf.int = TRUE
)
#create Wilcoxon label for figure
sex_p <- sex_wilcox$p.value

sex_test_label <- if (sex_p < 0.001) {"Wilcoxon*','~~italic(p)<0.001"
} else {
  paste0("Wilcoxon*','~~italic(p)==", formatC(sex_p,
                                              format = "f",
                                              digits = 3))
}





#determine annotation positions
y_min_sex <- min(sex_data$mean_d13C, na.rm = TRUE
)
y_max_sex <- max(sex_data$mean_d13C, na.rm = TRUE
)

y_n_sex <- y_min_sex - 0.35
y_test_sex <- y_max_sex + 0.35




#P6 d13C by sex
p6 <- ggplot(sex_data,
             aes(x = sex_label,
                 y = mean_d13C)
) +
  geom_violin(aes(fill = sex_label),       #violin distribution
              trim = FALSE,
              scale = "width",
              color = NA,
              alpha = 0.35
  ) +
  ggforce::geom_sina(aes(color = sex_label),              #individual participants
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.65,
                     seed = 42
  ) +
  geom_boxplot(width = 0.14,                  #median and interquartile range
               outlier.shape = NA,
               fill = "white",
               color = "gray20",
               linewidth = 0.45,
               alpha = 0.75
  ) +
  geom_point(data = sex_summary,                       #group mean shown as a diamond
             aes(x = sex_label,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             shape = 18,
             size = 3.5,
             color = "black"
  ) +
  geom_errorbar(data = sex_summary,              #95% confidence interval around group mean
                aes(x = sex_label,
                    ymin = ci_lower,
                    ymax = ci_upper),
                inherit.aes = FALSE,
                width = 0.08,
                linewidth = 0.6,
                color = "black"
  ) +
  geom_text(data = sex_summary,               #number of participants in each group
            aes(x = sex_label,
                y = y_n_sex,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.4,
            fontface = "bold",
            color = "grey35"
  ) +
  annotate(geom = "text",                               #wilcoxon test result
           x = 1.5,
           y = y_test_sex,
           label = sex_test_label,
           parse = TRUE,
           size = 4
  ) +
  scale_fill_brewer(palette = "Set2"
  ) +
  scale_color_brewer(palette = "Set2"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.12))
  ) +
  labs(title = "δ¹³C by Biological Sex at Birth",
       subtitle = paste("Diamond = mean | Error bar = 95% CI |", "Box = median and IQR | Points = participants"),
       x = "Biological sex at birth",
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10, color = "grey40")
  )


#save p6
ggsave(filename = file.path("outputs", "d13C_by_sex.png"),
       plot = p6,
       width = 7,
       height = 6,
       dpi = 300
)


