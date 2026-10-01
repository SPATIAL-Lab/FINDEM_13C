

#The sequence of figures that I'd like to show in the talk is:

#1. d13C vs. participant age, w/ test for group effect
#2. d13C vs. enamel formation age, w/ d13Catm curve overlaid
#3. d13Ccorr (Suess-corrected) vs. participant age, w/ test for any remaining group effect (presumably there is none?)
#4. d13Ccorr vs. childhood socioeconomic status (w/ test for group effect)
#5. d13Ccorr vs. sex at birth (w/ test)
#6. d13Ccorr vs. ancestry (w/ test)
#7. d13Ccorr vs. dietary condition (lactose, allergy; w/ test)
#8. d13Ccorr map (by state)

#I think this would show most of the cool results and provide a basis for some 
#discussion of what they might mean and how they could be useful. 
#Two things I could use your help with.

#First, we need to make sure all the analyses (after #1) are updated to use the Suess-corrected d13C values.

#Second, for a professional presentation it would be good to adopt a consistent 
#style and format across all the graphs (to the degree possible). 
#I like the violin plots you've been making. We should settle on a background 
#style (either the grey background used in some of your earlier plots or white...I feel like the violins 
#stand out better on the grey, though? The titles and headings you've been adding are good for a 
#research talk (we'll want to remove them later for the paper). For plots #1-3, try to make 
#sure that the data appear in the same order on each plot 
#(i.e., that older samples are always on the left and younger ones on the right).







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

#---------------------------------------------------------------------------------------------------------------------

#Read data

path <- "data/comp.xlsx"

iso_raw <- read_excel(path, sheet = "iso")
ind_raw <- read_excel(path, sheet = "ind")
res_raw <- read_excel(path, sheet = "res")


law_dome_file <- "data/Law_Dome_GHG_2000years.xlsx"
# year range 1020-1996
#Rubino, Mauro; Etheridge, David; Thornton, David; Allison, Colin; Francey, Roger; Langenfelds, Ray; Steele, Paul; 
#Trudinger, Cathy; Spencer, Darren; Curran, Mark; Van Ommen, Tas; & Smith, Andrew (2019): Law Dome Ice Core 2000-Year 
#CO2, CH4, N2O and d13C-CO2. v3. CSIRO. Data Collection. https://doi.org/10.25919/5bfe29ff807fb

#law dome annual spline
ld_raw <- read_excel(law_dome_file,
                     sheet = "Splines fits",
                     skip = 3
)

sp_monthly_file <- "data/monthly_flask_c13_spo.csv"
# year range 1977-2024
#C. D. Keeling, S. C. Piper, R. B. Bacastow, M. Wahlen, T. P. Whorf, M. Heimann, and H. A. Meijer, 
#Exchanges of atmospheric CO2 and 13CO2 with the terrestrial biosphere and oceans from 1978 to 2000. I. 
#Global aspects, SIO Reference Series, No. 01-06, Scripps Institution of Oceanography, San Diego, 88 pages, 2001.

#south pole monthly flask data
sp_raw <- read_csv(sp_monthly_file,
                   skip = 58,
                   show_col_types = FALSE
)
#---------------------------------------------------------------------------------------------------------------------

#PREP DATA

#prep iso
iso <- iso_raw %>%
  mutate(d18O = as.numeric(d18O),                           #ensure d180, d13c,d180.sd, and d13c.sd numeric
         d13C = as.numeric(d13C),
         d18O.sd = as.numeric(na_if(d18O.sd, "NA")),
         d13C.sd = as.numeric(na_if(d13C.sd, "NA")),
         n = as.integer(n),                                            #ensure n reads as integer
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
  summarise(mean_d13C = mean(d13C, na.rm = TRUE),       #calculate mean (d13C, d18O)
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
         sex = factor(sex,
                      levels = c(1, 2, 77, 3),
                      labels = c("Female",
                                 "Male",
                                 "Other",
                                 "Prefer not to answer")),
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


#prep res
res <- res_raw %>%
  mutate(state = str_trim(as.character(state)),              #assign to character, trim whitespace
         zipthree = str_pad(as.character(zipthree),          #zip to character
                            width = 3,                       #ensure pad 3 digits
                            side = "left",                   #pad left side
                            pad = "0"),                      #pad with zero rather than space
         nyears = as.numeric(nyears)                         #nyears as numeric
  ) 



#longest single residence period     
res_longest <- res %>%
  filter(!is.na(state), state != "NA", !is.na(nyears)) %>%        #filter out missing state, "na", and missing nyears
  group_by(participant_id) %>%
  slice_max(order_by = nyears, n = 1, with_ties = FALSE) %>%   #choose largest nyears value, arbitrary tie-breaker
  ungroup()








#YOF prep
#Note: the first age range of 18-25 represents a 7-year period whereas the rest of the age ranges represent 5-year

age_yof <- ind_raw %>%                                      #age_YOF will be used to associate sample with suess data
  transmute(participant_id,
            age_range = as.numeric(age_range),                     #convert age range to numeric
            age_min = case_when(age_range == 1  ~ 18,              #age min: lower range
                                age_range == 2  ~ 25,
                                age_range == 3  ~ 30,
                                age_range == 4  ~ 35,
                                age_range == 5  ~ 40,
                                age_range == 6  ~ 45,
                                age_range == 7  ~ 50,
                                age_range == 8  ~ 55,
                                age_range == 9  ~ 60,
                                age_range == 10 ~ 65,
                                age_range == 11 ~ 70,
                                TRUE ~ NA_real_),
            age_max = case_when(age_range == 1  ~ 24,                #age max:upper range
                                age_range == 2  ~ 29,
                                age_range == 3  ~ 34,
                                age_range == 4  ~ 39,
                                age_range == 5  ~ 44,
                                age_range == 6  ~ 49,
                                age_range == 7  ~ 54,
                                age_range == 8  ~ 59,
                                age_range == 9  ~ 64,
                                age_range == 10 ~ 69,   #changed upper to 69 rather than 70 which was listed in the key
                                age_range == 11 ~ 74,          #capped 70+ at 74 for equation
                                TRUE ~ NA_real_)) %>%
  mutate(age_mean = (age_min + age_max) / 2,                           #calculate mean age
         YOB = 2025 - age_mean,                                      # calculate YOB: Year of birth
         YOF = YOB + 10,                                             #calculate YOF: Year of enamel formation
         YOF_year = round(YOF)
  )

#---------------------------------------------------------------------------------------------------------------------


#suess data prep

# Note: d13C data from law dome represents annual spline values whereas the d13C data from the south pole represents
# monthly observations converted to annual mean values. 

#law dome: select required columns. Already annual values.
ld_d13C <- ld_raw %>%
  transmute(year = as.integer(`Year AD...10`),
            d13C_atm = as.numeric(`d13CO2 spline (50 yr, permil)`),
            source = "Law Dome") %>%
  filter(!is.na(year),
         !is.na(d13C_atm),
         year <= 1976                     #1976 because south pole data starts at 1977
  )





#south pole flask: select required columns. Given in monthly values so must calculate annual averages.
sp_annual <- sp_raw %>%
  transmute(year = as.integer(Yr),
            d13C_monthly = as.numeric(`13C filled per-mil`)) %>%  #"13C filled per-mil" provides a value for all months
  filter(!is.na(year),
         !is.na(d13C_monthly)) %>%
  group_by(year) %>%
  summarise(d13C_atm = mean(d13C_monthly, na.rm = TRUE),  #calculate annual mean from monthly values
            n_months = n(),
            .groups = "drop") %>%
  mutate(source = "South Pole flask")





#combine to get continuous atmospheric record for FINDEM YOF
atmos_d13C_merge <- bind_rows(ld_d13C,
                              sp_annual) %>%
  arrange(year)%>%
  filter(year <= 2024
         )




#assign atmospheric d13C value to estimated YOF
YOF_atmos <- age_yof %>%
  left_join(atmos_d13C_merge,
            by = c("YOF_year" = "year")) %>%
  rename(d13C_atm_YOF = d13C_atm
  )




#join d13Catmos at YOF to iso participant enamel data
iso_pid_atmos <- iso_pid %>%
  left_join(YOF_atmos %>%
              select(participant_id,
                     age_range,
                     age_min,
                     age_max,
                     age_mean,
                     YOB,
                     YOF,
                     YOF_year,
                     d13C_atm_YOF),
            by = "participant_id"
            )



#suess correction
iso_pid_suess <- iso_pid_atmos %>%
  mutate(mean_d13C_corrected = mean_d13C - (d13C_atm_YOF + 6.5)
  )




present_pid <- iso_pid_suess %>%
  inner_join(ind %>%
               select(participant_id,
                      sex,
                      ancestry,
                      child_ecostatus,
                      lactose_tolerance,
                      food_allergies),
             by = "participant_id") %>%
  mutate(age_group_plot = factor(age_range,
                                 levels = 11:1,
                                 labels = c("70–74",
                                            "65–69",
                                            "60–64",
                                            "55–59",
                                            "50–54",
                                            "45–49",
                                            "40–44",
                                            "35–39",
                                            "30–34",
                                            "25–29",
                                            "18–24"),
                                 ordered = TRUE)
         )


#create presentation theme for uniformity
present_theme <- theme_dark(base_size = 13) +
  theme(panel.background = element_rect(fill = "grey60",
                                        color = NA),
        plot.background = element_rect(fill = "white",
                                       color = NA),
        panel.grid.major = element_line(color = "gray40",
                                        linewidth = 0.5),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = 14),
        plot.subtitle = element_text(size = 10,
                                     color = "grey30"),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(size = 9),
        legend.background = element_rect(fill = "white",
                                         color = NA),
        legend.key = element_rect(fill = "white",
                                  color = NA)
  )
#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#1. d13C vs. participant age, w/ test for group effect
# UNCORRECTED d13C VS PARTICIPANT AGE

#p1 data
p1_data <- present_pid %>%
  filter(!is.na(mean_d13C),
         !is.na(age_group_plot)
         )


#Kruskal-Wallis test
p1_kw <- kruskal.test(mean_d13C ~ age_group_plot,
                      data = p1_data
                      )
p1_kw

#statistical label
p1_pval <- p1_kw$p.value

p1_label <- paste0("Kruskal-Wallis \np = ",
                   format(p1_pval,
                          scientific = TRUE,
                          digits = 3)
                   )


#summary stats
p1_summary <- p1_data %>%
  group_by(age_group_plot) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C, na.rm = TRUE),
            .groups = "drop"
            )


y_min_p1 <- min(p1_data$mean_d13C, na.rm = TRUE
                )
y_max_p1 <- max(p1_data$mean_d13C, na.rm = TRUE
                )
y_n_p1 <- y_min_p1 - 1.3


# plot1
p1 <- ggplot(p1_data,
             aes(x = age_group_plot,
                 y = mean_d13C)
             ) +
  geom_violin(fill = "turquoise3",
              color = "turquoise4",
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = "Individual participant"),
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42
                     ) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.55,
               alpha = 0.5
               ) +
  geom_point(data = p1_summary,
             aes(x = age_group_plot,
                 y = mean_d13C_group,
                 shape = "Group mean"),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
             ) +
  annotate("text",                  #add statistical values to plot
           x = 2.03,
           y = y_max_p1 + 0.7,
           label = p1_label,
           parse = FALSE,
           hjust = 0,
           size = 2.7,
           fontface = "bold"
           ) +
  geom_text(data = p1_summary,                                    #number of participants in each age group
            aes(x = age_group_plot,
                y = y_n_p1,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
  ) +
  scale_color_manual(name = NULL,
                     values = c("Individual participant" = "turquoise4")
                     )+
  labs(title = "\nEnamel δ¹³C by Participant Age\n  ",
       x = "\nParticipant Age Range",
       y = expression(delta^{13}*C[enamel]~("\u2030"))
       ) +
  guides(fill = "none",
         color = "none"
         ) +
  present_theme

p1

#save p1
ggsave(filename = file.path("outputs", "p1_present.png"),
       plot = p1,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)



#--------------------------------------------------------------------------------------------------------------------




# FIGURE 2
# ENAMEL d13C VS YEAR OF ENAMEL FORMATION
# WITH ATMOSPHERIC d13C CURVE
#-----------------------------------------------------------------------------------------------------------------------
#2. d13C vs. enamel formation age, w/ d13Catm curve overlaid
#p2_1 uses YOF, p2_2 uses 5-year YOF bins for violins

#p2_1 data
p2_1_data <- present_pid %>%
  filter(!is.na(YOF_year),
         !is.na(mean_d13C)
         )

#limit atmospheric curve to YOF range
p2_1_atmos <- atmos_d13C_merge %>%
  filter(year >= min(p2_1_data$YOF_year, na.rm = TRUE),
         year <= max(p2_1_data$YOF_year, na.rm = TRUE)
         )

# Figure 2 colors
p2_colors <- c("Law Dome" = "dodgerblue3",
               "South Pole flask" = "red4",
               "Participant enamel" = "turquoise4"
               )


p2_correlate_data <- present_pid %>%
  filter(!is.na(mean_d13C),
         !is.na(d13C_atm_YOF)
         )


#pearson correlation
p2_correlate <- cor.test(p2_correlate_data$d13C_atm_YOF,
                         p2_correlate_data$mean_d13C,
                         method = "pearson"
                         )
p2_correlate

#linear regression
p2_lreg <- lm(mean_d13C ~ d13C_atm_YOF,
              data = p2_correlate_data
              )

p2_r2 <- summary(p2_lreg)$r.squared

p2_pval <- coef(summary(p2_lreg))[
  "d13C_atm_YOF",
  "Pr(>|t|)"
]

p2_label <- paste0("R² = ", round(p2_r2, 2),
                   "\np = ",
                   format(p2_pval,
                          scientific = TRUE,
                          digits = 3)
                   )

p2_1 <- ggplot() +
  geom_line(data = p2_1_atmos,            #suess curve d13C
            aes(x = year,
                y = d13C_atm,
                color = source
                ),
            linewidth = 1.2
            ) +
  geom_point(data = p2_1_data,                    #enamel d13C
             aes(x = YOF_year,
                 y = mean_d13C,
                 color = "Participant enamel"),
             size = 2,
             alpha = 0.65,
             position = position_jitter(width = 0.35,
                                        height = 0)
             ) +
  annotate("text",
           x = 1963,
           y = -2,
           label = p2_label,
           hjust = 0,
           vjust = 1,
           size = 2.7,
           fontface = "bold"
  ) +
  scale_color_manual(name = NULL,
                     values = p2_colors,
                     breaks = c("Law Dome",
                                "South Pole flask"),
                     labels = c("Law Dome atmospheric δ¹³C          ",
                                "South Pole atmospheric δ¹³C")
                     ) +
  labs(title = "\nEnamel and Atmospheric δ¹³C Through Time\n  ",
       x = "\nEstimated year of enamel formation",
       y = expression(delta^{13}*C~("\u2030"))
       ) +
  present_theme


p2_1


#p2_2 data with 5-year bins for violin
p2_2_data <- present_pid %>%
  filter(!is.na(YOF_year),
         !is.na(mean_d13C)) %>%
  mutate(YOF_5yr = floor(YOF_year / 5) * 5)

p2_2_atmos <- atmos_d13C_merge %>%
  filter(year >= min(p2_2_data$YOF_5yr, na.rm = TRUE) - 4,
         year <= max(p2_2_data$YOF_5yr, na.rm = TRUE) + 4
         )

p2_2 <- ggplot() +
  geom_line(data = p2_2_atmos,             #suess curve d13C
            aes(x = year,
                y = d13C_atm,
                color = source),
            linewidth = 1.2
            ) +
  geom_violin(data = p2_2_data,               #enamel d13C by 5-year YOF group
              aes(x = YOF_5yr,
                  y = mean_d13C,
                  group = YOF_5yr),
              fill = "turquoise4",
              color = "grey30",
              alpha = 0.40,
              trim = FALSE,
              scale = "width",
              width = 4
              ) +
  ggforce::geom_sina(data = p2_2_data,             #individual participants
                     aes(x = YOF_5yr,
                         y = mean_d13C,
                         group = YOF_5yr),
                     color = "turquoise4",
                     size = 1.5,
                     alpha = 0.55,
                     maxwidth = 1.5,
                     seed = 42
                     ) +
  scale_x_continuous(
    breaks = sort(unique(p2_2_data$YOF_5yr)),
    labels = sort(unique(p2_2_data$YOF_5yr))
  ) +
  annotate("text",
           x = 1960.5,
           y = -2,
           label = p2_label,
           hjust = 0,
           vjust = 1,
           size = 2.7,
           fontface = "bold"
  ) +
  scale_color_manual(name = NULL,
                     values = p2_colors,
                     breaks = c(
                       "Law Dome",
                       "South Pole flask"),
                     labels = c("Law Dome atmospheric δ¹³C        ",
                                "South Pole atmospheric δ¹³C")
                     ) +
  labs(title = "\nEnamel and Atmospheric δ¹³C Through Time\n  ",
       x = "\nEstimated year of enamel formation",
       y = expression(delta^{13}*C~("\u2030"))
       ) +
  present_theme

p2_2

#save p2
ggsave(filename = file.path("outputs", "p2_present.png"),
       plot = p2_2,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#3. d13Ccorr (Suess-corrected) vs. participant age, w/ test for any remaining group effect (presumably there is none?)

#p3_data
p3_data <- present_pid %>%
  filter(!is.na(mean_d13C_corrected),
         !is.na(age_group_plot)
         )


#Kruskal-Wallis
p3_kw <- kruskal.test(mean_d13C_corrected ~ age_group_plot,
                      data = p3_data
                      )

p3_kw



#statistical label
p3_pval <- p3_kw$p.value

p3_label <- paste0("Kruskal-Wallis\np = ",
                   format(p3_pval,
                          scientific = TRUE,
                          digits = 3)
                   )


#summary stats
p3_summary <- p3_data %>%
  group_by(age_group_plot) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C_corrected,
                                   na.rm = TRUE),
            .groups = "drop"
            )

y_min_p3 <- min(p3_data$mean_d13C_corrected, 
                na.rm = TRUE
                )
y_max_p3 <- max(p3_data$mean_d13C_corrected,
                na.rm = TRUE
                )
y_n_p3 <- y_min_p3 - 1.3
  
  
  
#plot 3
p3 <- ggplot(p3_data,
             aes(x = age_group_plot,
                 y = mean_d13C_corrected)
             ) +
  geom_violin(fill = "turquoise3",
              color = "turquoise4",
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = "Individual participant"),
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42
                     ) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.55,
               alpha = 0.5
               ) +
  geom_point(data = p3_summary,
             aes(x = age_group_plot,
                 y = mean_d13C_group,
                 shape = "Group mean"),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
             ) +
  annotate("text",                  #add statistical values to plot
           x = 2.03,
           y = y_max_p3 + 0.7,
           label = p3_label,
           parse = FALSE,
           hjust = 0,
           size = 2.7,
           fontface = "bold"
           ) +
  geom_text(data = p3_summary,                                    #number of participants in each age group
            aes(x = age_group_plot,
                y = y_n_p3,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
            ) +
  scale_color_manual(name = NULL,
                     values = c("Individual participant" = "turquoise4")
                     )+
  labs(title = "\nSuess-Corrected Enamel δ¹³C by Participant Age\n ",
       x = "\nParticipant Age Range",
       y = expression(delta^{13}*C[enamel-corrected]~("\u2030"))
       ) +
  
  guides(fill = "none",
         color = "none"
         ) +
  
  present_theme


p3

#save p3
ggsave(filename = file.path("outputs", "p3_present.png"),
       plot = p3,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#4. d13Ccorr (Suess-corrected) vs. childhood socioeconomic status (w/ test for group effect)

#p4 data
p4_data <- present_pid %>%
  filter(!is.na(mean_d13C_corrected),
         !is.na(child_ecostatus)
         )


#Kruskal-Wallis
p4_kw <- kruskal.test(mean_d13C_corrected ~ child_ecostatus,
                      data = p4_data
                      )

p4_kw


p4_pval <- p4_kw$p.value

p4_label <- paste0("Kruskal-Wallis\n(all groups included)\np = ",
  format(p4_pval,
         scientific = FALSE,
         digits = 3)
  )

#summary stats
p4_summary <- p4_data %>%
  group_by(child_ecostatus) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C_corrected,
                                   na.rm = TRUE),
            .groups = "drop"
            )

y_min_p4 <- min(p4_data$mean_d13C_corrected, 
                na.rm = TRUE
                )
y_max_p4 <- max(p4_data$mean_d13C_corrected,
                na.rm = TRUE
                )
y_n_p4 <- y_min_p4 - 1.3



#plot 4
p4 <- ggplot(p4_data,
             aes(x = child_ecostatus,
                 y = mean_d13C_corrected)
             ) +
  
  geom_violin(fill = "turquoise3",
              color = "turquoise4",
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = "Individual participant"),
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42
                     ) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.55,
               alpha = 0.50
               ) +
  geom_point(data = p4_summary,
             aes(x = child_ecostatus,
                 y = mean_d13C_group,
                 shape = "Group mean"),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
             ) +
  annotate("text",
           x = 1.2,
           y = y_max_p4 + 0.5,
           label = p4_label,
           parse = FALSE,
           hjust = 0,
           size = 2.7,
           fontface = "bold"
           ) +
  geom_text(data = p4_summary,                         #number of participants in each socioeconomic-status group
            aes(x = child_ecostatus,
                y = y_n_p4,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
  ) +
  scale_color_manual(name = NULL,
                     values = c("Individual participant" = "turquoise4")
  )+
  labs(title = "\nSuess-Corrected Enamel δ¹³C by Childhood Socioeconomic Status\n ",
       x = "\nChildhood socioeconomic status",
       y = expression(delta^{13}*C[enamel-corrected]~("\u2030"))
       ) +
  guides(fill = "none",
         color = "none"
         ) +
  present_theme


p4

#save p4
ggsave(filename = file.path("outputs", "p4_present.png"),
       plot = p4,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)





#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------













#---------------------------------------------------------------------------------------------------------------------


























