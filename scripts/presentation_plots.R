

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
#https://scrippsco2.ucsd.edu/data/atmospheric-co2-data/sampling-station-records/south-pole/
#C. D. Keeling, S. C. Piper, R. B. Bacastow, M. Wahlen, T. P. Whorf, M. Heimann, and H. A. Meijer, 
#Exchanges of atmospheric CO2 and 13CO2 with the terrestrial biosphere and oceans from 1978 to 2000. I. 
#Global aspects, SIO Reference Series, No. 01-06, Scripps Institution of Oceanography, San Diego, 88 pages, 2001.

#south pole monthly flask data
sp_raw <- read_csv(sp_monthly_file,
                   skip = 58,
                   show_col_types = FALSE
)


ml_monthly_file <-"data/monthly_flask_c13_mlo.csv"
# year range 1980-2024
# https://scrippsco2.ucsd.edu/data/atmospheric-co2-data/sampling-station-records/mauna-loa-observatory-hawaii/
# C. D. Keeling, S. C. Piper, R. B. Bacastow, M. Wahlen, T. P. Whorf, M. Heimann, and  H. A. Meijer, 
#Exchanges of atmospheric CO2 and 13CO2 with the terrestrial biosphere and  oceans from 1978 to 2000.  
#I. Global aspects, SIO Reference Series, No. 01-06, Scripps  Institution of Oceanography, San Diego, 
#88 pages, 2001.     

#mauna loa monthly flask data
ml_raw <- read_csv(ml_monthly_file,
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

# Note: d13C data from law dome represents annual spline values whereas the d13C data from the south pole and mauna 
# loa represent monthly observations converted to annual mean values. 

#law dome: select required columns. Already annual values.
ld_d13C <- ld_raw %>%
  transmute(year = as.integer(`Year AD...10`),
            d13C_atm = as.numeric(`d13CO2 spline (50 yr, permil)`),
            source = "Law Dome") %>%
  filter(!is.na(year),
         !is.na(d13C_atm),
         year <= 1979                     #1979 because mauna loa data starts at 1980
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



#mauna loa flask: select required columns. Given in monthly values so must calculate annual averages.
ml_annual <- ml_raw %>%
  transmute(year = as.integer(Yr),
            d13C_monthly = as.numeric(`13C filled per-mil`))%>%
  filter(!is.na(year),
         !is.na(d13C_monthly))%>%
  group_by(year)%>%
  summarise(d13C_atm = mean(d13C_monthly, na.rm = TRUE),
            n_months = n(),
            .groups = "drop")%>%
  mutate(source = "Mauna Loa flask")


#Law dome and SP combine to get continuous atmospheric record for FINDEM YOF
#atmos_d13C_merge <- bind_rows(ld_d13C,
 #                             sp_annual) %>%
  #arrange(year)%>%
  #filter(year <= 2024
   #      )

#Law dome and ML combine to get continuous atmospheric record for FINDEM YOF
atmos_d13C_merge <- bind_rows(ld_d13C,
                              ml_annual) %>%
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
           x = 1.5,
           y = y_max_p1 + 0.7,
           label = p1_label,
           parse = FALSE,
           hjust = 0.5,
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
  coord_cartesian(ylim = c(-16, -3)) +
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









#possible color palette for age spectrum
#age_fill <- c("70–74" = "chocolate4",
 #             "65–69" = "chocolate3",
  #            "60–64" = "chocolate2",
   #           "55–59" = "darkorange2",
    #          "50–54" = "darkorange1",
     #        "40–44" = "goldenrod2",
      #        "35–39" = "goldenrod1",
       #       "30–34" = "gold1",
        #      "25–29" = "khaki1",
         #     "18–24" = "lemonchiffon"
#)

#age_color <- c("70–74" = "saddlebrown",
 #              "65–69" = "chocolate4",
  #             "60–64" = "chocolate3",
   #            "55–59" = "darkorange4",
    #           "50–54" = "darkorange3",
     #          "45–49" = "darkorange2",
      #         "40–44" = "darkgoldenrod4",
       #        "35–39" = "darkgoldenrod3",
        #       "30–34" = "darkgoldenrod2",
         #      "25–29" = "goldenrod3",
          #     "18–24" = "goldenrod2"
#)





#p1a <- ggplot(p1_data,
 #            aes(x = age_group_plot,
  #               y = mean_d13C)
#) +
 # geom_violin(aes(fill = age_group_plot,
  #            color = age_group_plot),
   #           trim = FALSE,
    #          scale = "width",
     #         alpha = 0.40
  #) +
  #ggforce::geom_sina(aes(color = age_group_plot),
   #                  size = 1.5,
    #                 alpha = 0.50,
     #                maxwidth = 0.70,
      #               seed = 42
  #) +
  #geom_boxplot(width = 0.13,
   #            outlier.shape = NA,
    #           fill = "white",
     #          color = "grey20",
      #         linewidth = 0.55,
       #        alpha = 0.5
  #) +
  #geom_point(data = p1_summary,
   #          aes(x = age_group_plot,
    #             y = mean_d13C_group,
     #            shape = "Group mean"),
      #       inherit.aes = FALSE,
       #      alpha = 0.3,
        #     shape = 23,
         #    size = 3.2,
          #   stroke = 0.7,
           #  fill = "white",
            # color = "black"
  #) +
  #annotate("text",                  #add statistical values to plot
   #        x = 1.5,
    #       y = y_max_p1 + 0.7,
     #      label = p1_label,
      #     parse = FALSE,
       #    hjust = 0.5,
        #   size = 2.7,
         #  fontface = "bold"
  #) +
  #geom_text(data = p1_summary,                                    #number of participants in each age group
   #         aes(x = age_group_plot,
    #            y = y_n_p1,
     #           label = paste0("n = ", n_participants)),
      #      inherit.aes = FALSE,
       #     size = 3.2,
        #    fontface = "bold",
         #   color = "black"
  #) +
  #labs(title = "\nEnamel δ¹³C by Participant Age\n  ",
   #    x = "\nParticipant Age Range",
    #   y = expression(delta^{13}*C[enamel]~("\u2030"))
  #) +
  #guides(fill = "none",
   #      color = "none"
  #) +
  #coord_cartesian(ylim = c(-16, -3)) +
  #scale_fill_manual(values = age_fill
  #) +
  #scale_color_manual(values = age_color
  #) +
  #present_theme

#p1a
#--------------------------------------------------------------------------------------------------------------------




# FIGURE 2
# ENAMEL d13C VS YEAR OF ENAMEL FORMATION
# WITH ATMOSPHERIC d13C CURVE
#-----------------------------------------------------------------------------------------------------------------------
#2. d13C vs. enamel formation age, w/ d13Catm curve overlaid
#p2a uses YOF, p2b uses 5-year YOF bins for violins

#p2a data
p2a_data <- present_pid %>%
  filter(!is.na(YOF_year),
         !is.na(mean_d13C)
         )

#limit atmospheric curve to YOF range
p2a_atmos <- atmos_d13C_merge %>%
  filter(year >= min(p2a_data$YOF_year, na.rm = TRUE),
         year <= max(p2a_data$YOF_year, na.rm = TRUE)
         )

# Figure 2 colors
p2_colors <- c("Law Dome" = "dodgerblue3",
               "Mauna Loa flask" = "red4",
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


#plot 2 option 1
p2a <- ggplot() +
  geom_line(data = p2a_atmos,            #suess curve d13C
            aes(x = year,
                y = d13C_atm,
                color = source
                ),
            linewidth = 1.2
            ) +
  geom_point(data = p2a_data,                    #enamel d13C
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
                                "Mauna Loa flask"),
                     labels = c("Law Dome atmospheric δ¹³C          ",
                                "Mauna Loa atmospheric δ¹³C")
                     ) +
  labs(title = "\nEnamel and Atmospheric δ¹³C Through Time\n  ",
       x = "\nEstimated year of enamel formation",
       y = expression(delta^{13}*C~("\u2030"))
       ) +
  present_theme


p2a




#p2b data with 5-year bins for violin
p2b_data <- present_pid %>%
  filter(!is.na(YOF_year),
         !is.na(mean_d13C)) %>%
  mutate(YOF_5yr = floor(YOF_year / 5) * 5)

p2b_atmos <- atmos_d13C_merge %>%
  filter(year >= min(p2b_data$YOF_5yr, na.rm = TRUE) - 4,
         year <= max(p2b_data$YOF_5yr, na.rm = TRUE) + 4
         )

#summary stats for 5-year YOF groups
p2b_summary <- p2b_data %>%
  group_by(YOF_5yr) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C, na.rm = TRUE),
            .groups = "drop"
            )




#plot 2 option 2
p2b <- ggplot() +
  geom_line(data = p2b_atmos,             #suess curve d13C
            aes(x = year,
                y = d13C_atm,
                color = source),
            linewidth = 1.2
            ) +
  geom_violin(data = p2b_data,               #enamel d13C by 5-year YOF group
              aes(x = YOF_5yr,
                  y = mean_d13C,
                  group = YOF_5yr),
              fill = "turquoise3",
              color = "turquoise4",
              alpha = 0.40,
              trim = FALSE,
              scale = "width",
              width = 4
              ) +
  ggforce::geom_sina(data = p2b_data,             #individual participants
                     aes(x = YOF_5yr,
                         y = mean_d13C,
                         group = YOF_5yr),
                     color = "turquoise4",
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 1.5,
                     seed = 42
                     ) +
  scale_x_continuous(
    breaks = sort(unique(p2b_data$YOF_5yr)),
    labels = sort(unique(p2b_data$YOF_5yr))
  ) +
  annotate("text",
           x = 1960.5,
           y = -5,
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
                       "Mauna Loa flask"),
                     labels = c("Law Dome atmospheric δ¹³C        ",
                                "Mauna Loa atmospheric δ¹³C")
                     ) +
  labs(title = "\nEnamel and Atmospheric δ¹³C Through Time\n  ",
       x = "\nEstimated year of enamel formation",
       y = expression(delta^{13}*C~("\u2030"))
       ) +
  coord_cartesian(ylim = c(-16, -3)) +
  present_theme

p2b

#save p2
ggsave(filename = file.path("outputs", "p2_present.png"),
       plot = p2b,
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
y_n_p3 <- y_min_p3 - 2.7
  
  
  
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
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
             ) +
  annotate("text",                  #add statistical values to plot
           x = 2.5,
           y = y_max_p3 - 0.7,
           label = p3_label,
           parse = FALSE,
           hjust = 0.5,
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
  coord_cartesian(ylim = c(-16, -3)) +
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



#fill colors
p4_fill <- c("Lower" = "darkolivegreen1",
             "Lower/Middle" = "darkolivegreen2",
             "Middle" = "darkolivegreen3",
             "Middle/Upper" = "darkolivegreen4",
             "Upper" = "darkolivegreen"
             )

#outline and dot color
p4_color <- c("Lower" = "yellow4",
              "Lower/Middle" = "olivedrab",
              "Middle" = "palegreen4",
              "Middle/Upper" = "seagreen4",
              "Upper" = "darkgreen"
              )

#plot 4
p4 <- ggplot(p4_data,
             aes(x = child_ecostatus,
                 y = mean_d13C_corrected)
             ) +
  
  geom_violin(aes(fill = child_ecostatus,
              color = child_ecostatus),
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = child_ecostatus),
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
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
             ) +
  annotate("text",
           x = 1.5,
           y = y_max_p4 + 0.5,
           label = p4_label,
           parse = FALSE,
           hjust = 0.5,
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
  scale_fill_manual(values = p4_fill
  ) +
  scale_color_manual(values = p4_color
  ) +
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
#5. d13Ccorr (Suess-corrected) vs. sex at birth (w/ test)

#p5 data
p5_data <- present_pid %>%
  filter(!is.na(mean_d13C_corrected),
         sex %in% c("Female", "Male")) %>%
  droplevels()


#wilcoxon rank-sum test
p5_wilcox <- wilcox.test(mean_d13C_corrected ~ sex,
                         data = p5_data,
                         exact = FALSE
                         )

p5_wilcox


p5_pval <- p5_wilcox$p.value

p5_pval


#p5 label
p5_label <- paste0("Wilcoxon\np = ",
                   format(p5_pval,
                          scientific = TRUE,
                          digits = 3)
)


#summary stats
p5_summary <- p5_data %>%
  group_by(sex) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C_corrected,
                                   na.rm = TRUE),
            .groups = "drop"
            )


y_max_p5 <- max(p5_data$mean_d13C_corrected,
                na.rm = TRUE
                )
y_min_p5 <- min(p5_data$mean_d13C_corrected, 
                na.rm = TRUE
                )
y_n_p5 <- y_min_p5 - 1.3



#fill colors
p5_fill <- c("Female" = "darkorchid3",
             "Male" = "dodgerblue3"
             )

#outline and dot color
p5_color <- c("Female" = "darkorchid4",
              "Male" = "dodgerblue4"
              )


#plot 5
p5 <- ggplot(p5_data,
             aes(x = sex,
                 y = mean_d13C_corrected)
             ) +
  geom_violin(aes(fill = sex,
                  color = sex),
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = sex),
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42
                     ) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.5,
               alpha = 0.50
               ) +
  geom_point(data = p5_summary,
             aes(x = sex,
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
           x = 0.6,
           y = y_max_p5 + 0.5,
           label = p5_label,
           parse = FALSE,
           hjust = 0,
           size = 3,
           fontface = "bold"
           ) +
  geom_text(data = p5_summary,
            aes(x = sex,
                y = y_n_p5,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
            ) +
  scale_fill_manual(values = p5_fill
                    ) +
  scale_color_manual(values = p5_color
                     ) +
  labs(title = "\nSuess-Corrected Enamel δ¹³C by Sex\n ",
       x = "\nSex at birth",
       y = expression(delta^{13}*C[enamel-corrected]~("\u2030"))
       ) +
  guides(fill = "none",
         color = "none"
         ) +
  present_theme


p5


#save p5
ggsave(filename = file.path("outputs", "p5_present.png"),
       plot = p5,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)





#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#6. d13Ccorr (Suess-corrected) vs. ancestry (w/ test)


#p6 data
p6_data <- present_pid %>%
  filter(!is.na(mean_d13C_corrected),
         !is.na(ancestry),
         ancestry != "No answer") %>%
  droplevels()


#Kruskal-Wallis
p6_kw <- kruskal.test(mean_d13C_corrected ~ ancestry,
                      data = p6_data
                      )

p6_kw


# Pairwise comparisons
# BH correction for multiple comparisons
p6_pairwise <- pairwise.wilcox.test(p6_data$mean_d13C_corrected,
                                    p6_data$ancestry,
                                    p.adjust.method = "BH",
                                    exact = FALSE
                                    )

p6_pairwise


#p6 p-value
p6_pval <- p6_kw$p.value

p6_label <- paste0("Kruskal-Wallis\n(all groups included)\np = ",
                   format(p6_pval,
                          scientific = TRUE,
                          digits = 3)
)


p6_summary <- p6_data %>%
  group_by(ancestry) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C_corrected,
                                   na.rm = TRUE),
            .groups = "drop"
            )



y_min_p6 <- min(p6_data$mean_d13C, 
                na.rm = TRUE
                )
y_max_p6 <- max(p6_data$mean_d13C_corrected,
                na.rm = TRUE
                )
y_n_p6 <- y_min_p6 - 1.3



#plot 6
p6 <- ggplot(p6_data,
             aes(x = ancestry,
                 y = mean_d13C_corrected)
             ) + 
  geom_violin(fill = "orangered3",
              color = "orangered4",
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(color = "orangered4",
                     size = 1.5,
                     alpha = 0.50,
                     maxwidth = 0.70,
                     seed = 42) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.55,
               alpha = 0.50
               ) +
  geom_point(data = p6_summary,
             aes(x = ancestry,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
             )+
  annotate("text",
           x = 1.5,
           y = y_max_p6 - 0.3,
           label = p6_label,
           parse = FALSE,
           hjust = 0.5,
           size = 2.5,
           fontface = "bold"
           ) +
  geom_text(data = p6_summary,                                    #number of participants in each ancestry group
            aes(x = ancestry,
                y = y_n_p6,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
  ) +
  labs(title = "\nSuess-Corrected Enamel δ¹³C by Ancestry\n  ",
       x = "\nAncestry",
       y = expression(delta^{13}*C[enamel-corrected]~("\u2030"))
       ) +
  guides(fill = "none",
         color = "none"
         ) +
  present_theme 


p6



#save p6
ggsave(filename = file.path("outputs", "p6_present.png"),
       plot = p6,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#7. d13Ccorr (Suess-corrected) vs. dietary condition (a=lactose, b=allergy; w/ test)

#p7a data
p7a_data <- present_pid %>%
  filter(!is.na(mean_d13C_corrected),
         !is.na(lactose_tolerance)
         )

#p7a wilcoxon test
p7a_wilcox <- wilcox.test(mean_d13C_corrected ~ lactose_tolerance,
                          data = p7a_data,
                          exact = FALSE
                          )

p7a_wilcox


p7a_pval <- p7a_wilcox$p.value

#p7a label
p7a_label <- paste0("Wilcoxon\np = ",
                    format(p7a_pval,
                           scientific = FALSE,
                           digits = 3)
)


p7a_summary <- p7a_data %>%
  group_by(lactose_tolerance) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C_corrected, 
                                   na.rm = TRUE),
            .groups = "drop"
            )


y_max_p7a <- max(p7a_data$mean_d13C_corrected,
                 na.rm = TRUE
                 )
y_min_p7a <- min(p7a_data$mean_d13C_corrected, 
                na.rm = TRUE
                )
y_n_p7a <- y_min_p7a - 1.3



#fill colors
p7_fill <- c("No" = "chartreuse4",
             "Yes" = "chocolate3"
)

#outline and dot color
p7_color <- c("No" = "darkolivegreen",
              "Yes" = "coral4"
)



#plot p7a
p7a <- ggplot(p7a_data,
              aes(x = lactose_tolerance,
                  y = mean_d13C_corrected)
              ) +
  geom_violin(aes(fill = lactose_tolerance,
                  color = lactose_tolerance),
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = lactose_tolerance),
                     size = 1.5,
                     alpha = 0.55,
                     maxwidth = 0.70,
                     seed = 42
                     ) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.5,
               alpha = 0.50
               ) +
  geom_point(data = p7a_summary,
             aes(x = lactose_tolerance,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
  ) +
  annotate("text",
           x = 1.1,
           y = y_max_p7a + 0.5,
           label = p7a_label,
           parse = FALSE,
           hjust = 0,
           size = 3,
           fontface = "bold"
  ) +
  geom_text(data = p7a_summary,
            aes(x = lactose_tolerance,
                y = y_n_p7a,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
  ) +
  scale_fill_manual(values = p7_fill
  ) +
  scale_color_manual(values = p7_color
  ) +
  labs(x = "\nLactose Intolerant",
       y = expression(delta^{13}*C[enamel-corrected]~("\u2030"))
       ) +
  guides(fill = "none",
         color = "none"
         ) +
  present_theme



p7a
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

p7b_data <- present_pid %>%
  filter(!is.na(mean_d13C_corrected),
         !is.na(food_allergies)
         )


p7b_wilcox <- wilcox.test(mean_d13C_corrected ~ food_allergies,
                          data = p7b_data,
                          exact = FALSE
                          )

p7b_wilcox


p7b_pval <- p7b_wilcox$p.value



p7b_label <- paste0("Wilcoxon\np = ",
                    format(p7b_pval,
                           scientific = FALSE,
                           digits = 3)
)


p7b_summary <- p7b_data %>%
  group_by(food_allergies) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C_corrected, na.rm = TRUE),
            .groups = "drop"
            )


y_max_p7b <- max(p7b_data$mean_d13C_corrected,
                 na.rm = TRUE
                 )
y_min_p7b <- min(p7b_data$mean_d13C_corrected, 
                 na.rm = TRUE
)
y_n_p7b <- y_min_p7b - 1.3

#plot 7b
p7b <- ggplot(p7b_data,
              aes(x = food_allergies,
                  y = mean_d13C_corrected)
              ) + 
  geom_violin(aes(fill = food_allergies,
                  color = food_allergies),
              trim = FALSE,
              scale = "width",
              alpha = 0.40
              ) +
  ggforce::geom_sina(aes(color = food_allergies),
                     size = 1.5,
                     alpha = 0.55,
                     maxwidth = 0.70,
                     seed = 42
                     ) +
  geom_boxplot(width = 0.13,
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.5,
               alpha = 0.80
               ) +
  geom_point(data = p7b_summary,
             aes(x = food_allergies,
                 y = mean_d13C_group),
             inherit.aes = FALSE,
             alpha = 0.3,
             shape = 23,
             size = 3.2,
             stroke = 0.7,
             fill = "white",
             color = "black"
  ) +
  annotate("text",
           x = 1.1,
           y = y_max_p7b + 0.5,
           label = p7b_label,
           parse = FALSE,
           hjust = 0,
           size = 3,
           fontface = "bold"
           ) +
  geom_text(data = p7b_summary,
            aes(x = food_allergies,
                y = y_n_p7b,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.2,
            fontface = "bold",
            color = "black"
  ) +
  scale_fill_manual(values = p7_fill
  ) +
  scale_color_manual(values = p7_color
  ) +
  labs(x = "\nFood Allergies",
       y = NULL
       ) +
  guides(fill = "none",
         color = "none"
         ) +
  present_theme


p7b




# Combine dietary plots

p7 <- p7a + p7b +
  plot_annotation(
    title = "\nSuess-Corrected Enamel δ¹³C by Dietary Condition\n "
  )


p7



#save p7
ggsave(filename = file.path("outputs", "p7_present.png"),
       plot = p7,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#8. d13Ccorr map (by state)

#p8 summary stats
state_d13C_corrected <- present_pid %>%
  select(participant_id,
         mean_d13C_corrected) %>%
  inner_join(res_longest %>%
      select(participant_id,
             state),
      by = "participant_id") %>%
  filter(!is.na(mean_d13C_corrected),
         !is.na(state),
         !state %in% c("Puerto Rico",
                       "Guam",
                       "PR",
                       "GU")) %>%
  group_by(state) %>%
  summarise(state_mean_d13C = mean(mean_d13C_corrected, na.rm = TRUE),
            n_participants = n(),
            .groups = "drop") %>%
  mutate(state_mean_d13C = if_else(n_participants >= 3,   #at least 3 participants per state
                                   state_mean_d13C,
                                   NA_real_)
         )


#plot p8
p8 <- plot_usmap(data = state_d13C_corrected,
                 values = "state_mean_d13C",
                 color = "black",
                 linewidth = 0.3,
                 labels = FALSE
                 ) +
  scale_fill_gradientn(colors = c("wheat2",
                                  "tan",
                                  "goldenrod",
                                  "darkorange3",
                                  "chocolate",
                                  "brown4"),
                       na.value = "grey50",
                       name = "Mean corrected δ¹³C (‰)",
                       guide = guide_colorbar(barwidth = 12,
                                              barheight = 0.8,
                                              title.position = "top",
                                              title.hjust = 0.5)
                       ) +
  labs(title = "\nMean Suess-Corrected Enamel δ¹³C by State of Longest Residence",
       subtitle = "Grey = insufficient data"
       ) +
  theme(plot.title = element_text(face = "bold",
                                  size = 14,
                                  hjust = 0.5),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40",
                                     hjust = 0.5),
        legend.position = "bottom",
        legend.direction = "horizontal"
        )


p8





#save p8
ggsave(filename = file.path("outputs", "p8_present.png"),
       plot = p8,
       width = 10,
       height = 6,
       units = "in",
       dpi = 300
)





























