

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

#------------------------------------------------------------------------------------------------------------------------




#DIETARY PRACTICES

#prep plot3
#define dietary variables and labels
diet_vars <- c("vegan", 
               "vegetarian", 
               "pescatarian", 
               "local_foods",
               "trad_foods", 
               "lactose_tolerance", 
               "food_allergies")

diet_labels <- c("vegan" = "Vegan",
                 "vegetarian" = "Vegetarian",
                 "pescatarian" = "Pescatarian  ",
                 "local_foods" = "Local \nFoods",
                 "trad_foods" = "Traditional \nFoods",
                 "lactose_tolerance" = "Lactose \nIntolerant",
                 "food_allergies" = "Food \nAllergies")




#prep dietary data
#wide format for summary table, one row per participant
joined_diet_wide <- iso_pid %>%
  inner_join(ind_num %>%
               select(participant_id,
                      all_of(diet_vars)),
             by = "participant_id") %>%
  filter(!is.na(mean_d13C)
  )

# long format for violin plot, one row per participant × dietary variable
joined_diet <- joined_diet_wide %>%
  pivot_longer(cols = all_of(diet_vars),
               names_to = "diet_var",
               values_to = "response") %>%
  filter(!is.na(response), response %in% c(0, 1)) %>%
  mutate(response = factor(response, 
                           levels = c(0, 1),
                           labels = c("No", "Yes")),
         diet_var = unname(diet_labels[diet_var]),
         diet_var = factor(diet_var, levels = unname(diet_labels))
  )


#dietary prevalence
diet_prevalence <- joined_diet_wide %>%
  summarise(across(all_of(diet_vars), ~ mean(.x == 1, na.rm = TRUE) * 100)
  ) %>%
  pivot_longer(cols = everything(),
               names_to = "diet_variable",
               values_to = "percent_yes"
  ) %>%
  mutate(diet_variable = unname(diet_labels[diet_variable])
  )





#by d13C quartile for the proportion plot
joined_quart <- joined_diet_wide %>%
  filter(if_all(all_of(diet_vars), ~ !is.na(.x))
  ) %>%
  mutate(d13C_quartile = ntile(mean_d13C, 4),
         d13C_quartile = factor(d13C_quartile,
                                labels = c("Q1 (lowest)",
                                           "Q2",
                                           "Q3",
                                           "Q4 (highest)"))
  )



#proportion answering Yes in each quartile
summary_quart <- joined_quart %>%
  group_by(d13C_quartile
  ) %>%
  summarise(across(all_of(diet_vars), ~ mean(.x == 1, na.rm = TRUE)),
            .groups = "drop"
  ) %>%
  pivot_longer(cols = -d13C_quartile,
               names_to = "diet_var",
               values_to = "proportion"
  ) %>%
  mutate(diet_var = unname(diet_labels[diet_var]),
         diet_var = factor(diet_var, levels = unname(diet_labels))
  )



#plot 3 dietary-practice quartile figure
p3 <- ggplot(summary_quart,
             aes(x = diet_var,
                 y = proportion,
                 color = d13C_quartile,
                 group = d13C_quartile)
) +
  geom_line(linewidth = 1.2,
            alpha = 0.8
  ) +
  geom_point(size = 3
  ) +
  scale_color_manual(values = c("olivedrab",
                                "skyblue3",
                                "coral",
                                "purple"),
                     name = "δ¹³C quartile"
  ) +
  scale_y_continuous(labels = percent_format(),
                     limits = c(0, 1)
  ) +
  labs(title = "Dietary Practice by δ¹³C Quartile",
       subtitle = "Each line = one participant δ¹³C quartile",
       x = NULL,
       y = "Proportion answering Yes"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40"
        ),
        legend.position = "right"
  )


#save p3
ggsave(filename = file.path("outputs", "diet_d13C_quart.png"),
       plot = p3,
       width = 9,
       height = 6,
       dpi = 300
)



#p3 option 2 violin, box, sina
#d13C distributions by dietary practice
#summary stats
summary_diet <- joined_diet %>%
  group_by(diet_var,
           response
  ) %>%
  summarise(n_participants = n(),
            median_d13C = median(
              mean_d13C,
              na.rm = TRUE),
            .groups = "drop"
  )



#wilcoxon tests
wilcox_results <- joined_diet %>%
  group_by(diet_var
  ) %>%
  filter(all(c("Yes", "No") %in% response)
  ) %>%
  summarise(n_yes = sum(response == "Yes"),
            n_no = sum(response == "No"),
            median_yes = median(mean_d13C[response == "Yes"],na.rm = TRUE),
            median_no = median(mean_d13C[response == "No"],na.rm = TRUE),
            p_value = wilcox.test(mean_d13C[response == "Yes"],
                                  mean_d13C[response == "No"],
                                  exact = FALSE)$p.value,
            .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value,
                          method = "BH")
  )



#labels for adjusted p-values
sig_labels <- wilcox_results %>%
  select(diet_var,
         p_adj
  ) %>%
  mutate(response = factor("Yes",levels = c("No", "Yes")),
         label = case_when(p_adj < 0.001 ~ "p(adj) < 0.001", 
                           TRUE ~ paste0("p(adj) = ", formatC(p_adj,
                                                              format = "f",
                                                              digits = 3)))
  )



#p3_2 violin plot for all dietary variables
diet_y_min <- min(joined_diet$mean_d13C,na.rm = TRUE
)
diet_y_max <- max(joined_diet$mean_d13C,na.rm = TRUE
)

p3_2 <- ggplot(joined_diet,
               aes(x = response,
                   y = mean_d13C)
) +
  geom_violin(aes(fill = response),              #violin distribution
              color = NA,
              alpha = 0.25,
              trim = FALSE,
              scale = "width"
  ) +
  geom_boxplot(width = 0.12,                    #median and IQR
               outlier.shape = NA,
               fill = "white",
               color = "grey20",
               linewidth = 0.5,
               alpha = 0.6
  ) +
  ggforce::geom_sina(aes(color = response),            #individual participants
                     size = 1.2,
                     alpha = 0.45,
                     maxwidth = 0.6,
                     seed = 42
  ) +
  geom_text(data = summary_diet,                        #participant counts
            aes(x = response,
                y = diet_y_min - 0.7,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 2.8,
            color = "grey30"
  ) +
  geom_text(data = sig_labels,aes(x = 1.5,                        #adjusted p-values
                                  y = diet_y_max + 0.2,
                                  label = label),
            inherit.aes = FALSE,
            size = 2.8,
            color = "grey30"
  ) +
  scale_fill_manual(values = c("No" = "skyblue2",
                               "Yes" = "salmon2")
  ) +
  scale_color_manual(values = c("No" = "dodgerblue2",
                                "Yes" = "chocolate2")
  ) +
  facet_wrap(~ diet_var,
             nrow = 2,
             scales = "free_x"
  ) +
  labs(title = "δ¹³C by Dietary Practice",
       subtitle = "Violin = density | Box = median and IQR | Points = individual participants",
       x = NULL,
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 12
  ) +
  theme(legend.position = "none",
        panel.grid.major.x =element_blank(),
        panel.grid.minor =element_blank(),
        strip.text = element_text(face = "bold",
                                  size = 11),
        strip.background = element_rect(fill = "grey95",
                                        color = NA),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40")
  )


#save p3_2
ggsave(filename = file.path("outputs", "vio_d13C_diet.png"),
       plot = p3_2,
       width = 14,
       height = 8,
       dpi = 300
)




#focused, lactose intolerance and food allergies
#prep focused dietary data
focus_diet_labels <- c(lactose_tolerance = "Lactose Intolerance",
                       food_allergies = "Food Allergies"
)
focus_diet_dat <- iso_pid %>%
  inner_join(ind_num %>%
               select(participant_id,
                      lactose_tolerance,
                      food_allergies),
             by = "participant_id"
  ) %>%
  select(participant_id,
         mean_d13C,
         lactose_tolerance,
         food_allergies
  ) %>%
  pivot_longer(cols = c(lactose_tolerance,
                        food_allergies),
               names_to = "diet_variable",
               values_to = "response"
  ) %>%
  filter(!is.na(mean_d13C),
         !is.na(response),
         response %in% c(0, 1)
  ) %>%
  mutate(response = factor(response,
                           levels = c(0, 1),
                           labels = c("No", "Yes")),
         diet_variable = recode(diet_variable, !!!focus_diet_labels),
         diet_variable = factor(diet_variable, levels = unname(focus_diet_labels))
  )



#focused dietary summary stats
focus_diet_summary <- focus_diet_dat %>%
  group_by(diet_variable,
           response
  ) %>%
  summarise(n_participants = n(),
            mean_d13C_group = mean(mean_d13C, na.rm = TRUE),
            median_d13C = median(mean_d13C, na.rm = TRUE),
            sd_d13C = sd(mean_d13C, na.rm = TRUE),
            .groups = "drop"
  )



#wilcoxon tests for focused dietary variables
focus_wilcox <- focus_diet_dat %>%
  group_by(diet_variable
  ) %>%
  summarise(test = list(wilcox.test(mean_d13C ~ response, exact = FALSE)),
            .groups = "drop"
  ) %>%
  mutate(W = map_dbl(test, ~ unname(.x$statistic)),
         p_value = map_dbl(test, ~ .x$p.value),
         p_adjusted = p.adjust(p_value,
                               method = "BH"),
         label = case_when(p_adjusted < 0.001 ~ "italic(p)[adj] < 0.001",
                           TRUE ~ paste0("italic(p)[adj] == ", format(round(p_adjusted,3),
                                                                      nsmall = 3)))
  )



#focused dietary figure
focus_y_min <- min(focus_diet_dat$mean_d13C, na.rm = TRUE
)
focus_y_max <- max(focus_diet_dat$mean_d13C, na.rm = TRUE
)


p_focus_diet <- ggplot(focus_diet_dat, 
                       aes(x = response,
                           y = mean_d13C)
) +
  geom_violin(aes(fill = response),              #distribution
              trim = FALSE,
              scale = "width",
              alpha = 0.35,
              color = NA
  ) +
  ggforce::geom_sina(color = "dodgerblue3",          #individual participants
                     size = 1.4,
                     alpha = 0.45,
                     maxwidth = 0.70,
                     seed = 42
  ) +
  geom_boxplot( width = 0.14,                         #median and IQR
                outlier.shape = NA,
                fill = "white",
                color = "grey20",
                linewidth = 0.5,
                alpha = 0.75
  ) +
  geom_text(data = focus_diet_summary,                   #participant counts
            aes(x = response,
                y = focus_y_min - 0.25,
                label = paste0("n = ", n_participants)),
            inherit.aes = FALSE,
            size = 3.4,
            fontface = "bold",
            color = "grey30"
  ) +
  geom_text(data = focus_wilcox,                    #adjusted p-values
            aes(x = 1.5,
                y = focus_y_max + 0.2,
                label = label),
            inherit.aes = FALSE,
            parse = TRUE,
            size = 4
  ) +
  facet_wrap(~ diet_variable, nrow = 1
  ) +
  scale_fill_manual(values = c("No" = "skyblue2",
                               "Yes" = "salmon2")
  ) +
  labs(title = "δ¹³C by Food Allergies and Lactose Intolerance",
       subtitle = paste("Violin = density | Box = median and IQR |",
                        "Points = individual participants"),
       x = NULL,
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold",
                                  size = 12),
        strip.background = element_rect(fill = "grey95",
                                        color = NA),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40")
  )


#save
ggsave(filename = file.path("outputs", "violin_d13C_lactose_allergies.png"),
       plot = p_focus_diet,
       width = 9,
       height = 6,
       dpi = 300
)



#dietary summary table
diet_group_summary <- joined_diet_wide %>%
  pivot_longer(cols = all_of(diet_vars),
               names_to = "diet_var",
               values_to = "response"
  ) %>%
  filter(!is.na(response),
         response %in% c(0, 1)
  ) %>%
  mutate(response = if_else(response == 1, "Yes", "No"),
         diet_var = unname(diet_labels[diet_var]),
         diet_var = factor(diet_var, 
                           levels = unname(diet_labels))
  ) %>%
  group_by(diet_var,
           response
  ) %>%
  summarise(n_participants = n(),
            mean_d13C_group = round(mean(mean_d13C, na.rm = TRUE), 2),
            sd_d13C = round(sd(mean_d13C, na.rm = TRUE), 2),
            median_d13C = round(median(mean_d13C, na.rm = TRUE), 2),
            IQR_d13C = round(IQR(mean_d13C, na.rm = TRUE), 2),
            .groups = "drop"
  )
