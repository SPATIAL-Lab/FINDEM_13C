

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



#unique states and number of moves per participant
mobility <- res %>%
  filter(!is.na(state), state != "NA") %>%        #filter out missing state, and "na"
  group_by(participant_id) %>%
  summarise(n_states = n_distinct(state),         #number of distinct states
            n_moves = n() - 1,                    #convert number of states to transitions
            .groups = "drop")                     #remove grouping structure


#-----------------------------------------------------------------------------------------------------------------------




#RESIDENTIAL MOBILITY


#prepare mobility data
joined_mob <- iso_pid %>%
  inner_join(mobility,
             by = "participant_id"
  )


#correlation between mean participant d13C and number of unique states lived in
cor_states <- cor.test(joined_mob$mean_d13C,
                       joined_mob$n_states
)


#correlation between mean participant d13C and number of residential moves
cor_moves <- cor.test(joined_mob$mean_d13C,
                      joined_mob$n_moves
)



#plot 2: d13C vs number of states lived in
p2 <- ggplot(joined_mob,
             aes(x = n_states,
                 y = mean_d13C)
) +
  geom_jitter(aes(color = cohort),           #individual participants
              width = 0.15,
              height = 0,
              alpha = 0.65,
              shape = 16
  ) +
  geom_smooth(method = "lm",                                 #linear regression line with 95% CI
              se = TRUE,
              color = "grey30",
              linewidth = 0.8,
              fill = "grey80",
              alpha = 0.3
  ) +
  annotate("text",                                                 #correlation annotation
           x = max(joined_mob$n_states,na.rm = TRUE) * 0.95,
           y = max(joined_mob$mean_d13C,na.rm = TRUE),
           label = paste0("r = ",round(cor_states$estimate,2),
                          "\np = ",round(cor_states$p.value,3)),
           hjust = 1,
           size = 3.5,
           color = "grey30"
  ) +
  scale_color_brewer(palette = "Dark2",
                     name = "Cohort"
  ) +
  scale_x_continuous(breaks = 1:max(joined_mob$n_states,na.rm = TRUE)
  ) +
  labs(title = "Mean δ¹³C vs Geographic Mobility",
       subtitle = "Each point = one participant | line = linear fit ± 95% CI",
       x = "Number of unique states lived in",
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40"),
        panel.grid.minor = element_blank()
  )


#save p2
ggsave(filename = file.path("outputs","d13C_vs_states.png"),
       plot = p2,
       width = 9,
       height = 6,
       dpi = 300
)




#plot2_2: d13C vs number of residential moves
p2_2 <- ggplot(joined_mob,
               aes(x = n_moves,
                   y = mean_d13C)
) +
  geom_jitter(aes(color = cohort),                         #individual participants
              width = 0.15,
              height = 0,
              alpha = 0.65,
              shape = 16
  ) +
  geom_smooth(method = "lm",            #linear regression line with 95% CI
              se = TRUE,
              color = "grey30",
              linewidth = 0.8,
              fill = "grey80",
              alpha = 0.3
  ) +
  annotate("text",                                                #correlation annotation
           x = max(joined_mob$n_moves,na.rm = TRUE) * 0.95,
           y = max(joined_mob$mean_d13C,na.rm = TRUE),
           label = paste0("r = ",round(cor_moves$estimate,2),
                          "\np = ",round(cor_moves$p.value,3)),
           hjust = 1,
           size = 3.5,
           color = "grey30"
  ) +
  scale_color_brewer(palette = "Dark2",
                     name = "Cohort"
  ) +
  scale_x_continuous(breaks = 0:max(joined_mob$n_moves,na.rm = TRUE)
  ) +
  labs(title = "Mean δ¹³C vs Residential Mobility",
       subtitle ="Each point = one participant | line = linear fit ± 95% CI",
       x = "Number of residential moves",
       y = "Mean participant δ¹³C (‰)"
  ) +
  theme_bw(base_size = 13
  ) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 10,
                                     color = "grey40"),
        panel.grid.minor = element_blank()
  )



#save p2_2
ggsave(filename = file.path("outputs", "d13C_vs_moves.png"),
       plot = p2_2,
       width = 9,
       height = 6,
       dpi = 300
)


