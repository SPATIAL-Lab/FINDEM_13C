

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
#------------------------------------------------------------------------------------------------------------------------

#GEOGRAPHIC

#keep only participants in res_longest and iso_pid, organize by id
combined <- iso_pid %>%
  inner_join(res_longest, by = "participant_id")

#Mean participant-level d13C for each state, each participant weighted equally regardless of number of samples
state_d13C <- combined %>%
  filter(!is.na(state), !state %in% c("Puerto Rico", "Guam")) %>%
  group_by(state) %>%                                                     #summarize isotope values by state
  summarise(state_mean_d13C = mean(mean_d13C, na.rm = TRUE),
            state_median_d13C = median(mean_d13C, na.rm = TRUE),
            state_sd_d13C = sd(mean_d13C, na.rm = TRUE),
            n_participants = n(),
            .groups = "drop") %>%
  filter(n_participants >= 3)                                          #filter to keep only states with enough data


#δ¹³C values grouped by 3-digit ZIP
zip_d13C <- combined %>%
  filter(!is.na(zipthree)) %>%
  group_by(zipthree) %>%
  summarise(zip_mean_d13C = mean(mean_d13C, na.rm = TRUE),
            n_participants = n(),
            .groups = "drop") %>%
  filter(n_participants >= 2)


#ZIP prefix centroids from zipcodeR
zip_centroids <- zipcodeR::zip_code_db %>%
  as_tibble() %>%                                     #dataframe
  select(zipcode, lat, lng) %>%
  filter(!is.na(lat), !is.na(lng)) %>%
  mutate(zip3 = str_sub(zipcode, 1, 3)) %>%           #use first 3 digits of zip
  group_by(zip3) %>%
  summarise(lat = mean(lat, na.rm = TRUE),            #average zip lat
            lon = mean(lng, na.rm = TRUE),            #average zip long
            .groups = "drop")



#join and project to usmap coordinates
zip_plot_data <- zip_d13C %>%
  left_join(zip_centroids, by = c("zipthree" = "zip3")) %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  as.data.frame() %>%
  usmap::usmap_transform(input_names = c("lon", "lat"))%>%
  mutate(x = st_coordinates(geometry)[, 1],
         y = st_coordinates(geometry)[, 2]) %>%
  st_drop_geometry()




#map1 
#state choropleth, mean δ¹³C per state
map1 <- plot_usmap(data = state_d13C,
                   values = "state_mean_d13C",
                   color = "black",
                   linewidth = 0.3,
                   labels = FALSE) +
  scale_fill_gradientn(colors = c("wheat2", "tan", "goldenrod", "darkorange3", "chocolate", "brown4"),
                       na.value = "grey45",
                       name = "Mean δ¹³C (‰)",
                       guide = guide_colorbar(barwidth = 12, 
                                              barheight = 0.8,
                                              title.position = "top", 
                                              title.hjust = 0.5)) +
  labs(title = "Mean δ¹³C by State of Longest Residence",
       subtitle = "Grey = insufficient data (<3 participants)")+
  theme(plot.title = element_text(face = "bold", 
                                  size = 14, 
                                  hjust = 0.6,
                                  color = "black"),
        plot.subtitle = element_text(size = 10, 
                                     color = "grey45", 
                                     hjust = 0.6),
        legend.position = "bottom",
        legend.direction = "horizontal"
        )


#save map1
ggsave(
  filename = file.path("outputs", "map1_state_d13C.png"),
  plot = map1,
  width = 11,
  height = 7,
  dpi = 300
)


#map2
#map aggregated by 3-digit zip prefix
us_base <- plot_usmap(color = "black", 
                      linewidth = 0.3, 
                      fill = "grey30") +
  theme(panel.background = element_rect(fill = "white", color = NA))

map2 <- us_base +
  geom_point(data = zip_plot_data,
             aes(x = x, 
                 y = y,
                 color = zip_mean_d13C,
                 size = n_participants),
             alpha = 0.80,
             shape = 16) +
  scale_color_gradientn(colors = c("wheat2", "tan", "goldenrod", "darkorange3", "chocolate", "red4"),
                        name = "Mean δ¹³C (‰)",
                        guide = guide_colorbar(barwidth = 10, barheight = 0.8,
                                               title.position = "top", title.hjust = 0.5)) +
  scale_size_continuous(name = "Participants (n)",
                        range = c(2, 10),
                        breaks = c(2, 5, 10, 20)) +
  labs(title = "δ¹³C by 3-Digit ZIP Prefix of Longest Residence",
       subtitle = paste0("ZIP prefixes with ≥2 participants (n prefixes = ", nrow(zip_plot_data), ")"),
       caption = "Bubble color = mean δ¹³C; size = participant count") +
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 2, override.aes = list(color = "black"))) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40"),
        plot.caption = element_text(size = 8,  color = "grey50"),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal")


ggsave(
  filename = file.path("outputs", "map2_zip_d13C.png"),
  plot = map2,
  width = 11,
  height = 7,
  dpi = 300
)







