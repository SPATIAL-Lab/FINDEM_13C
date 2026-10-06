#!/usr/bin/env Rscript

# Aggregate the county-level fast-food outcome underlying Figure 2 of:
# Althoff T, Nilforoshan H, Hua J, Leskovec J. Large-scale diet tracking data
# reveal disparate associations between food environment and diet. Nat Commun.
# 2022;13:267. https://doi.org/10.1038/s41467-021-27522-y
#
# The paper's outcome is the mean number of FAST-FOOD FOOD-LOG ENTRIES per
# participant per active day, not the number of meals. It is derived from
# MyFitnessPal users observed from 2010-2016 and is not population-representative.
# Figure 2 displays counties with more than 30 participants. This script uses
# those same county-level data and produces a participant-weighted state mean.
#
# Outputs:
#   output/dietdisparities_fast_food_entries_2010_2016_by_state.csv
#   output/dietdisparities_fast_food_entries_2010_2016_p8_masked_state_map.png

required_packages <- c("curl", "dplyr", "ggplot2", "readr", "stringr", "usmap")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the required packages first:\ninstall.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))", call. = FALSE)
}

library(curl)
library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(usmap)

archive_url <- "https://snap.stanford.edu/dietdisparities/DietDisparities.zip"
data_dir <- "data/dietdisparities"
output_dir <- "outputs"
archive_path <- file.path(data_dir, "DietDisparities.zip")
extract_dir <- file.path(data_dir, "DietDisparities")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# The archive is the authors' public data/code release. A locally downloaded
# archive at archive_path is reused, which also provides a manual fallback if a
# network or institutional firewall blocks the Stanford download.
if (!file.exists(archive_path)) {
  download_handle <- new_handle()
  handle_setheaders(
    download_handle,
    "User-Agent" = paste(
      "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
    )
  )
  response <- curl_fetch_disk(archive_url, path = archive_path, handle = download_handle)
  if (response$status_code != 200L) {
    unlink(archive_path)
    stop(
      "The Stanford archive download returned HTTP ", response$status_code, ". ",
      "Download DietDisparities.zip from https://snap.stanford.edu/dietdisparities/ ",
      "and save it as ", archive_path, ", then run this script again.",
      call. = FALSE
    )
  }
}

if (!dir.exists(extract_dir)) {
  dir.create(extract_dir, recursive = TRUE)
  unzip(archive_path, exdir = extract_dir)
}

# The archive contains ZIP- and county-level files. Figure 2 uses the county
# outcome data, so select exactly one CSV whose name identifies it as county data.
csv_files <- list.files(extract_dir, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
county_files <- csv_files[str_detect(basename(csv_files), regex("county", ignore_case = TRUE))]
if (length(county_files) != 1L) {
  stop(
    "Could not identify exactly one county-level CSV in the archive. Candidates: ",
    paste(county_files, collapse = "; "), call. = FALSE
  )
}

county_data <- read_csv(county_files, show_col_types = FALSE, name_repair = "universal")
normal_names <- names(county_data) |>
  str_to_lower() |>
  str_replace_all("[^a-z0-9]+", "_") |>
  str_replace_all("(^_|_$)", "")

county_fips_col <- "STCOUNTYFP"
fast_food_col <- "fastfood"
participant_col <- "count"

county_values <- county_data |>
  transmute(
    county_fips = str_pad(as.character(.data[[county_fips_col]]), width = 5, pad = "0"),
    fips = str_sub(county_fips, 1, 2),
    fast_food_entries_per_participant_per_day = as.numeric(.data[[fast_food_col]]),
    n_participants = as.numeric(.data[[participant_col]])
  ) |>
  filter(
    str_detect(county_fips, "^[0-9]{5}$"),
    n_participants > 30,
    !is.na(fast_food_entries_per_participant_per_day)
  )

if (!nrow(county_values)) {
  stop(
    "No county records remained after applying Figure 2's >30-participant threshold. ",
    "Inspect the detected columns and the archive version.", call. = FALSE
  )
}

# Weight county means by their Figure 2 participant counts. This gives counties
# with more logged participants proportionally more influence; it is not a
# survey-weighted or population-representative estimate.
state_estimates <- county_values |>
  group_by(fips) |>
  summarise(
    n_counties = n(),
    fast_food_entries_per_participant_per_day = weighted.mean(
      fast_food_entries_per_participant_per_day, w = n_participants
    ),
    n_participants = sum(n_participants),
    .groups = "drop"
  ) |>
  arrange(fips)

# States shown in gray in p8_present.png, retained only for visual comparison.
p8_masked_fips <- c("02", "05", "09", "23", "32", "34", "38", "46", "56")
map_values <- state_estimates |>
  mutate(
    fast_food_entries_for_map = if_else(
      fips %in% p8_masked_fips, NA_real_, fast_food_entries_per_participant_per_day
    )
  )

# The exported table remains complete and unmasked.
write_csv(
  state_estimates,
  file.path(output_dir, "dietdisparities_fast_food_entries_2010_2016_by_state.csv")
)

map <- plot_usmap(
  regions = "states",
  data = map_values,
  values = "fast_food_entries_for_map",
  color = "black",
  linewidth = 0.45
) +
  scale_fill_gradientn(
    colours = c("#EDD6A6", "#E4B640", "#E8A400", "#D5660C", "#9C2527"),
    name = "Fast-food entries\nper participant per day",
    na.value = "#888888",
    labels = scales::label_number(accuracy = 0.01),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = grid::unit(8, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  ) +
  labs(
    title = "Fast-food entries logged by MyFitnessPal participants",
    subtitle = "Participant-weighted state aggregation of Figure 2 county outcomes (2010–2016)",
    caption = paste(
      "Gray = states masked for visual comparability with P8_present.png (AK, AR, CT, ME, NV, NJ, ND, SD, WY), not necessarily missing data.",
      "The paper's measure is mean fast-food food-log entries per participant per active day, not meals or servings.",
      "Counties with >30 participants are included; state means are weighted by county participant count and are not population-representative estimates.",
      "Source: Althoff et al. (2022), Nature Communications 13:267; authors' public DietDisparities archive.",
      sep = "\n"
    )
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    plot.subtitle = element_text(size = 11, color = "black"),
    plot.caption = element_text(size = 8, hjust = 0, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(size = 11, hjust = 0.5, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.margin = margin(t = 4, b = 2),
    plot.margin = margin(10, 18, 10, 18)
  )

ggsave(
  file.path(output_dir, "p8b_present.png"),
  map, width = 11, height = 7, dpi = 300, bg = "white"
)

message("Used county file: ", county_files)
message("Wrote state table and P8-masked map to ", output_dir)

# Clean up
unlink(data_dir, recursive = TRUE)
