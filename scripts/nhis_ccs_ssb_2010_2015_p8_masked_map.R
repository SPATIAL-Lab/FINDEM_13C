#!/usr/bin/env Rscript

# Map CDC's published state estimates of daily sugar-sweetened-beverage (SSB)
# intake among US adults, from combined 2010 and 2015 NHIS Cancer Control
# Supplement (CCS) data.
#
# This is a reproduction of the published state table, not a re-analysis of
# public-use NHIS microdata. CDC's original state analysis required restricted
# state identifiers accessed through the NCHS Research Data Center.
#
# Measure: percentage of adults (aged >=18 y) reporting SSB intake >=1 time/day.
# SSBs include regular soda, sweetened fruit drinks, sports/energy drinks, and
# sweetened coffee/tea drinks. This is NOT a soda-only measure and NOT a mean
# number of drinks per person per day.
#
# Source: Chevinsky JR, Lee SH, Blanck HM, Park S. Prev Chronic Dis.
# 2021;18:200434. https://www.cdc.gov/pcd/issues/2021/20_0434.htm
# Published state table: https://www.cdc.gov/pcd/issues/2021/20_0434a.htm
#
# Outputs:
#   output/nhis_ccs_ssb_daily_prevalence_2010_2015_by_state.csv
#   output/nhis_ccs_ssb_daily_prevalence_2010_2015_p8_masked_warm_style_map.png

required_packages <- c("dplyr", "ggplot2", "readr", "usmap")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the required packages first:\ninstall.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))", call. = FALSE)
}

library(dplyr)
library(ggplot2)
library(readr)
library(usmap)

output_dir <- "outputs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Values transcribed from CDC's published tabular figure.
state_estimates <- read_csv(I(
  "fips,state,daily_ssb_prevalence\n01,Alabama,65.0\n02,Alaska,44.5\n04,Arizona,64.5\n05,Arkansas,74.2\n06,California,62.7\n08,Colorado,59.4\n09,Connecticut,72.2\n10,Delaware,68.0\n11,District of Columbia,64.8\n12,Florida,67.2\n13,Georgia,68.1\n15,Hawaii,76.4\n16,Idaho,58.8\n17,Illinois,62.7\n18,Indiana,65.7\n19,Iowa,50.5\n20,Kansas,54.9\n21,Kentucky,67.2\n22,Louisiana,68.7\n23,Maine,65.5\n24,Maryland,65.4\n25,Massachusetts,66.8\n26,Michigan,59.0\n27,Minnesota,50.4\n28,Mississippi,64.5\n29,Missouri,59.1\n30,Montana,64.9\n31,Nebraska,58.0\n32,Nevada,63.8\n33,New Hampshire,69.7\n34,New Jersey,69.5\n35,New Mexico,68.5\n36,New York,65.6\n37,North Carolina,62.7\n38,North Dakota,59.2\n39,Ohio,57.2\n40,Oklahoma,66.0\n41,Oregon,51.5\n42,Pennsylvania,65.9\n44,Rhode Island,65.7\n45,South Carolina,70.2\n46,South Dakota,72.5\n47,Tennessee,66.4\n48,Texas,62.5\n49,Utah,53.6\n50,Vermont,67.3\n51,Virginia,59.6\n53,Washington,55.0\n54,West Virginia,59.4\n55,Wisconsin,50.4\n56,Wyoming,73.2"
), col_types = cols(
  fips = col_character(), state = col_character(), daily_ssb_prevalence = col_double()
), show_col_types = FALSE)

stopifnot(nrow(state_estimates) == 51L, !anyDuplicated(state_estimates$fips))

# States shown in gray in p8_present.png.  The mask is only for visual
# comparability: CDC published NHIS estimates for every one of these states.
p8_masked_fips <- c("02", "05", "09", "23", "32", "34", "38", "46", "56")

map_values <- state_estimates |>
  mutate(
    daily_ssb_prevalence_for_map = if_else(
      fips %in% p8_masked_fips, NA_real_, daily_ssb_prevalence
    )
  )

# The exported estimates remain complete and unmasked.
write_csv(
  state_estimates,
  file.path(output_dir, "nhis_ccs_ssb_daily_prevalence_2010_2015_by_state.csv")
)

map <- plot_usmap(
  regions = "states",
  data = map_values,
  values = "daily_ssb_prevalence_for_map",
  color = "black",
  linewidth = 0.45
) +
  # Palette and boundary treatment match the supplied reference figure.
  scale_fill_gradientn(
    colours = c("#EDD6A6", "#E4B640", "#E8A400", "#D5660C", "#9C2527"),
    name = "Adults reporting SSB intake ≥1 time/day (%)",
    na.value = "#888888",
    labels = scales::label_number(accuracy = 1, suffix = "%"),
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = grid::unit(8, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  ) +
  labs(
    title = "Daily sugar-sweetened beverage intake among US adults",
    subtitle = "Published state estimates from combined 2010 and 2015 NHIS Cancer Control Supplement data",
    caption = paste(
      "Gray = states masked for visual comparability with P8_present.png (AK, AR, CT, ME, NV, NJ, ND, SD, WY), not missing NHIS data.",
      "Outcome: percent of adults reporting total SSB intake ≥1 time/day. SSBs include regular soda, sweetened fruit drinks,",
      "sports/energy drinks, and sweetened coffee/tea drinks; this is not a soda-only measure or a common single-year estimate.",
      "Values reproduce CDC's published state table (Chevinsky et al., 2021; combined n = 56,260).",
      "CDC used restricted NHIS state identifiers; this script maps the published estimates and does not recalculate them from public-use microdata.",
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

message("Wrote the published NHIS CCS state estimates and P8-masked map to ", output_dir)
