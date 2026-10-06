#!/usr/bin/env Rscript

# Download the USDA ERS Feed Grains Yearbook Tables CSV and extract annual
# Annual per-capita corn use for high-fructose corn syrup, glucose and dextrose, starch,
# and cereals and foods/other products (Yearbook Table 31). The output has one
# row per marketing year and one million-bushel column for each use category.
# The analysis is capped at marketing year 2025. It downloads Census Bureau population estimates and adds annual average
# population and bushels-per-person columns for the corresponding marketing year.
#
# Output:
#   output/usda_ers_corn_selected_food_industrial_uses_per_capita.csv
#
# Source: U.S. Department of Agriculture, Economic Research Service,
# Feed Grains Database, Feed Grains: Yearbook Tables.
# https://www.ers.usda.gov/data-products/feed-grains-database/feed-grains-yearbook-tables
# Population source: U.S. Census Bureau, National Population Estimates,
# distributed through FRED series POP (monthly, thousands of persons).
# https://fred.stlouisfed.org/series/POP

required_packages <- c("curl", "dplyr", "readr", "rvest", "stringr", "tidyr", "xml2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install the required packages first:\ninstall.packages(c(\"",
       paste(missing_packages, collapse = "\", \""), "\"))", call. = FALSE)
}

library(curl)
library(dplyr)
library(readr)
library(rvest)
library(stringr)
library(tidyr)
library(xml2)

page_url <- paste0(
  "https://www.ers.usda.gov/data-products/feed-grains-database/",
  "feed-grains-yearbook-tables"
)
data_dir <- "data/usda_ers_feed_grains"
output_dir <- "outputs"
# The ERS Yearbook can include projected future marketing years. Limit the
# per-capita analysis to the latest completed population year.
analysis_end_year <- 2025L
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use a standard browser user agent because the ERS site may reject bare
# command-line HTTP requests. The CSV URL itself is discovered from the
# landing page, avoiding a hard-coded, versioned media URL.
ers_handle <- new_handle()
handle_setheaders(
  ers_handle,
  "User-Agent" = paste(
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
    "(KHTML, like Gecko) Chrome/120.0 Safari/537.36"
  ),
  "Referer" = page_url
)

get_csv_url <- function(page_url, handle) {
  response <- curl_fetch_memory(page_url, handle = handle)
  if (response$status_code != 200L) {
    stop("ERS landing page returned HTTP ", response$status_code, ".", call. = FALSE)
  }

  page <- read_html(response$content)
  hrefs <- page |>
    html_elements("a") |>
    html_attr("href") |>
    na.omit()
  csv_hrefs <- hrefs[str_detect(hrefs, regex("\\.csv(?:[?]|$)", ignore_case = TRUE))]

  if (length(csv_hrefs) != 1L) {
    stop(
      "Could not identify exactly one Yearbook CSV link on the ERS page. ",
      "Found: ", paste(csv_hrefs, collapse = "; "), call. = FALSE
    )
  }
  url_absolute(csv_hrefs, page_url)
}

csv_url <- get_csv_url(page_url, ers_handle)
raw_csv_path <- file.path(data_dir, "feed_grains_yearbook_tables_all_years.csv")

download_response <- curl_fetch_disk(csv_url, path = raw_csv_path, handle = ers_handle)
if (download_response$status_code != 200L) {
  unlink(raw_csv_path)
  stop(
    "ERS CSV download returned HTTP ", download_response$status_code, ". ",
    "Try opening the Yearbook Tables page in a browser, then run this script again.",
    call. = FALSE
  )
}

raw <- read_csv(raw_csv_path, show_col_types = FALSE, name_repair = "universal")

# The corn marketing year runs September-August. Use the mean of Census
# population estimates across those 12 months as the denominator, so the
# population reference period matches each marketing-year numerator.
population_url <- "https://fred.stlouisfed.org/graph/fredgraph.csv?id=POP"
population_path <- file.path(data_dir, "us_census_population_fred_pop.csv")
population_handle <- new_handle()
handle_setheaders(
  population_handle,
  "User-Agent" = "Mozilla/5.0 (compatible; USDA-ERS-corn-use-script/1.0)"
)
population_response <- curl_fetch_disk(
  population_url, path = population_path, handle = population_handle
)
if (population_response$status_code != 200L) {
  unlink(population_path)
  stop("Population download returned HTTP ", population_response$status_code, ".", call. = FALSE)
}

population_raw <- read_csv(population_path, show_col_types = FALSE)
if (!all(c("observation_date", "POP") %in% names(population_raw))) {
  stop(
    "The population download did not contain the expected observation_date and POP columns.",
    call. = FALSE
  )
}
population_by_marketing_year <- population_raw |>
  transmute(
    date = as.Date(observation_date),
    population_people = as.numeric(POP) * 1000,
    population_year = if_else(
      as.integer(format(date, "%m")) >= 9L,
      as.integer(format(date, "%Y")),
      as.integer(format(date, "%Y")) - 1L
    )
  ) |>
  filter(!is.na(population_people)) |>
  group_by(population_year) |>
  summarise(
    us_population_people = mean(population_people),
    months_in_population_average = n(),
    .groups = "drop"
  )

# Find required fields defensively; ERS can rename fields while retaining the
# same data. The script stops instead of silently extracting an unintended row.
normal_names <- names(raw) |>
  str_to_lower() |>
  str_replace_all("[^a-z0-9]+", "_") |>
  str_replace_all("(^_|_$)", "")

find_column <- function(candidates, label) {
  matches <- which(normal_names %in% candidates)
  if (length(matches) != 1L) {
    stop(
      "Could not uniquely identify the ", label, " column. Available columns: ",
      paste(names(raw), collapse = ", "), call. = FALSE
    )
  }
  names(raw)[matches]
}

commodity_col <- find_column(c("commodity", "fgyt_commodity"), "commodity")
attribute_col <- find_column(c("attribute", "fgyt_attribute", "data_item", "item"), "attribute")
year_col <- find_column(c("year", "market_year", "marketing_year", "fgyt_year"), "year")
value_col <- find_column(c("value", "fgyt_value", "amount"), "value")
table_col <- {
  matches <- which(normal_names %in% c(
    "table", "table_number", "table_no", "table_name", "fgyt_table", "fgyt_table_name"
  ))
  if (length(matches) == 1L) names(raw)[matches] else NULL
}
timeperiod_col <- find_column(
  c("timeperiod", "fgyt_timeperiod", "period_desc", "period_description"),
  "time-period"
)

corn <- raw |>
  filter(str_detect(.data[[commodity_col]], regex("^corn$", ignore_case = TRUE)))

# Table 31 contains the detailed food, seed, and industrial uses of corn.
if (!is.null(table_col)) {
  corn <- corn |>
    filter(str_detect(as.character(.data[[table_col]]), "(^|[^0-9])31([^0-9]|$)"))
}

corn <- corn |>
  mutate(
    use_category = case_when(
      str_detect(.data[[attribute_col]], regex("high[- ]fructose corn syrup|hfcs", ignore_case = TRUE)) ~
        "High-fructose corn syrup",
      str_detect(.data[[attribute_col]], regex("glucose and dextrose", ignore_case = TRUE)) ~
        "Glucose and dextrose",
      str_detect(.data[[attribute_col]], regex("^starch( use)?$", ignore_case = TRUE)) ~
        "Starch",
      str_detect(.data[[attribute_col]], regex("cereals? and (foods?|other products)", ignore_case = TRUE)) ~
        "Cereals and foods",
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    !is.na(use_category),
    str_detect(
      .data[[timeperiod_col]],
      regex("^marketing year\\s+sep[-–]aug", ignore_case = TRUE)
    )
  )

if (!nrow(corn)) {
  available_attributes <- raw |>
    filter(str_detect(.data[[commodity_col]], regex("^corn$", ignore_case = TRUE))) |>
    pull(.data[[attribute_col]]) |>
    unique()
  stop(
    "No requested corn-use rows were found. Corn attributes available in this download: ",
    paste(available_attributes, collapse = "; "), call. = FALSE
  )
}

result_long <- corn |>
  transmute(
    Year = as.character(.data[[year_col]]),
    use_category,
    use_column = case_when(
      use_category == "High-fructose corn syrup" ~ "high_fructose_corn_syrup_million_bushels",
      use_category == "Glucose and dextrose" ~ "glucose_and_dextrose_million_bushels",
      use_category == "Starch" ~ "starch_million_bushels",
      use_category == "Cereals and foods" ~ "cereals_and_foods_million_bushels"
    ),
    corn_use = suppressWarnings(as.numeric(.data[[value_col]]))
  ) |>
  arrange(Year, factor(
    use_category,
    levels = c("High-fructose corn syrup", "Glucose and dextrose", "Starch", "Cereals and foods")
  ))

expected_categories <- c(
  "High-fructose corn syrup", "Glucose and dextrose", "Starch", "Cereals and foods"
)
missing_categories <- setdiff(expected_categories, unique(result_long$use_category))
if (length(missing_categories)) {
  stop(
    "The ERS download did not contain all requested categories. Missing: ",
    paste(missing_categories, collapse = ", "), call. = FALSE
  )
}

duplicate_rows <- result_long |>
  count(Year, use_column) |>
  filter(n > 1L)
if (nrow(duplicate_rows)) {
  stop("More than one annual record was found for at least one year/category.", call. = FALSE)
}

result <- result_long |>
  select(Year, use_column, corn_use) |>
  pivot_wider(names_from = use_column, values_from = corn_use) |>
  mutate(population_year = as.integer(str_extract(Year, "^[0-9]{4}"))) |>
  left_join(population_by_marketing_year, by = "population_year") |>
  filter(population_year <= analysis_end_year) |>
  arrange(Year)

if (anyNA(result[, c(
  "high_fructose_corn_syrup_million_bushels",
  "glucose_and_dextrose_million_bushels",
  "starch_million_bushels",
  "cereals_and_foods_million_bushels"
)])) {
  stop("At least one annual total is missing after reshaping the requested categories.", call. = FALSE)
}

if (anyNA(result$us_population_people) || any(result$months_in_population_average != 12L)) {
  missing_or_incomplete <- result |>
    filter(is.na(us_population_people) | months_in_population_average != 12L) |>
    pull(Year)
  stop(
    "A complete 12-month population denominator was not available for: ",
    paste(missing_or_incomplete, collapse = ", "), call. = FALSE
  )
}

if (!nrow(result)) {
  stop("No rows remain after applying analysis_end_year = ", analysis_end_year, ".", call. = FALSE)
}

result <- result |>
  mutate(
    across(
      ends_with("_million_bushels"),
      ~ .x * 1e6 / us_population_people,
      .names = "{.col}_per_capita_bushels"
    )
  ) |>
  select(-population_year, -months_in_population_average)

write_csv(result, file.path(output_dir, "usda_ers_corn_selected_food_industrial_uses_per_capita.csv"))
message("Downloaded: ", csv_url)
message("Downloaded population: ", population_url)
message("Wrote: ", file.path(output_dir, "usda_ers_corn_selected_food_industrial_uses_per_capita.csv"))

# Clean up
unlink(data_dir, recursive = TRUE)
