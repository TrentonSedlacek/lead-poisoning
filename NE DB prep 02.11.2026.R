library(tidyverse)
library(sf)
library(tigris)
library(tidycensus)

options(tigris_use_cache = TRUE)
select <- dplyr::select
filter <- dplyr::filter

DATA_PATH   <- "K:/CLPPP/TrentonS"
OUTPUT_PATH <- file.path(DATA_PATH, "shiny_data")
BLL_CUTOFF  <- 3.5
dir.create(OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)

# --- LHD-County Mapping ---

lhd_map <- tribble(
  ~county, ~lhd,
  "Adams", "South Heartland District", "Antelope", "North Central District",
  "Arthur", "West Central District", "Banner", "Panhandle", "Blaine", "Loup Basin",
  "Boone", "East Central District", "Box Butte", "Panhandle", "Boyd", "North Central District",
  "Brown", "North Central District", "Buffalo", "Two Rivers", "Burt", "Elkhorn Logan Valley",
  "Butler", "Four Corners", "Cass", "Sarpy/Cass", "Cedar", "Northeast Nebraska",
  "Chase", "Southwest Nebraska", "Cherry", "North Central District", "Cheyenne", "Panhandle",
  "Clay", "South Heartland District", "Colfax", "East Central District",
  "Cuming", "Elkhorn Logan Valley", "Custer", "Loup Basin", "Dakota", "Dakota County",
  "Dawes", "Panhandle", "Dawson", "Two Rivers", "Deuel", "Panhandle",
  "Dixon", "Northeast Nebraska", "Dodge", "Three Rivers", "Douglas", "Douglas County",
  "Dundy", "Southwest Nebraska", "Fillmore", "Public Health Solutions",
  "Franklin", "Two Rivers", "Frontier", "Southwest Nebraska", "Furnas", "Southwest Nebraska",
  "Gage", "Public Health Solutions", "Garden", "Panhandle", "Garfield", "Loup Basin",
  "Gosper", "Two Rivers", "Grant", "Panhandle", "Greeley", "Loup Basin",
  "Hall", "Central District", "Hamilton", "Central District", "Harlan", "Two Rivers",
  "Hayes", "Southwest Nebraska", "Hitchcock", "Southwest Nebraska",
  "Holt", "North Central District", "Hooker", "West Central District", "Howard", "Loup Basin",
  "Jefferson", "Public Health Solutions", "Johnson", "Southeast District",
  "Kearney", "Two Rivers", "Keith", "Southwest Nebraska", "Keya Paha", "North Central District",
  "Kimball", "Panhandle", "Knox", "North Central District", "Lancaster", "Lincoln-Lancaster County",
  "Lincoln", "West Central District", "Logan", "West Central District", "Loup", "Loup Basin",
  "Madison", "Elkhorn Logan Valley", "McPherson", "West Central District",
  "Merrick", "Central District", "Morrill", "Panhandle", "Nance", "East Central District",
  "Nemaha", "Southeast District", "Nuckolls", "South Heartland District",
  "Otoe", "Southeast District", "Pawnee", "Southeast District", "Perkins", "Southwest Nebraska",
  "Phelps", "Two Rivers", "Pierce", "North Central District", "Platte", "East Central District",
  "Polk", "Four Corners", "Red Willow", "Southwest Nebraska", "Richardson", "Southeast District",
  "Rock", "North Central District", "Saline", "Public Health Solutions", "Sarpy", "Sarpy/Cass",
  "Saunders", "Three Rivers", "Scotts Bluff", "Panhandle", "Seward", "Four Corners",
  "Sheridan", "Panhandle", "Sherman", "Loup Basin", "Sioux", "Panhandle",
  "Stanton", "Elkhorn Logan Valley", "Thayer", "Public Health Solutions",
  "Thomas", "West Central District", "Thurston", "Northeast Nebraska", "Valley", "Loup Basin",
  "Washington", "Three Rivers", "Wayne", "Northeast Nebraska", "Webster", "South Heartland District",
  "Wheeler", "Loup Basin", "York", "Four Corners"
)
saveRDS(lhd_map, file.path(OUTPUT_PATH, "lhd_map.rds"))

# --- Load and Clean Lead Data ---

raw <- read_csv(
  file.path(DATA_PATH, "MatchedGeo2010-2023.csv"),
  col_types = cols(.default = col_guess()),
  show_col_types = FALSE
)

lead <- raw %>%
  filter(
    Age_yr <= 6,
    !is.na(Latitude), Latitude != 0,
    !is.na(Longitude),
    FeatureMatchingResultType != "Unmatchable",
    !is.na(result),
    !is.na(PATIENT_LOCAL_ID),
    !is.na(CensusTract2020), as.numeric(CensusTract2020) > 0,
    !is.na(CensusCountyFips2020),
    !is.na(CensusBlockGroup2020)
  ) %>%
  mutate(
    elevated = as.integer(result >= BLL_CUTOFF),
    county = str_to_title(str_trim(county)),
    county = if_else(county == "Mcpherson", "McPherson", county),
    street_address = paste(Street, City, State, Zip, sep = ", "),

    # Build GEOIDs from component columns (NOT GeoLocationID2020, which is
    # corrupted by scientific notation in the CSV — e.g. "3.11E+11" loses digits).
    # CensusTract2020 is a decimal like 105.05 or integer like 40.
    # Convert to 6-digit TRACTCE: floor gives base, remainder * 100 gives suffix.
    tract_num = as.numeric(CensusTract2020),
    tract_int = floor(tract_num),
    tract_dec = round((tract_num - tract_int) * 100),
    TRACTCE     = sprintf("%04d%02d", tract_int, tract_dec),
    tract_geoid = paste0("31", sprintf("%03d", as.integer(CensusCountyFips2020)), TRACTCE),
    bg_geoid    = paste0(tract_geoid, as.integer(CensusBlockGroup2020))
  ) %>%
  filter(
    nchar(tract_geoid) == 11,
    nchar(bg_geoid) == 12
  ) %>%
  select(-tract_num, -tract_int, -tract_dec, -TRACTCE) %>%
  left_join(lhd_map, by = "county")

# --- First Test Per Child ---

first_test <- lead %>%
  arrange(PATIENT_LOCAL_ID, sample_year) %>%
  group_by(PATIENT_LOCAL_ID) %>%
  slice(1) %>%
  ungroup()
saveRDS(first_test, file.path(OUTPUT_PATH, "first_test.rds"))

# --- Block Group Stats ---

bg_stats <- first_test %>%
  filter(!is.na(bg_geoid)) %>%
  group_by(bg_geoid) %>%
  summarize(
    n_children = n(),
    n_elevated = sum(elevated, na.rm = TRUE),
    mean_bll   = mean(result, na.rm = TRUE),
    max_bll    = max(result, na.rm = TRUE),
    county     = first(county),
    lhd        = first(lhd),
    .groups = "drop"
  ) %>%
  mutate(pct_elevated = round(100 * n_elevated / n_children, 2))

bg_addr <- first_test %>%
  filter(!is.na(bg_geoid)) %>%
  group_by(bg_geoid, address_id) %>%
  summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
  group_by(bg_geoid) %>%
  summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")

bg_stats <- left_join(bg_stats, bg_addr, by = "bg_geoid")
saveRDS(bg_stats, file.path(OUTPUT_PATH, "bg_stats.rds"))

# --- Tract Stats ---

tract_stats <- first_test %>%
  filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid) %>%
  summarize(
    n_children = n(),
    n_elevated = sum(elevated, na.rm = TRUE),
    mean_bll   = mean(result, na.rm = TRUE),
    max_bll    = max(result, na.rm = TRUE),
    county     = first(county),
    lhd        = first(lhd),
    .groups = "drop"
  ) %>%
  mutate(pct_elevated = round(100 * n_elevated / n_children, 2))

tract_addr <- first_test %>%
  filter(!is.na(tract_geoid)) %>%
  group_by(tract_geoid, address_id) %>%
  summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
  group_by(tract_geoid) %>%
  summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")

tract_stats <- left_join(tract_stats, tract_addr, by = "tract_geoid")
saveRDS(tract_stats, file.path(OUTPUT_PATH, "tract_stats.rds"))

# --- County Stats ---

county_stats <- first_test %>%
  group_by(county) %>%
  summarize(
    n_children = n(),
    n_elevated = sum(elevated, na.rm = TRUE),
    mean_bll   = mean(result, na.rm = TRUE),
    max_bll    = max(result, na.rm = TRUE),
    lhd        = first(lhd),
    .groups = "drop"
  ) %>%
  mutate(pct_elevated = round(100 * n_elevated / n_children, 2))

county_addr <- first_test %>%
  group_by(county, address_id) %>%
  summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
  group_by(county) %>%
  summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")

county_stats <- left_join(county_stats, county_addr, by = "county")
saveRDS(county_stats, file.path(OUTPUT_PATH, "county_stats.rds"))

# --- LHD Stats ---

lhd_stats <- first_test %>%
  filter(!is.na(lhd)) %>%
  group_by(lhd) %>%
  summarize(
    n_children = n(),
    n_elevated = sum(elevated, na.rm = TRUE),
    mean_bll   = mean(result, na.rm = TRUE),
    max_bll    = max(result, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pct_elevated = round(100 * n_elevated / n_children, 2))

lhd_addr <- first_test %>%
  filter(!is.na(lhd)) %>%
  group_by(lhd, address_id) %>%
  summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
  group_by(lhd) %>%
  summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")

lhd_counties <- lhd_map %>% count(lhd, name = "n_counties")
lhd_stats <- lhd_stats %>%
  left_join(lhd_addr, by = "lhd") %>%
  left_join(lhd_counties, by = "lhd")
saveRDS(lhd_stats, file.path(OUTPUT_PATH, "lhd_stats.rds"))

# --- Address Stats ---
# For addresses with 2+ children tested

address_stats <- first_test %>%
  group_by(address_id) %>%
  summarize(
    n_children     = n(),
    n_elevated     = sum(elevated, na.rm = TRUE),
    mean_bll       = mean(result, na.rm = TRUE),
    max_bll        = max(result, na.rm = TRUE),
    first_bll      = result[which.min(sample_year)[1]],
    first_year     = min(sample_year, na.rm = TRUE),
    last_year      = max(sample_year, na.rm = TRUE),
    lat            = first(Latitude),
    lng            = first(Longitude),
    county         = first(county),
    lhd            = first(lhd),
    bg_geoid       = first(bg_geoid),
    tract_geoid    = first(tract_geoid),
    street_address = first(street_address),
    .groups = "drop"
  ) %>%
  filter(n_children >= 2)
saveRDS(address_stats, file.path(OUTPUT_PATH, "address_stats.rds"))

# --- Model Data ---
# First child at each address, then check if subsequent children were elevated

first_child <- first_test %>%
  arrange(address_id, sample_year) %>%
  group_by(address_id) %>%
  slice(1) %>%
  ungroup() %>%
  select(address_id, first_id = PATIENT_LOCAL_ID, first_year = sample_year,
         first_bll = result, bg_geoid, tract_geoid, county, lhd,
         lat = Latitude, lng = Longitude, street_address)

subsequent <- first_test %>%
  inner_join(first_child %>% select(address_id, first_id, first_year), by = "address_id") %>%
  filter(PATIENT_LOCAL_ID != first_id, sample_year >= first_year) %>%
  group_by(address_id) %>%
  summarize(n_subsequent = n(), n_sub_elevated = sum(elevated, na.rm = TRUE), .groups = "drop")

model_data <- first_child %>%
  inner_join(subsequent, by = "address_id") %>%
  filter(!is.na(bg_geoid), nchar(bg_geoid) == 12) %>%
  mutate(outcome = as.integer(n_sub_elevated > 0))
saveRDS(model_data, file.path(OUTPUT_PATH, "model_data.rds"))

# --- Geographic Boundaries ---

ne_counties <- counties(state = "NE", cb = TRUE, year = 2021) %>% st_transform(4326)
ne_tracts   <- tracts(state = "NE", cb = TRUE, year = 2021)   %>% st_transform(4326)
ne_bgs      <- block_groups(state = "NE", cb = TRUE, year = 2021) %>% st_transform(4326)

lhd_bounds <- ne_counties %>%
  left_join(lhd_map, by = c("NAME" = "county")) %>%
  filter(!is.na(lhd)) %>%
  group_by(lhd) %>%
  summarize(geometry = st_union(geometry), .groups = "drop")

saveRDS(ne_counties, file.path(OUTPUT_PATH, "ne_counties.rds"))
saveRDS(ne_tracts,   file.path(OUTPUT_PATH, "ne_tracts.rds"))
saveRDS(ne_bgs,      file.path(OUTPUT_PATH, "ne_bgs.rds"))
saveRDS(lhd_bounds,  file.path(OUTPUT_PATH, "lhd_bounds.rds"))

# --- Census (ACS 2018-2022) ---

census_vars <- c(
  tot_hu = "B25034_001", pre39 = "B25034_011", y40_49 = "B25034_010",
  y50_59 = "B25034_009", y60_69 = "B25034_008", y70_79 = "B25034_007",
  med_inc = "B19013_001", pov_tot = "B17001_001", pov_below = "B17001_002",
  ten_tot = "B25003_001", renter = "B25003_003",
  pop_tot = "B02001_001", white = "B02001_002",
  age_tot = "B01001_001", m_u5 = "B01001_003", m_5 = "B01001_004",
  f_u5 = "B01001_027", f_5 = "B01001_028"
)

make_census <- function(acs_raw, areas_df, id_col) {
  acs_raw %>%
    transmute(
      geoid        = GEOID,
      pct_pre1980  = 100 * (pre39E + y40_49E + y50_59E + y60_69E + y70_79E) / pmax(tot_huE, 1),
      median_income = med_incE,
      pct_poverty  = 100 * pov_belowE / pmax(pov_totE, 1),
      pct_renter   = 100 * renterE / pmax(ten_totE, 1),
      pct_nonwhite = 100 * (pop_totE - whiteE) / pmax(pop_totE, 1),
      housing_units = tot_huE,
      pct_children_u6 = 100 * (m_u5E + m_5E + f_u5E + f_5E) / pmax(age_totE, 1)
    ) %>%
    left_join(areas_df, by = "geoid") %>%
    mutate(
      housing_density     = housing_units / pmax(area_sq_mi, 0.01),
      log_housing_density = log1p(housing_density)
    ) %>%
    select(-area_sq_mi) %>%
    rename(!!id_col := geoid)
}

bg_areas <- ne_bgs %>%
  mutate(area_sq_mi = as.numeric(st_area(geometry)) / 2589988) %>%
  st_drop_geometry() %>%
  select(geoid = GEOID, area_sq_mi)

bg_census <- get_acs(
  geography = "block group", state = "NE", year = 2022,
  variables = census_vars, output = "wide"
) %>%
  make_census(bg_areas, "bg_geoid")
saveRDS(bg_census, file.path(OUTPUT_PATH, "bg_census.rds"))

tract_areas <- ne_tracts %>%
  mutate(area_sq_mi = as.numeric(st_area(geometry)) / 2589988) %>%
  st_drop_geometry() %>%
  select(geoid = GEOID, area_sq_mi)

tract_census <- get_acs(
  geography = "tract", state = "NE", year = 2022,
  variables = census_vars, output = "wide"
) %>%
  make_census(tract_areas, "tract_geoid")
saveRDS(tract_census, file.path(OUTPUT_PATH, "tract_census.rds"))

# --- Verification ---

bg_census_check   <- readRDS(file.path(OUTPUT_PATH, "bg_census.rds"))
model_data_check  <- readRDS(file.path(OUTPUT_PATH, "model_data.rds"))
first_test_check  <- readRDS(file.path(OUTPUT_PATH, "first_test.rds"))

message("\n=== Verification ===")
message("Files saved: ", length(list.files(OUTPUT_PATH, "\\.rds$")))
message("model_data bg_geoid matches bg_census: ",
        sum(model_data_check$bg_geoid %in% bg_census_check$bg_geoid), "/", nrow(model_data_check))
message("first_test bg_geoid matches ne_bgs GEOID: ",
        sum(first_test_check$bg_geoid %in% ne_bgs$GEOID), "/", nrow(first_test_check))
message("first_test tract_geoid matches ne_tracts GEOID: ",
        sum(first_test_check$tract_geoid %in% ne_tracts$GEOID), "/", nrow(first_test_check))
message("Done.")
