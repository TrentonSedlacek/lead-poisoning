###############################################################################
# NE_DB_prep.R — Nebraska Childhood Lead Surveillance Dashboard
# Produces ALL .rds files consumed by NE_DB_app_OPTIMIZED.R
###############################################################################
library(tidyverse); library(sf); library(tigris); library(tidycensus)
options(tigris_use_cache = TRUE)

DATA_PATH   <- "K:/CLPPP/TrentonS"
OUTPUT_PATH <- file.path(DATA_PATH, "shiny_data")
BLL_CUTOFF  <- 3.5
dir.create(OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)
save_rds <- function(obj, name) saveRDS(obj, file.path(OUTPUT_PATH, paste0(name, ".rds")))

# --- 1. LHD-County Mapping ---------------------------------------------------
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
  "Wheeler", "Loup Basin", "York", "Four Corners")
save_rds(lhd_map, "lhd_map")

# --- 2. Load and Clean Lead Data ---------------------------------------------
lead <- read_csv(file.path(DATA_PATH, "MatchedGeo2010-2023.csv"),
                 col_types = cols(.default = col_guess()), show_col_types = FALSE) |>
  filter(Age_yr <= 6, !is.na(Latitude), Latitude != 0, !is.na(Longitude),
         FeatureMatchingResultType != "Unmatchable", !is.na(result),
         !is.na(PATIENT_LOCAL_ID),
         !is.na(CensusTract2020), as.numeric(CensusTract2020) > 0,
         !is.na(CensusCountyFips2020), !is.na(CensusBlockGroup2020)) |>
  mutate(elevated       = as.integer(result >= BLL_CUTOFF),
         county         = str_to_title(str_trim(county)),
         county         = if_else(county == "Mcpherson", "McPherson", county),
         street_address = paste(Street, City, State, Zip, sep = ", "),
         tract_num      = as.numeric(CensusTract2020),
         TRACTCE        = sprintf("%04d%02d", floor(tract_num), round((tract_num %% 1) * 100)),
         tract_geoid    = paste0("31", sprintf("%03d", as.integer(CensusCountyFips2020)), TRACTCE),
         bg_geoid       = paste0(tract_geoid, as.integer(CensusBlockGroup2020))) |>
  filter(nchar(tract_geoid) == 11, nchar(bg_geoid) == 12) |>
  select(-tract_num, -TRACTCE) |>
  left_join(lhd_map, by = "county")

# --- 3. All Tests / First Test / Max Test -------------------------------------
save_rds(lead |> select(address_id, PATIENT_LOCAL_ID, street_address, SAMPLE_DATE, result,
                        elevated, AGE_YR = Age_yr, AGE_MO, SEX = patient_SEX,
                        RACE_ETH = RACE, SAMPLE_TYPE = sample_type, sample_year, county, lhd),
         "all_tests")
first_test <- lead |> slice_min(sample_year, by = PATIENT_LOCAL_ID, with_ties = FALSE)
save_rds(first_test, "first_test")
# Priority-based max: prefer V (venous), then C (capillary), then unknown
max_test <- lead |>
  mutate(.st_priority = case_when(
    toupper(sample_type) == "V" ~ 1L,
    toupper(sample_type) == "C" ~ 2L,
    TRUE ~ 3L)) |>
  slice_min(.st_priority, by = PATIENT_LOCAL_ID, with_ties = TRUE) |>
  slice_max(result, by = PATIENT_LOCAL_ID, with_ties = FALSE) |>
  select(all_of(names(first_test)))
save_rds(max_test, "max_test")

# --- 4. Address Stats --------------------------------------------------------
address_stats_all <- first_test |>
  summarize(n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
            mean_bll = mean(result, na.rm = TRUE), max_bll = max(result, na.rm = TRUE),
            first_bll = result[which.min(sample_year)], last_bll = result[which.max(sample_year)],
            first_year = min(sample_year, na.rm = TRUE), last_year = max(sample_year, na.rm = TRUE),
            lat = first(Latitude), lng = first(Longitude), county = first(county), lhd = first(lhd),
            bg_geoid = first(bg_geoid), tract_geoid = first(tract_geoid),
            street_address = first(street_address), .by = address_id)
address_stats <- filter(address_stats_all, n_children >= 2)
save_rds(address_stats_all, "address_stats_all"); save_rds(address_stats, "address_stats")

# --- 5. Model Data (first child at each 2+ child address + outcome) ----------
re_flag <- function(re, pattern, exclude_hisp = TRUE) {
  hit <- grepl(pattern, re) & !is.na(re)
  if (exclude_hisp) hit <- hit & !grepl("hispanic|latino", re)
  as.integer(hit)
}
first_child <- first_test |>
  slice_min(sample_year, by = address_id, with_ties = FALSE) |>
  mutate(re          = tolower(trimws(RACE)),
         re_hispanic = re_flag(re, "hispanic|latino", exclude_hisp = FALSE),
         re_white_nh = re_flag(re, "white"),
         re_black_nh = re_flag(re, "black|african"),
         re_aian_nh  = re_flag(re, "indian|alaska|aian|native am"),
         re_asian_nh = re_flag(re, "asian"),
         re_nhopi_nh = re_flag(re, "pacific|hawaiian|nhopi"),
         re_multi_nh = re_flag(re, "multi|two|more"),
         re_other_nh = re_flag(re, "other"),
         across(starts_with("re_"), \(x) if_else(is.na(RACE), NA_integer_, x))) |>
  select(address_id, first_id = PATIENT_LOCAL_ID, first_year = sample_year,
         first_bll = result, first_age_mo = AGE_MO, first_sex = patient_SEX,
         first_sample_type = sample_type, bg_geoid, tract_geoid, county, lhd,
         lat = Latitude, lng = Longitude, street_address, starts_with("re_"))

subsequent <- first_test |>
  inner_join(select(first_child, address_id, first_id, first_year), by = "address_id") |>
  filter(PATIENT_LOCAL_ID != first_id, sample_year >= first_year) |>
  summarize(n_subsequent = n(), n_sub_elevated = sum(elevated, na.rm = TRUE), .by = address_id)

model_data <- first_child |>
  inner_join(subsequent, by = "address_id") |>
  filter(!is.na(bg_geoid), nchar(bg_geoid) == 12) |>
  mutate(outcome = as.integer(n_sub_elevated > 0))
save_rds(model_data, "model_data")

# --- 6. Geographic Boundaries -------------------------------------------------
ne_counties <- counties(state = "NE", cb = TRUE, year = 2021)     |> st_transform(4326)
ne_tracts   <- tracts(state = "NE", cb = TRUE, year = 2021)       |> st_transform(4326)
ne_bgs      <- block_groups(state = "NE", cb = TRUE, year = 2021) |> st_transform(4326)
lhd_bounds  <- ne_counties |>
  left_join(lhd_map, by = c("NAME" = "county")) |>
  filter(!is.na(lhd)) |> summarize(geometry = st_union(geometry), .by = lhd)
walk(c("ne_counties","ne_tracts","ne_bgs","lhd_bounds"), \(nm) save_rds(get(nm), nm))

# --- 7. Census (ACS 2019-2023) -----------------------------------------------
census_vars <- c(
  tot_hu = "B25034_001", pre39 = "B25034_011", y40_49 = "B25034_010",
  y50_59 = "B25034_009", y60_69 = "B25034_008", y70_79 = "B25034_007",
  med_inc = "B19013_001", pov_tot = "B17001_001", pov_below = "B17001_002",
  ten_tot = "B25003_001", renter = "B25003_003",
  pop_tot = "B02001_001", white = "B02001_002", age_tot = "B01001_001",
  m_u5 = "B01001_003", m_5 = "B01001_004", f_u5 = "B01001_027", f_5 = "B01001_028",
  occ_tot = "B25002_001", vacant = "B25002_003",
  rb_tot = "B25070_001", rb_30 = "B25070_007", rb_35 = "B25070_008",
  rb_40 = "B25070_009", rb_50 = "B25070_010")

areas_from <- function(sf_obj) {
  sf_obj |> mutate(area_sq_mi = as.numeric(st_area(geometry)) / 2589988) |>
    st_drop_geometry() |> select(geoid = GEOID, area_sq_mi)
}
make_census <- function(acs_raw, areas_df, id_col) {
  acs_raw |>
    transmute(
      geoid           = GEOID,
      pct_pre1980     = 100 * (pre39E + y40_49E + y50_59E + y60_69E + y70_79E) / pmax(tot_huE, 1),
      pct_pre1950     = 100 * (pre39E + y40_49E) / pmax(tot_huE, 1),
      median_income   = med_incE,
      pct_poverty     = if_else(pov_totE >= 10, 100 * pov_belowE / pov_totE, NA_real_),
      pct_renter      = 100 * renterE / pmax(ten_totE, 1),
      pct_nonwhite    = 100 * (pop_totE - whiteE) / pmax(pop_totE, 1),
      housing_units   = tot_huE,
      pct_children_u6 = 100 * (m_u5E + m_5E + f_u5E + f_5E) / pmax(age_totE, 1),
      n_children_u6   = m_u5E + m_5E + f_u5E + f_5E,
      pct_vacant      = if_else(occ_totE >= 5, 100 * vacantE / occ_totE, NA_real_),
      pct_rent_burden = if_else(rb_totE >= 10,
                                100 * (rb_30E + rb_35E + rb_40E + rb_50E) / rb_totE, NA_real_),
      race_cat = factor(case_when(pct_nonwhite >= 50 ~ "Majority Non-White",
                                  pct_nonwhite >= 25 ~ "Mixed", TRUE ~ "Majority White"),
                        levels = c("Majority White", "Mixed", "Majority Non-White")),
      income_cat = factor(case_when(is.na(med_incE) ~ NA_character_, med_incE < 35000 ~ "Low (<$35k)",
                                    med_incE < 55000 ~ "Moderate ($35-55k)",
                                    med_incE < 80000 ~ "Middle ($55-80k)", TRUE ~ "High (>$80k)"),
                          levels = c("Low (<$35k)","Moderate ($35-55k)","Middle ($55-80k)","High (>$80k)"))
    ) |>
    left_join(areas_df, by = "geoid") |>
    mutate(housing_density = housing_units / pmax(area_sq_mi, 0.01),
           log_housing_density = log1p(housing_density)) |>
    select(-area_sq_mi) |> rename(!!id_col := geoid)
}

bg_census <- get_acs(geography = "block group", state = "NE", year = 2023,
                     variables = census_vars, output = "wide") |>
  make_census(areas_from(ne_bgs), "bg_geoid")
tract_census <- get_acs(geography = "tract", state = "NE", year = 2023,
                        variables = census_vars, output = "wide") |>
  make_census(areas_from(ne_tracts), "tract_geoid")
save_rds(bg_census, "bg_census"); save_rds(tract_census, "tract_census")

# --- 8. Pre-join Census -------------------------------------------------------
save_rds(left_join(model_data, bg_census, by = "bg_geoid"),          "model_data_bg")
save_rds(left_join(model_data, tract_census, by = "tract_geoid"),    "model_data_tract")
save_rds(left_join(address_stats, bg_census, by = "bg_geoid"),       "address_stats_bg")
save_rds(left_join(address_stats, tract_census, by = "tract_geoid"), "address_stats_tract")

# --- 9. ACS Under-6 + GEOID Lookups ------------------------------------------
county_lookup <- ne_counties |> st_drop_geometry() |> select(COUNTYFP, County = NAME)
bg_county <- ne_bgs |> st_drop_geometry() |>
  transmute(bg_geoid = GEOID, COUNTYFP = substr(GEOID, 3, 5)) |>
  left_join(county_lookup, by = "COUNTYFP")

county_u6 <- bg_census |> select(bg_geoid, n_children_u6) |>
  left_join(bg_county, by = "bg_geoid") |>
  summarize(n_acs_u6 = sum(n_children_u6, na.rm = TRUE), .by = County) |>
  filter(!is.na(County)) |> rename(geo_id = County)

acs_u6 <- list(
  bg    = select(bg_census, geo_id = bg_geoid, n_acs_u6 = n_children_u6),
  tract = select(tract_census, geo_id = tract_geoid, n_acs_u6 = n_children_u6),
  county = county_u6,
  lhd = county_u6 |> left_join(lhd_map, by = c("geo_id" = "county")) |>
    summarize(n_acs_u6 = sum(n_acs_u6, na.rm = TRUE), .by = lhd) |>
    filter(!is.na(lhd)) |> rename(geo_id = lhd))
total_acs_u6  <- sum(bg_census$n_children_u6, na.rm = TRUE)
n_elevated_u6 <- sum(first_test$elevated[first_test$Age_yr <= 5], na.rm = TRUE)
save_rds(acs_u6, "acs_u6"); save_rds(total_acs_u6, "total_acs_u6"); save_rds(n_elevated_u6, "n_elevated_u6")

geo_county_lookup <- list(ne_tracts, ne_bgs) |>
  map(\(sf) sf |> st_drop_geometry() |>
        transmute(GEOID, COUNTYFP = substr(GEOID, 3, 5)) |>
        left_join(county_lookup, by = "COUNTYFP") |> select(GEOID, County)) |>
  list_rbind()
geo_lhd_lookup <- geo_county_lookup |>
  left_join(lhd_map, by = c("County" = "county")) |> select(GEOID, County, LHD = lhd)
save_rds(geo_county_lookup, "geo_county_lookup"); save_rds(geo_lhd_lookup, "geo_lhd_lookup")

# --- 10. Address Counts by Year + Startup Summaries ---------------------------
save_rds(first_test |>
           summarize(total_children = n_distinct(PATIENT_LOCAL_ID),
                     sample_year = min(sample_year, na.rm = TRUE), .by = address_id) |>
           mutate(type = if_else(total_children >= 2, "2+ Children", "1 Child")) |>
           summarize(n_addr = n(), .by = c(sample_year, type)), "addr_by_year")

ssf <- function(d, grp = character(0), by_addr = FALSE) {
  if (by_addr) {
    d |> summarize(e = max(elevated, na.rm = TRUE), bll = mean(result, na.rm = TRUE),
                   .by = all_of(c(grp, "address_id"))) |>
      summarize(mean_bll = mean(bll, na.rm = TRUE), n_children = n(),
                n_elevated = sum(e > 0, na.rm = TRUE),
                pct_elevated = 100 * mean(e > 0, na.rm = TRUE), .by = all_of(grp))
  } else {
    d |> summarize(mean_bll = mean(result, na.rm = TRUE), n_children = n(),
                   n_elevated = sum(elevated, na.rm = TRUE),
                   pct_elevated = 100 * mean(elevated, na.rm = TRUE), .by = all_of(grp))
  }
}
save_rds(list(
  state_avgs = ssf(first_test, "sample_year"), state_avg_all = ssf(first_test),
  state_avgs_addr = ssf(first_test, "sample_year", by_addr = TRUE),
  state_avg_all_addr = ssf(first_test, by_addr = TRUE),
  addr_elevated_summary = summarize(first_test, e = max(elevated, na.rm = TRUE), .by = address_id),
  total_acs_u6 = total_acs_u6, n_elevated_u6 = n_elevated_u6,
  years = sort(unique(first_test$sample_year)),
  lhd_names = sort(unique(first_test$lhd[!is.na(first_test$lhd)]))), "startup_summaries")

# --- 11. Verification ---------------------------------------------------------
expected <- c("first_test","max_test","all_tests","model_data","address_stats",
              "address_stats_all","ne_counties","ne_tracts","ne_bgs","lhd_bounds","lhd_map",
              "bg_census","tract_census","model_data_bg","model_data_tract",
              "address_stats_bg","address_stats_tract","acs_u6","total_acs_u6",
              "n_elevated_u6","geo_county_lookup","geo_lhd_lookup","addr_by_year","startup_summaries")
message(sprintf("\n=== Verification: %d files, all present: %s ===",
                length(list.files(OUTPUT_PATH, "\\.rds$")),
                all(file.exists(file.path(OUTPUT_PATH, paste0(expected, ".rds"))))))
message(sprintf("first_test: %d | max_test: %d | model_data: %d | total_acs_u6: %d",
                nrow(first_test), nrow(max_test), nrow(model_data), total_acs_u6))
message("Done — all files needed by NE_DB_app_OPTIMIZED.R")