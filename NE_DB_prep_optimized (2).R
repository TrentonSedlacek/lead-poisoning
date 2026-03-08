###############################################################################
# NE_DB_prep.R — Nebraska Childhood Lead Surveillance Dashboard
# Produces ALL .rds files consumed by NE_DB_app.R
# ACS modes: year-matched, non-overlapping fixed, endpoint (runtime)
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
         tract_num      = as.numeric(CensusTract2020),
         TRACTCE        = sprintf("%04d%02d", floor(tract_num), round((tract_num %% 1) * 100)),
         tract_geoid    = paste0("31", sprintf("%03d", as.integer(CensusCountyFips2020)), TRACTCE),
         bg_geoid       = paste0(tract_geoid, as.integer(CensusBlockGroup2020)),
         AddressIdentifier = gsub('^"+|"+$', '', AddressIdentifier)) |>
  filter(nchar(tract_geoid) == 11, nchar(bg_geoid) == 12) |>
  select(-tract_num, -TRACTCE) |>
  left_join(lhd_map, by = "county")

# --- 3. All Tests / First Test / Max Test ------------------------------------
save_rds(lead |> select(AddressIdentifier, PATIENT_LOCAL_ID, SAMPLE_DATE, result,
                        elevated, AGE_YR = Age_yr, AGE_MO, SEX = patient_SEX,
                        RACE_ETH = RACE, SAMPLE_TYPE = sample_type, sample_year, county, lhd),
         "all_tests")
first_test <- lead |> slice_min(sample_year, by = PATIENT_LOCAL_ID, with_ties = FALSE)
save_rds(first_test, "first_test")
max_test <- lead |> slice_max(result, by = PATIENT_LOCAL_ID, with_ties = FALSE) |>
  select(all_of(names(first_test)))
save_rds(max_test, "max_test")

# --- 4. Address Stats (with geometric mean pre-computed) ---------------------
calc_geo_mean <- function(x) {
  x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else exp(mean(log(pmax(x, 0.1))))
}

make_addr_stats <- function(ft) {
  ft |> summarize(
    n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
    mean_bll = mean(result, na.rm = TRUE),
    mean_bll_geo = calc_geo_mean(result),
    max_bll = max(result, na.rm = TRUE),
    first_bll = result[which.min(sample_year)], last_bll = result[which.max(sample_year)],
    first_year = min(sample_year, na.rm = TRUE), last_year = max(sample_year, na.rm = TRUE),
    lat = first(Latitude), lng = first(Longitude), county = first(county), lhd = first(lhd),
    bg_geoid = first(bg_geoid), tract_geoid = first(tract_geoid),
    .by = AddressIdentifier)
}

address_stats_all <- make_addr_stats(first_test)
address_stats <- filter(address_stats_all, n_children >= 2)
save_rds(address_stats_all, "address_stats_all"); save_rds(address_stats, "address_stats")

# --- 5. Pre-aggregated Geo Stats ---------------------------------------------
# Pre-compute per (geo_level x year x geo_id) for both test sources
# so the app can just filter instead of recomputing.
message("=== Pre-aggregating geo stats ===")

preaggregate_geo <- function(ft, geo_col) {
  child <- ft |>
    filter(!is.na(.data[[geo_col]])) |>
    summarize(
      n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
      mean_bll = mean(result, na.rm = TRUE),
      mean_bll_geo = calc_geo_mean(result),
      county = first(county),
      .by = c(sample_year, geo_id = .data[[geo_col]]))
  addr <- ft |>
    filter(!is.na(.data[[geo_col]])) |>
    summarize(e = max(elevated, na.rm = TRUE),
              .by = c(sample_year, geo_id = .data[[geo_col]], AddressIdentifier)) |>
    summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE),
              .by = c(sample_year, geo_id))
  left_join(child, addr, by = c("sample_year", "geo_id")) |>
    mutate(pct_elevated = 100 * n_elevated / n_children)
}

geo_levels <- c(lhd = "lhd", county = "county", tract = "tract_geoid", bg = "bg_geoid")
geo_stats_first <- list(); geo_stats_max <- list()
for (nm in names(geo_levels)) {
  geo_stats_first[[nm]] <- preaggregate_geo(first_test, geo_levels[nm]) |> mutate(geo_level = nm)
  geo_stats_max[[nm]]   <- preaggregate_geo(max_test,   geo_levels[nm]) |> mutate(geo_level = nm)
}
save_rds(bind_rows(geo_stats_first), "geo_stats_first")
save_rds(bind_rows(geo_stats_max),   "geo_stats_max")
message("  geo_stats_first: ", nrow(bind_rows(geo_stats_first)), " rows")
message("  geo_stats_max:   ", nrow(bind_rows(geo_stats_max)), " rows")

# --- 6. Model Data -----------------------------------------------------------
re_flag <- function(re, pattern, exclude_hisp = TRUE) {
  hit <- grepl(pattern, re) & !is.na(re)
  if (exclude_hisp) hit <- hit & !grepl("hispanic|latino", re)
  as.integer(hit)
}

add_race_flags_and_select <- function(sliced) {
  sliced |>
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
    select(AddressIdentifier, first_id = PATIENT_LOCAL_ID, first_year = sample_year,
           first_bll = result, first_age_mo = AGE_MO, first_sex = patient_SEX,
           first_sample_type = sample_type, bg_geoid, tract_geoid, county, lhd,
           lat = Latitude, lng = Longitude, starts_with("re_"))
}

build_addr_features <- function(ft, mode = "first") {
  sliced <- if (mode == "max") {
    ft |> slice_max(result, by = AddressIdentifier, with_ties = FALSE)
  } else {
    ft |> slice_min(sample_year, by = AddressIdentifier, with_ties = FALSE)
  }
  add_race_flags_and_select(sliced)
}

first_child   <- build_addr_features(first_test, "first")
max_bll_child <- build_addr_features(first_test, "max")

# 6a. Model data using first child as feature source
subsequent <- first_test |>
  inner_join(select(first_child, AddressIdentifier, first_id, first_year), by = "AddressIdentifier") |>
  filter(PATIENT_LOCAL_ID != first_id, sample_year >= first_year) |>
  summarize(n_subsequent = n(), n_sub_elevated = sum(elevated, na.rm = TRUE), .by = AddressIdentifier)

model_data <- first_child |>
  inner_join(subsequent, by = "AddressIdentifier") |>
  filter(!is.na(bg_geoid), nchar(bg_geoid) == 12) |>
  mutate(outcome = as.integer(n_sub_elevated > 0))
save_rds(model_data, "model_data")

# 6b. Model data using max-BLL child as feature source
subsequent_max <- first_test |>
  inner_join(select(max_bll_child, AddressIdentifier, first_id, first_year), by = "AddressIdentifier") |>
  filter(PATIENT_LOCAL_ID != first_id) |>
  summarize(n_subsequent = n(), n_sub_elevated = sum(elevated, na.rm = TRUE), .by = AddressIdentifier)

model_data_max <- max_bll_child |>
  inner_join(subsequent_max, by = "AddressIdentifier") |>
  filter(!is.na(bg_geoid), nchar(bg_geoid) == 12) |>
  mutate(outcome = as.integer(n_sub_elevated > 0))
save_rds(model_data_max, "model_data_max")

# 6c. All-address model data (outcome computed at runtime)
model_data_all     <- first_child   |> filter(!is.na(bg_geoid), nchar(bg_geoid) == 12)
model_data_all_max <- max_bll_child |> filter(!is.na(bg_geoid), nchar(bg_geoid) == 12)
save_rds(model_data_all, "model_data_all"); save_rds(model_data_all_max, "model_data_all_max")

# --- 7. Geographic Boundaries ------------------------------------------------
ne_counties <- counties(state = "NE", cb = TRUE, year = 2021)     |> st_transform(4326)
ne_tracts   <- tracts(state = "NE", cb = TRUE, year = 2021)       |> st_transform(4326)
ne_bgs      <- block_groups(state = "NE", cb = TRUE, year = 2021) |> st_transform(4326)
lhd_bounds  <- ne_counties |>
  left_join(lhd_map, by = c("NAME" = "county")) |>
  filter(!is.na(lhd)) |> summarize(geometry = st_union(geometry), .by = lhd)
walk(c("ne_counties","ne_tracts","ne_bgs","lhd_bounds"), \(nm) save_rds(get(nm), nm))

ne_tracts_2010 <- tracts(state = "NE", cb = TRUE, year = 2019) |> st_transform(4326)

# --- 8. Census (ACS 5-Year Estimates, All Vintages 2010-2023) ----------------
census_vars_stable <- c(
  tot_hu  = "B25034_001",
  med_inc = "B19013_001", pov_tot = "B17001_001", pov_below = "B17001_002",
  ten_tot = "B25003_001", renter  = "B25003_003",
  pop_tot = "B02001_001", white   = "B02001_002", age_tot  = "B01001_001",
  m_u5    = "B01001_003", m_5     = "B01001_004", f_u5    = "B01001_027", f_5 = "B01001_028",
  occ_tot = "B25002_001", vacant  = "B25002_003",
  rb_tot  = "B25070_001", rb_30   = "B25070_007", rb_35   = "B25070_008",
  rb_40   = "B25070_009", rb_50   = "B25070_010",
  med_hv  = "B25077_001")

b25034_new <- c(y70_79 = "B25034_007", y60_69 = "B25034_008",
                y50_59 = "B25034_009", y40_49 = "B25034_010", pre39 = "B25034_011")
b25034_old <- c(y70_79 = "B25034_006", y60_69 = "B25034_007",
                y50_59 = "B25034_008", y40_49 = "B25034_009", pre39 = "B25034_010")

get_census_vars <- function(yr) c(census_vars_stable, if (yr >= 2015) b25034_new else b25034_old)

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
      median_income   = med_incE, median_house_value = med_hvE,
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

# --- 8a. Tract crosswalk: 2010 GEOID <-> 2020 GEOID -------------------------
message("=== Building 2010-to-2020 tract crosswalk ===")
xwalk_rds_path <- file.path(OUTPUT_PATH, "tract_xwalk.rds")
xwalk_txt_path <- file.path(DATA_PATH, "tract_xwalk_2020_2010.txt")

if (file.exists(xwalk_rds_path)) {
  message("Loading cached crosswalk from .rds...")
  tract_xwalk <- readRDS(xwalk_rds_path)
} else {
  xwalk_url <- "https://www2.census.gov/geo/docs/maps-data/data/rel2020/tract/tab20_tract20_tract10_natl.txt"
  xwalk_ok <- FALSE
  if (file.exists(xwalk_txt_path) && file.size(xwalk_txt_path) < 15e6) {
    message("Existing crosswalk file truncated. Removing..."); file.remove(xwalk_txt_path)
  }
  if (!file.exists(xwalk_txt_path)) {
    message("Downloading tract relationship file (18 MB)...")
    old_dl <- getOption("timeout"); options(timeout = 600)
    tryCatch({
      download.file(xwalk_url, xwalk_txt_path, mode = "wb", quiet = FALSE)
      xwalk_ok <- file.exists(xwalk_txt_path) && file.size(xwalk_txt_path) > 15e6
      if (!xwalk_ok) { message("Download incomplete."); if (file.exists(xwalk_txt_path)) file.remove(xwalk_txt_path) }
    }, error = function(e) {
      message(sprintf("Download failed: %s", e$message))
      if (file.exists(xwalk_txt_path)) file.remove(xwalk_txt_path)
    })
    options(timeout = old_dl)
  } else xwalk_ok <- TRUE

  if (xwalk_ok) {
    message("Parsing Census relationship file...")
    tract_xwalk_raw <- read_delim(xwalk_txt_path, delim = "|", show_col_types = FALSE)
    needed <- c("STATE_2020","GEOID_TRACT_20","GEOID_TRACT_10","AREALAND_PART","AREALAND_TRACT_20")
    if (all(needed %in% names(tract_xwalk_raw))) {
      tract_xwalk <- tract_xwalk_raw |>
        filter(STATE_2020 == "31") |>
        transmute(geoid_2020 = GEOID_TRACT_20, geoid_2010 = GEOID_TRACT_10,
                  aland_pct = AREALAND_PART / pmax(AREALAND_TRACT_20, 1)) |>
        mutate(aland_pct = replace_na(aland_pct, 0)) |>
        group_by(geoid_2020) |> mutate(weight = aland_pct / pmax(sum(aland_pct), 1e-10)) |> ungroup()
    } else { message("Columns not as expected."); xwalk_ok <- FALSE }
  }
  if (!xwalk_ok) {
    message("Building crosswalk via spatial intersection...")
    t20 <- ne_tracts |> st_transform(32104); t10 <- ne_tracts_2010 |> st_transform(32104)
    xsf <- suppressWarnings(st_intersection(t20 |> select(geoid_2020=GEOID), t10 |> select(geoid_2010=GEOID)))
    xsf$inter_area <- as.numeric(st_area(xsf))
    t20a <- t20 |> mutate(total_area=as.numeric(st_area(geometry))) |> st_drop_geometry() |> select(geoid_2020=GEOID, total_area)
    tract_xwalk <- xsf |> st_drop_geometry() |> left_join(t20a, by="geoid_2020") |>
      mutate(aland_pct = inter_area/pmax(total_area,1)) |>
      group_by(geoid_2020) |> mutate(weight = aland_pct/pmax(sum(aland_pct),1e-10)) |> ungroup() |>
      select(geoid_2020, geoid_2010, aland_pct, weight) |> filter(weight > 0.001)
  }
}
message(sprintf("Crosswalk: %d rows, %d unique 2020 tracts, %d unique 2010 tracts",
                nrow(tract_xwalk), n_distinct(tract_xwalk$geoid_2020), n_distinct(tract_xwalk$geoid_2010)))
save_rds(tract_xwalk, "tract_xwalk")

# --- 8b. Pull ACS 2010-2023 and build time series ----------------------------
count_cols  <- c("tot_huE","pre39E","y40_49E","y50_59E","y60_69E","y70_79E",
                 "pov_totE","pov_belowE","ten_totE","renterE","pop_totE","whiteE","age_totE",
                 "m_u5E","m_5E","f_u5E","f_5E","occ_totE","vacantE",
                 "rb_totE","rb_30E","rb_35E","rb_40E","rb_50E")
median_cols <- c("med_incE", "med_hvE")
areas_tracts_2020 <- areas_from(ne_tracts)
areas_bgs_2020    <- areas_from(ne_bgs)

crosswalk_to_2020 <- function(acs_2010_raw, xwalk) {
  acs_data <- acs_2010_raw |> select(GEOID, all_of(c(count_cols, median_cols)))
  xwalk |> left_join(acs_data, by = c("geoid_2010" = "GEOID")) |>
    group_by(GEOID = geoid_2020) |>
    summarize(
      across(all_of(count_cols), \(x) sum(x * weight, na.rm = TRUE)),
      across(all_of(median_cols), \(x) {
        w <- weight[!is.na(x)]; v <- x[!is.na(x)]
        if (length(v) == 0) NA_real_ else sum(v * w) / pmax(sum(w), 1e-10)
      }), .groups = "drop")
}

old_timeout <- getOption("timeout"); options(timeout = 300)
tract_census_list <- list(); bg_census_list <- list()

for (yr in 2010:2023) {
  vars <- get_census_vars(yr)
  message(sprintf("--- Pulling ACS %d (tract) ---", yr))
  acs_tract <- tryCatch(
    get_acs(geography = "tract", state = "NE", year = yr, variables = vars, output = "wide"),
    error = function(e) { message(sprintf("  ERROR tract %d: %s", yr, e$message)); NULL })
  Sys.sleep(1)
  if (!is.null(acs_tract)) {
    if (yr <= 2019) {
      message(sprintf("  Crosswalking %d records...", nrow(acs_tract)))
      acs_tract_2020 <- crosswalk_to_2020(acs_tract, tract_xwalk)
      tract_df <- make_census(acs_tract_2020, areas_tracts_2020, "tract_geoid") |> mutate(acs_year = yr)
    } else {
      tract_df <- make_census(acs_tract, areas_tracts_2020, "tract_geoid") |> mutate(acs_year = yr)
    }
    tract_census_list[[as.character(yr)]] <- tract_df
    message(sprintf("  Tract %d: %d rows", yr, nrow(tract_df)))
  }
  if (yr >= 2020) {
    message(sprintf("--- Pulling ACS %d (block group) ---", yr))
    acs_bg <- tryCatch(
      get_acs(geography = "block group", state = "NE", year = yr, variables = vars, output = "wide"),
      error = function(e) { message(sprintf("  ERROR bg %d: %s", yr, e$message)); NULL })
    Sys.sleep(1)
    if (!is.null(acs_bg)) {
      bg_df <- make_census(acs_bg, areas_bgs_2020, "bg_geoid") |> mutate(acs_year = yr)
      bg_census_list[[as.character(yr)]] <- bg_df
      message(sprintf("  BG %d: %d rows", yr, nrow(bg_df)))
    }
  }
}
options(timeout = old_timeout)

tract_census_ts <- bind_rows(tract_census_list)
message(sprintf("\ntract_census_ts: %d rows (%d tracts x %d years)",
                nrow(tract_census_ts), n_distinct(tract_census_ts$tract_geoid),
                n_distinct(tract_census_ts$acs_year)))
save_rds(tract_census_ts, "tract_census_ts")

message("Building BG time series (pre-2020 from tract, 2020+ native)...")
bg_to_tract <- ne_bgs |> st_drop_geometry() |>
  transmute(bg_geoid = GEOID, tract_geoid = substr(GEOID, 1, 11))
bg_from_tract_list <- list()
for (yr in sort(unique(tract_census_ts$acs_year[tract_census_ts$acs_year <= 2019]))) {
  tract_yr <- tract_census_ts |> filter(acs_year == yr)
  if (nrow(tract_yr) == 0) next
  bg_from_tract_list[[as.character(yr)]] <- bg_to_tract |>
    left_join(tract_yr |> select(-any_of("acs_year")), by = "tract_geoid") |>
    mutate(acs_year = yr)
}
bg_census_ts <- bind_rows(bind_rows(bg_from_tract_list), bind_rows(bg_census_list)) |>
  arrange(bg_geoid, acs_year)
message(sprintf("bg_census_ts: %d rows (%d BGs x up to %d years)",
                nrow(bg_census_ts), n_distinct(bg_census_ts$bg_geoid), n_distinct(bg_census_ts$acs_year)))
save_rds(bg_census_ts, "bg_census_ts")

# Static snapshots (2023) for backward compat
bg_census    <- bg_census_ts    |> filter(acs_year == 2023) |> select(-acs_year)
tract_census <- tract_census_ts |> filter(acs_year == 2023) |> select(-acs_year)
save_rds(bg_census, "bg_census"); save_rds(tract_census, "tract_census")

# --- 9. Pre-join Census to Model Data ----------------------------------------
ACS_MIN <- min(tract_census_ts$acs_year); ACS_MAX <- max(tract_census_ts$acs_year)
acs_year_map <- function(yr) pmax(ACS_MIN, pmin(as.integer(yr), ACS_MAX))
NONOVERLAP_FIXED <- c(2013L, 2018L, 2023L)
assign_fixed_acs <- function(yr) {
  yr <- as.integer(yr)
  NONOVERLAP_FIXED[which.min(abs(yr - NONOVERLAP_FIXED))]
}

base_variants <- list(
  model_data = model_data, model_data_max = model_data_max,
  model_data_all = model_data_all, model_data_all_max = model_data_all_max)

geo_census <- list(
  tract = list(ts = tract_census_ts, col = "tract_geoid"),
  bg    = list(ts = bg_census_ts,    col = "bg_geoid"))

for (vname in names(base_variants)) {
  md <- base_variants[[vname]]
  for (gname in names(geo_census)) {
    gl <- geo_census[[gname]]
    # Year-matched
    joined_ym <- md |>
      mutate(acs_year = acs_year_map(first_year)) |>
      left_join(gl$ts, by = c(setNames("acs_year", "acs_year"), setNames(gl$col, gl$col))) |>
      select(-acs_year)
    save_rds(joined_ym, paste0(vname, "_", gname, "_ym"))
    # Non-overlapping fixed (2013, 2018, 2023)
    joined_nof <- md |>
      mutate(acs_year = sapply(first_year, assign_fixed_acs)) |>
      left_join(gl$ts, by = c(setNames("acs_year", "acs_year"), setNames(gl$col, gl$col))) |>
      select(-acs_year)
    save_rds(joined_nof, paste0(vname, "_", gname, "_nof"))
    message(sprintf("  Saved %s_%s_ym (%d) and _nof (%d)", vname, gname,
                    nrow(joined_ym), nrow(joined_nof)))
  }
}

# --- 10. ACS Under-6 Population Denominators ---------------------------------
county_lookup <- ne_counties |> st_drop_geometry() |> select(COUNTYFP, County = NAME)
bg_county <- ne_bgs |> st_drop_geometry() |>
  transmute(bg_geoid = GEOID, COUNTYFP = substr(GEOID, 3, 5)) |>
  left_join(county_lookup, by = "COUNTYFP")

tract_u6_ts  <- tract_census_ts |> select(geo_id = tract_geoid, acs_year, n_children_u6)
bg_u6_ts     <- bg_census_ts    |> select(geo_id = bg_geoid, acs_year, n_children_u6)
county_u6_ts <- bg_census_ts |> select(bg_geoid, acs_year, n_children_u6) |>
  left_join(bg_county, by = "bg_geoid") |>
  summarize(n_children_u6 = sum(n_children_u6, na.rm = TRUE), .by = c(County, acs_year)) |>
  filter(!is.na(County)) |> rename(geo_id = County)
lhd_u6_ts <- county_u6_ts |>
  left_join(lhd_map, by = c("geo_id" = "county")) |>
  summarize(n_children_u6 = sum(n_children_u6, na.rm = TRUE), .by = c(lhd, acs_year)) |>
  filter(!is.na(lhd)) |> rename(geo_id = lhd)
state_u6_ts <- bg_census_ts |>
  summarize(n_children_u6 = sum(n_children_u6, na.rm = TRUE), .by = acs_year)

acs_u6_ts <- list(bg = bg_u6_ts, tract = tract_u6_ts, county = county_u6_ts,
                  lhd = lhd_u6_ts, state = state_u6_ts)
save_rds(acs_u6_ts, "acs_u6_ts")

geo_county_lookup <- list(ne_tracts, ne_bgs) |>
  map(\(sf) sf |> st_drop_geometry() |>
        transmute(GEOID, COUNTYFP = substr(GEOID, 3, 5)) |>
        left_join(county_lookup, by = "COUNTYFP") |> select(GEOID, County)) |>
  list_rbind()
geo_lhd_lookup <- geo_county_lookup |>
  left_join(lhd_map, by = c("County" = "county")) |> select(GEOID, County, LHD = lhd)
save_rds(geo_county_lookup, "geo_county_lookup"); save_rds(geo_lhd_lookup, "geo_lhd_lookup")

# --- 11. Summaries -----------------------------------------------------------
save_rds(first_test |>
           summarize(total_children = n_distinct(PATIENT_LOCAL_ID),
                     sample_year = min(sample_year, na.rm = TRUE), .by = AddressIdentifier) |>
           mutate(type = if_else(total_children >= 2, "2+ Children", "1 Child")) |>
           summarize(n_addr = n(), .by = c(sample_year, type)), "addr_by_year")

total_acs_u6  <- sum(bg_census$n_children_u6, na.rm = TRUE)
n_elevated_u6 <- sum(first_test$elevated[first_test$Age_yr <= 5], na.rm = TRUE)

ssf <- function(d, grp = character(0), by_addr = FALSE) {
  if (by_addr) {
    d |> summarize(e = max(elevated, na.rm = TRUE), bll = mean(result, na.rm = TRUE),
                   .by = all_of(c(grp, "AddressIdentifier"))) |>
      summarize(mean_bll = mean(bll, na.rm = TRUE), n_children = n(),
                n_elevated = sum(e > 0, na.rm = TRUE),
                pct_elevated = 100 * mean(e > 0, na.rm = TRUE), .by = all_of(grp))
  } else {
    d |> summarize(mean_bll = mean(result, na.rm = TRUE), n_children = n(),
                   n_elevated = sum(elevated, na.rm = TRUE),
                   pct_elevated = 100 * mean(elevated, na.rm = TRUE), .by = all_of(grp))
  }
}
multi_child_addr <- first_test |>
  summarize(n = n(), .by = AddressIdentifier) |>
  filter(n >= 2) |> pull(AddressIdentifier)
save_rds(list(
  state_avgs = ssf(first_test, "sample_year"), state_avg_all = ssf(first_test),
  state_avgs_addr = ssf(first_test, "sample_year", by_addr = TRUE),
  state_avg_all_addr = ssf(first_test, by_addr = TRUE),
  addr_elevated_summary_all = summarize(first_test, e = max(elevated, na.rm = TRUE), .by = AddressIdentifier),
  addr_elevated_summary = first_test |>
    filter(AddressIdentifier %in% multi_child_addr) |>
    summarize(e = max(elevated, na.rm = TRUE), .by = AddressIdentifier),
  total_acs_u6 = total_acs_u6, n_elevated_u6 = n_elevated_u6,
  years = sort(unique(first_test$sample_year)),
  lhd_names = sort(unique(first_test$lhd[!is.na(first_test$lhd)]))), "startup_summaries")

# --- 12. Verification --------------------------------------------------------
expected <- c("first_test","max_test","all_tests","model_data","model_data_max",
              "model_data_all","model_data_all_max",
              "address_stats","address_stats_all",
              "ne_counties","ne_tracts","ne_bgs","lhd_bounds","lhd_map",
              "bg_census","tract_census",
              "geo_county_lookup","geo_lhd_lookup","addr_by_year","startup_summaries",
              "tract_xwalk","tract_census_ts","bg_census_ts",
              "geo_stats_first","geo_stats_max",
              paste0(rep(names(base_variants), each = 4), "_",
                     rep(c("tract_ym","bg_ym","tract_nof","bg_nof"), 4)),
              "acs_u6_ts")
present <- file.exists(file.path(OUTPUT_PATH, paste0(expected, ".rds")))
message(sprintf("\n=== Verification: %d/%d files present ===", sum(present), length(expected)))
if (!all(present)) message("Missing: ", paste(expected[!present], collapse = ", "))
message(sprintf("first_test: %d | model_data: %d | tract_census_ts: %d | bg_census_ts: %d",
                nrow(first_test), nrow(model_data), nrow(tract_census_ts), nrow(bg_census_ts)))
message(sprintf("geo_stats_first: %d | geo_stats_max: %d",
                nrow(bind_rows(geo_stats_first)), nrow(bind_rows(geo_stats_max))))
message("Done.")
