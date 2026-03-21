###############################################################################
# NE_DB_app.R — Nebraska Childhood Lead Surveillance Dashboard
# Requires: NE_DB_prep.R output in K:/CLPPP/TrentonS/shiny_data
###############################################################################

# Steps
# 1. CTRL + A
# 2. CTRL + Enter (or click Run)
# 3. Use app. Wait ~60 seconds upon start for default risk model to run
# If maps aren't loading or appear stuck, click a different setting to refresh

packages <- c("shiny", "shinydashboard", "leaflet", "sf", "tidyverse",
              "DT", "lme4", "broom.mixed", "plotly", "xgboost", "scales")
install.packages(setdiff(packages, rownames(installed.packages())))
invisible(lapply(packages, library, character.only = TRUE))

'%||%' <- function(a, b) if (!is.null(a)) a else b

DATA_PATH  <- "K:/CLPPP/TrentonS/shiny_data"
BLL_CUTOFF <- 3.5

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

calc_mean_bll <- function(x, use_geo = FALSE) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  if (use_geo) exp(mean(log(pmax(x, 0.1)))) else mean(x)
}

bll_mean_label <- function(use_geo) if (use_geo) "(Geometric)" else "(Arithmetic)"

summarize_lead <- function(data, geo_col, sup = 0, use_geo = FALSE) {
  data <- data %>% filter(!is.na(.data[[geo_col]]))
  child <- data %>%
    group_by(geo_id = .data[[geo_col]]) %>%
    summarize(n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
              mean_bll = calc_mean_bll(result, use_geo), county = first(county), .groups = "drop")
  addr <- data %>%
    group_by(geo_id = .data[[geo_col]], address_id) %>%
    summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
    group_by(geo_id) %>%
    summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")
  left_join(child, addr, by = "geo_id") %>%
    mutate(pct_elevated = 100 * n_elevated / n_children, suppressed = n_children < sup)
}

state_summary <- function(data, by_year = FALSE, by_address = FALSE, use_geo = FALSE) {
  if (by_address) {
    grp <- if (by_year) group_by(data, sample_year, address_id) else group_by(data, address_id)
    addr <- grp %>%
      summarize(e = max(elevated, na.rm = TRUE),
                bll = calc_mean_bll(result, use_geo), .groups = "drop")
    if (by_year) addr <- addr %>% group_by(sample_year)
    addr %>% summarize(mean_bll = mean(bll, na.rm = TRUE), n_children = n(),
                       n_elevated = sum(e > 0, na.rm = TRUE),
                       pct_elevated = 100 * mean(e > 0, na.rm = TRUE), .groups = "drop")
  } else {
    if (by_year) data <- data %>% group_by(sample_year)
    data %>% summarize(mean_bll = calc_mean_bll(result, use_geo), n_children = n(),
                       n_elevated = sum(elevated, na.rm = TRUE),
                       pct_elevated = 100 * mean(elevated, na.rm = TRUE), .groups = "drop")
  }
}

format_table <- function(df, id_name, unit = "Children", show_county = FALSE,
                         show_lhd = FALSE, mean_label = "") {
  tested   <- if (unit == "Addresses") df$n_addresses  else df$n_children
  elevated <- if (unit == "Addresses") df$addr_elevated else df$n_elevated
  bll_col_name <- if (nchar(mean_label) > 0) paste("Average BLL", mean_label) else "Average BLL"
  out <- tibble(!!id_name := df$geo_id, Tested = tested, Elevated = elevated,
                `% Elevated (elevated/tested)` = round(100 * elevated / pmax(tested, 1), 2),
                !!bll_col_name := round(df$mean_bll, 2), Suppressed = df$suppressed)
  if (show_lhd && "lhd" %in% names(df))
    out <- bind_cols(tibble(LHD = df$lhd), out)
  if (show_county) out <- bind_cols(tibble(County = df$county), out)
  out
}

compute_roc <- function(labels, scores) {
  ord <- order(scores, decreasing = TRUE)
  labels <- labels[ord]; scores <- scores[ord]
  n_pos <- sum(labels == 1); n_neg <- sum(labels == 0)
  tp <- cumsum(labels == 1); fp <- cumsum(labels == 0)
  sens <- c(0, tp / max(n_pos, 1)); spec <- c(1, 1 - fp / max(n_neg, 1))
  fpr <- 1 - spec
  auc <- abs(sum(diff(fpr) * (sens[-1] + sens[-length(sens)]) / 2))
  youden_idx <- which.max(sens + spec - 1)
  list(sensitivities = sens, specificities = spec, thresholds = c(Inf, scores), auc = auc,
       youden = list(sens = sens[youden_idx], spec = spec[youden_idx]))
}

build_quintile_palette <- function(vals, colors = c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026")) {
  q <- quantile(vals, probs = seq(0, 1, 0.2), na.rm = TRUE)
  list(breaks = q, colors = colors,
       labels = paste0("Q", 1:5, ": ", round(q[1:5], 1), "-", round(q[2:6], 1)))
}

assign_fill_colors <- function(geo, display_var, color_func) {
  v <- geo[[display_var]]
  sup <- if (!is.null(geo$suppressed)) !is.na(geo$suppressed) & geo$suppressed else rep(FALSE, length(v))
  ifelse(sup, "#666666", ifelse(is.na(v) | !is.finite(v), "#cccccc", sapply(v, color_func)))
}

is_year_var <- function(v) grepl("year", v, ignore.case = TRUE)

VAR_LABELS <- c(
  pct_pre1980 = "Pre-1980 Housing (%)", pct_pre1950 = "Pre-1950 Housing (%)",
  pct_poverty = "Poverty Rate (%)", pct_rent_burden = "Rent Burden \u226530% (%)",
  pct_renter = "Renter-Occupied (%)",
  median_income = "Median Income ($)", income_cat = "Income Category",
  pct_nonwhite = "Non-White (%)", race_cat = "Race Category",
  log_housing_density = "Housing Density (log)",
  pct_vacant = "Vacant Housing (%)",
  first_bll = "First Child's BLL", first_year = "First Test Year",
  first_sex = "First Child Sex", first_age_mo = "First Child Age (months)",
  first_sample_type = "Sample Type (Capillary/Venous)",
  re_hispanic = "Hispanic", re_white_nh = "White (Non-Hispanic)",
  re_black_nh = "Black (Non-Hispanic)", re_aian_nh = "AIAN (Non-Hispanic)",
  re_asian_nh = "Asian (Non-Hispanic)", re_nhopi_nh = "NHOPI (Non-Hispanic)",
  re_multi_nh = "Two or More Races (Non-Hispanic)", re_other_nh = "Other Race (Non-Hispanic)")

METRIC_LABELS <- c(
  n_elevated = "Elevated Count", n_children = "Tested Count",
  pct_elevated = "% Elevated (elevated/tested)", mean_bll = "Average BLL (\u00b5g/dL)",
  n_addresses = "Addresses Tested", addr_elevated = "Addresses Elevated")

get_label <- function(v) {
  if (v %in% names(VAR_LABELS))    return(unname(VAR_LABELS[v]))
  if (v %in% names(METRIC_LABELS)) return(unname(METRIC_LABELS[v]))
  gsub("_", " ", tools::toTitleCase(v))
}

get_bll_color <- function(bll) {
  case_when(bll >= 20 ~ "#67000d", bll >= 10 ~ "#a50f15", bll >= 5 ~ "#d73027",
            bll >= 3.5 ~ "#fc8d59", TRUE ~ "#c9a84c")
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

load_rds <- function(name) {
  path <- file.path(DATA_PATH, paste0(name, ".rds"))
  if (file.exists(path)) readRDS(path) else NULL
}

first_test    <- load_rds("first_test")
model_data    <- load_rds("model_data")
address_stats <- load_rds("address_stats")
ne_counties   <- load_rds("ne_counties")
lhd_bounds    <- load_rds("lhd_bounds")
ne_bgs        <- load_rds("ne_bgs")
ne_tracts     <- load_rds("ne_tracts")
bg_census     <- load_rds("bg_census")
tract_census  <- load_rds("tract_census")
lhd_map_df    <- load_rds("lhd_map")
all_tests     <- load_rds("all_tests")

model_data_bg    <- load_rds("model_data_bg")
model_data_tract <- load_rds("model_data_tract")
addr_stats_bg    <- load_rds("address_stats_bg")
addr_stats_tract <- load_rds("address_stats_tract")

if (!is.null(first_test) && "county" %in% names(first_test))
  first_test <- first_test %>% mutate(county = if_else(county == "Mcpherson", "McPherson", county))

max_test <- load_rds("max_test")
if (is.null(max_test) && !is.null(all_tests)) {
  max_test <- all_tests %>%
    arrange(PATIENT_LOCAL_ID, desc(result), sample_year) %>%
    group_by(PATIENT_LOCAL_ID) %>% dplyr::slice(1) %>% ungroup()
  shared_cols <- intersect(names(first_test), names(max_test))
  max_test <- max_test %>% select(all_of(shared_cols))
}
if (!is.null(max_test) && "county" %in% names(max_test))
  max_test <- max_test %>% mutate(county = if_else(county == "Mcpherson", "McPherson", county))

address_stats_all <- load_rds("address_stats_all")
if (is.null(address_stats_all)) {
  address_stats_all <- first_test %>%
    group_by(address_id) %>%
    summarize(
      n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
      mean_bll = mean(result, na.rm = TRUE), max_bll = max(result, na.rm = TRUE),
      first_bll = result[which.min(sample_year)][1],
      first_year = min(sample_year, na.rm = TRUE), last_year = max(sample_year, na.rm = TRUE),
      street_address = first(street_address), county = first(county), lhd = first(lhd),
      lat = first(Latitude), lng = first(Longitude),
      tract_geoid = first(tract_geoid), bg_geoid = first(bg_geoid), .groups = "drop")
}

sums <- load_rds("startup_summaries")
if (!is.null(sums)) {
  state_avgs            <- sums$state_avgs
  state_avg_all         <- sums$state_avg_all
  state_avgs_addr       <- sums$state_avgs_addr
  state_avg_all_addr    <- sums$state_avg_all_addr
  addr_elevated_summary <- sums$addr_elevated_summary
  addr_elevated_summary_all <- sums$addr_elevated_summary_all %||% sums$addr_elevated_summary
  years     <- sums$years     %||% (if (!is.null(first_test)) sort(unique(first_test$sample_year)) else 2010:2023)
  lhd_names <- sums$lhd_names %||% character(0)
} else {
  years     <- if (!is.null(first_test)) sort(unique(first_test$sample_year)) else 2010:2023
  lhd_names <- if (!is.null(first_test) && "lhd" %in% names(first_test))
    sort(unique(first_test$lhd[!is.na(first_test$lhd)])) else character(0)
  state_avgs         <- state_summary(first_test, by_year = TRUE)
  state_avg_all      <- state_summary(first_test)
  state_avgs_addr    <- state_summary(first_test, by_year = TRUE, by_address = TRUE)
  state_avg_all_addr <- state_summary(first_test, by_address = TRUE)
  addr_elevated_summary <- first_test %>%
    group_by(address_id) %>% summarize(e = max(elevated, na.rm = TRUE), .groups = "drop")
  addr_elevated_summary_all <- addr_elevated_summary
}

addr_by_year <- load_rds("addr_by_year")
if (is.null(addr_by_year)) {
  addr_child_counts <- first_test %>%
    group_by(address_id) %>%
    summarize(total_children = n_distinct(PATIENT_LOCAL_ID),
              first_year = min(sample_year, na.rm = TRUE), .groups = "drop")
  addr_by_year <- addr_child_counts %>%
    mutate(type = if_else(total_children >= 2, "2+ Children", "1 Child")) %>%
    group_by(first_year, type) %>% summarize(n_addr = n(), .groups = "drop") %>%
    rename(sample_year = first_year)
}

geo_county_lookup <- load_rds("geo_county_lookup")
geo_lhd_lookup    <- load_rds("geo_lhd_lookup")
if (is.null(geo_county_lookup) && !is.null(ne_tracts) && !is.null(ne_counties)) {
  tract_county <- ne_tracts %>% st_drop_geometry() %>%
    mutate(COUNTYFP = substr(GEOID, 3, 5)) %>%
    left_join(ne_counties %>% st_drop_geometry() %>% select(COUNTYFP, NAME), by = "COUNTYFP") %>%
    select(GEOID, County = NAME)
  bg_county <- if (!is.null(ne_bgs)) {
    ne_bgs %>% st_drop_geometry() %>%
      mutate(COUNTYFP = substr(GEOID, 3, 5)) %>%
      left_join(ne_counties %>% st_drop_geometry() %>% select(COUNTYFP, NAME), by = "COUNTYFP") %>%
      select(GEOID, County = NAME)
  } else NULL
  geo_county_lookup <- bind_rows(tract_county, bg_county)
  if (!is.null(lhd_map_df))
    geo_lhd_lookup <- geo_county_lookup %>%
    left_join(lhd_map_df, by = c("County" = "county")) %>%
    select(GEOID, County, LHD = lhd)
}

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "NE Lead Dashboard"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview",           tabName = "overview", icon = icon("home")),
      menuItem("Surveillance Map",   tabName = "map",      icon = icon("map")),
      menuItem("Priority Addresses", tabName = "ranked",   icon = icon("list-ol")),
      menuItem("Census Data",        tabName = "census",   icon = icon("chart-bar")),
      menuItem("Data Tables",        tabName = "tables",   icon = icon("table")),
      menuItem("Risk Model",         tabName = "model",    icon = icon("chart-line"))
    ),
    hr(),
    radioButtons("test_selection", "Child Test Selection:",
                 c("First Test per Child" = "first", "Maximum BLL per Child" = "max"), inline = FALSE),
    p(em("First = earliest test. Maximum = highest BLL recorded."),
      style = "padding: 0 15px; font-size: 10px; color: #aaa;"),
    hr(),
    numericInput("suppress_n", "Suppress areas with fewer than:", 10, min = 0, max = 50, step = 5),
    radioButtons("use_geo_mean", "BLL Average Type:",
                 c("Arithmetic Mean" = "arith", "Geometric Mean" = "geo"), inline = TRUE),
    radioButtons("addr_pop", "Address Population:",
                 c("All Tested (1+ child)" = "all", "2+ Children Tested" = "multi"), inline = TRUE),
    p(em("Affects address counts. Switch 'Count By' to 'Addresses' on maps.
         Also applies to Priority Addresses and Data Tables."),
      style = "padding: 0 15px; font-size: 10px; color: #aaa;")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .section-title { font-size: 20px; font-weight: 600; margin-bottom: 6px; }
      .section-desc { color: #444; margin-bottom: 15px; font-size: 13.5px; line-height: 1.6; }
      .metric-box { background: #f8f9fa; border-left: 3px solid #2c7fb8;
                    padding: 14px; margin: 8px 0; border-radius: 4px; }
      .state-avg-box { background: #e8f4f8; border: 1px solid #2c7fb8; padding: 10px;
                       border-radius: 4px; margin-top: 10px; font-size: 13px; }
      .compare-box { background: #f0f7f0; border: 1px solid #4a9; padding: 15px;
                     border-radius: 6px; margin: 10px 0; }
      .warning-box { background: #fff3cd; border: 1px solid #ffc107; padding: 10px;
                     border-radius: 4px; margin: 5px 0; font-size: 12px; }
      .toggle-group { border: 1px solid #ddd; border-radius: 4px; padding: 6px 8px;
                      margin-bottom: 6px; background: #fafafa; }
      .cor-high { background-color: #f8d7da !important; }
      .cor-borderline { background-color: #fff3cd !important; }
      .content-wrapper { font-size: 13.5px; }
      .box-title { font-size: 15px; }
      .tab-content .row > .col-sm-3 { width: 29.5% !important; }
      .tab-content .row > .col-sm-9 { width: 70.5% !important; }
    "))),

    tabItems(

      # ===== OVERVIEW =====
      tabItem(tabName = "overview",
              h3(class = "section-title", "Nebraska Childhood Lead Surveillance Dashboard"),
              p(class = "section-desc",
                "This dashboard provides the Nebraska Lead Team with interactive tools to analyze childhood blood lead
           level (BLL) testing data from 2010-2023. Use it to identify geographic areas and specific addresses
           where children may be at elevated risk for lead exposure."),
              p(class = "section-desc",
                strong("Surveillance Map:"), " View testing metrics by LHD, county, census tract,
          or block group. Filter by year and toggle between child-level and address-level counts.", br(),
                strong("Priority Addresses:"), " Identify and rank addresses with multiple elevated children for
          targeted outreach and remediation referrals.", br(),
                strong("Census Data:"), " Explore neighborhood characteristics (housing age, poverty, etc.) associated
          with lead risk at the tract or block group level.", br(),
                strong("Data Tables:"), " Download complete datasets for LHDs, counties, geographies, and addresses.", br(),
                strong("Risk Model:"), " Apply predictive modeling to identify addresses most likely to have future
          elevated children, enabling proactive outreach before testing occurs."),
              p(em("Elevated BLL defined as \u22653.5 \u00b5g/dL per 2021 CDC reference value. Data includes first test per child, ages 0-6.")),
              hr(),
              fluidRow(
                valueBox(format(nrow(first_test), big.mark = ","), "Children Tested",
                         icon = icon("child"), color = "light-blue", width = 6),
                valueBox(paste0(format(sum(first_test$elevated, na.rm = TRUE), big.mark = ","), " (",
                                round(100 * mean(first_test$elevated, na.rm = TRUE), 2), "%)"),
                         "Children Elevated", icon = icon("exclamation-triangle"), color = "blue", width = 6)),
              fluidRow(
                valueBox(format(nrow(address_stats_all), big.mark = ","), "Addresses Tested",
                         icon = icon("home"), color = "green", width = 6),
                valueBox({
                  paste0(format(sum(addr_elevated_summary_all$e > 0, na.rm = TRUE), big.mark = ","),
                         " (", round(100 * mean(addr_elevated_summary_all$e > 0, na.rm = TRUE), 2), "%)")
                }, "Addresses with Elevated Child", icon = icon("home"), color = "olive", width = 6)),
              fluidRow(
                valueBox(format(nrow(address_stats), big.mark = ","), "Addresses with 2+ Children Tested",
                         icon = icon("users"), color = "yellow", width = 6),
                valueBox({
                  paste0(format(sum(address_stats$n_elevated > 0, na.rm = TRUE), big.mark = ","), " (",
                         round(100 * mean(address_stats$n_elevated > 0, na.rm = TRUE), 2), "%)")
                }, "2+ Child Addresses with Elevated", icon = icon("users"), color = "orange", width = 6))
      ),

      # ===== SURVEILLANCE MAP =====
      tabItem(tabName = "map",
              h3(class = "section-title", "Geographic Surveillance"),
              p(class = "section-desc",
                "This map displays lead testing metrics aggregated to different geographic levels. Smaller geographies
           (tracts, block groups) provide more local detail but may have more suppressed values. Click any region
           to see statistics. The trend chart shows state-level changes over time."),
              fluidRow(
                box(width = 9,
                    leafletOutput("surv_map", height = 450),
                    hr(),
                    fluidRow(
                      column(10, plotlyOutput("state_avg_plot", height = 220)),
                      column(2, checkboxInput("show_trend_ci", "Show 95% CI", FALSE))),
                    hr(),
                    h5("Annual Summary Statistics"), DTOutput("trend_summary_table")),
                box(width = 3, title = "Map Options",
                    sliderInput("year_filter", "Years:",
                                min = min(years), max = max(years),
                                value = c(min(years), max(years)), sep = ""),
                    selectInput("geo_level", "Geographic Level:",
                                c("LHD" = "LHD", "County" = "County",
                                  "Census Tract" = "Tract", "Block Group" = "Block Group")),
                    selectInput("fill_var", "Measure:",
                                c("Average BLL (\u00b5g/dL)" = "mean_bll",
                                  "Elevated Count" = "n_elevated",
                                  "% Elevated (elevated/tested)" = "pct_elevated",
                                  "Tested Count" = "n_children")),
                    radioButtons("unit_type", "Count By:", c("Children", "Addresses"), inline = TRUE),
                    hr(),
                    checkboxInput("use_quintile", "Display by Quintile", FALSE),
                    checkboxInput("show_map_ci", "Show 95% CIs in Tooltips", FALSE),
                    hr(), h5("Selected Area:"), verbatimTextOutput("click_info"),
                    hr(),
                    div(class = "state-avg-box",
                        h5("Nebraska", style = "margin: 0 0 5px 0;"), uiOutput("state_avg_text"))))
      ),

      # ===== PRIORITY ADDRESSES =====
      tabItem(tabName = "ranked",
              h3(class = "section-title", "Priority Addresses for Outreach"),
              p(class = "section-desc",
                "Addresses ranked by elevated count, average BLL, or predicted risk. Use this to prioritize outreach,
           inspection referrals, or remediation assistance. Colors indicate BLL severity.
           Numbers show rank within the selected LHD."),
              fluidRow(
                box(width = 9,
                    leafletOutput("ranked_map", height = 500), hr(),
                    downloadButton("dl_ranked", "Download Displayed Addresses"),
                    br(), br(),
                    div(style = "height: 300px; overflow-y: auto;", DTOutput("ranked_table"))),
                box(width = 3, title = "Options",
                    sliderInput("ranked_years", "Year Range:",
                                min = min(years), max = max(years),
                                value = c(min(years), max(years)), sep = ""),
                    selectInput("ranked_lhd", "LHD:", c("Statewide" = "All", lhd_names)),
                    sliderInput("top_n", "Show Top N per LHD:", 5, 50, 10, 5),
                    radioButtons("rank_by", "Rank By:",
                                 c("Elevated Count" = "n_elevated",
                                   "Average BLL" = "mean_bll",
                                   "Predicted Risk" = "pred_prob")),
                    hr(), h5("Selected Address:"),
                    verbatimTextOutput("ranked_click_info"),
                    hr(),
                    radioButtons("map_tiles", "Base Map:", c("Street" = "street", "Satellite" = "sat"), inline = TRUE)))
      ),

      # ===== CENSUS DATA =====
      tabItem(tabName = "census",
              h3(class = "section-title", "Census Neighborhood Characteristics"),
              p(class = "section-desc",
                "American Community Survey (ACS) data linked to lead testing outcomes. Areas with older housing,
           higher poverty, or more rental units may have elevated lead risk. Click any area for details."),
              fluidRow(
                box(width = 9,
                    leafletOutput("census_map", height = 380),
                    checkboxInput("census_quintile", "Display by Quintile", FALSE),
                    hr(),
                    h5("Lead Outcomes by Quintile"),
                    fluidRow(
                      column(6, plotlyOutput("quintile_pct_plot", height = 320)),
                      column(6, plotlyOutput("quintile_bll_plot", height = 320))),
                    hr(),
                    h5("Quintile Summary"),
                    downloadButton("dl_quintile", "Download Quintile Summary"),
                    div(style = "max-height: 350px; overflow-y: auto;", DTOutput("census_quintile_table")),
                    hr(),
                    h5("Area Details"),
                    downloadButton("dl_area_detail", "Download Area Details")),
                box(width = 3, title = "Options",
                    selectInput("census_var", "Map Variable:",
                                c("Pre-1980 Housing (%)" = "pct_pre1980",
                                  "Pre-1950 Housing (%)" = "pct_pre1950",
                                  "Poverty Rate (%)" = "pct_poverty",
                                  "Rent Burden \u226530% (%)" = "pct_rent_burden",
                                  "Renter-Occupied (%)" = "pct_renter",
                                  "Median Income ($)" = "median_income",
                                  "Non-White (%)" = "pct_nonwhite",
                                  "Vacant Housing (%)" = "pct_vacant",
                                  "Housing Density (log)" = "log_housing_density",
                                  "Children Under 6 (%)" = "pct_children_u6")),
                    radioButtons("census_geo", "Level:", c("Tract" = "tract", "Block Group" = "bg"), inline = TRUE),
                    hr(),
                    h5("Selected Area:"), uiOutput("census_info_ui"),
                    hr(),
                    h5("Nebraska Overall:"),
                    uiOutput("census_state_overall"),
                    hr(),
                    p(em("Data sourced from Nebraska CLPPP blood lead testing records joined with ACS 5-year estimates at the selected geography level."),
                      style = "font-size: 10px; color: #888;")))
      ),

      # ===== DATA TABLES =====
      tabItem(tabName = "tables",
              h3(class = "section-title", "Data Tables"),
              p(class = "section-desc", "Browse and download datasets. Choose to download displayed rows only or all rows."),
              fluidRow(
                column(3, sliderInput("table_years", "Years:", min(years), max(years), c(min(years), max(years)), sep = "")),
                column(3, radioButtons("table_unit", "Unit:", c("Children", "Addresses"), inline = TRUE)),
                column(3, radioButtons("table_geo", "Small Area Geography:", c("Tract" = "tract", "Block Group" = "bg"), inline = TRUE)),
                column(3, selectInput("table_rows", "Display Rows:", c(10, 25, 50, 100), selected = 25))),
              tabsetPanel(
                tabPanel("LHDs", br(),
                         fluidRow(column(6, downloadButton("dl_lhd_disp", "Download Displayed")),
                                  column(6, downloadButton("dl_lhd_all", "Download All"))),
                         br(), DTOutput("lhd_table")),
                tabPanel("Counties", br(),
                         fluidRow(column(6, downloadButton("dl_county_disp", "Download Displayed")),
                                  column(6, downloadButton("dl_county_all", "Download All"))),
                         br(), DTOutput("county_table")),
                tabPanel("Small Area", br(),
                         fluidRow(column(6, downloadButton("dl_geo_disp", "Download Displayed")),
                                  column(6, downloadButton("dl_geo_all", "Download All"))),
                         br(), DTOutput("geo_table")),
                tabPanel("Addresses", br(),
                         fluidRow(column(6, downloadButton("dl_addr_disp", "Download Displayed")),
                                  column(6, downloadButton("dl_addr_all", "Download All"))),
                         br(), DTOutput("addr_table")))
      ),

      # ===== RISK MODEL =====
      tabItem(tabName = "model",
              h3(class = "section-title", "Predictive Risk Model"),
              p(class = "section-desc",
                "This model predicts which addresses are most likely to have future children with elevated BLL.
           It uses mixed-effects logistic regression with geographic random effects to account for neighborhood clustering.
           An XGBoost comparison is always run alongside for benchmarking."),
              fluidRow(
                box(width = 4, title = "Model Settings",
                    sliderInput("train_years", "Training Years:", min(years), max(years), c(2010, 2016), sep = ""),
                    checkboxInput("train_all", "Train on all years", FALSE),
                    div(class = "warning-box",
                        em("Train on all years: Use when you want predictions for future addresses but cannot evaluate
                    model performance (no held-out test set). Useful for deployment, not validation.")),
                    hr(),
                    radioButtons("random_effect", "Random Intercept Level:",
                                 c("Census Tract" = "tract_geoid", "Block Group" = "bg_geoid",
                                   "None (Fixed Effects Only)" = "none"),
                                 selected = "tract_geoid", inline = FALSE),
                    p(em("Adds a neighborhood-specific adjustment. Tracts are larger (more stable);
                  block groups are smaller (more local). 'None' runs standard logistic regression."),
                      style = "font-size: 11px; color: #666;"),
                    hr(),
                    checkboxInput("scale_vars", "Standardize variables", TRUE),
                    div(class = "warning-box",
                        em("Standardization: Makes odds ratios comparable across variables with different scales.
                    Note: Categorical, race/ethnicity indicator, and time variables are never standardized.")),
                    hr(),
                    h5("Neighborhood Variables:"),
                    checkboxGroupInput("neighborhood_covars", NULL,
                                       choices = c("Pre-1980 Housing (%)" = "pct_pre1980",
                                                   "Pre-1950 Housing (%)" = "pct_pre1950",
                                                   "Poverty Rate (%)" = "pct_poverty",
                                                   "Rent Burden \u226530% (%)" = "pct_rent_burden",
                                                   "Renter-Occupied (%)" = "pct_renter",
                                                   "Median Income ($)" = "median_income",
                                                   "Income Category" = "income_cat",
                                                   "Non-White (%)" = "pct_nonwhite",
                                                   "Housing Density (log)" = "log_housing_density",
                                                   "Vacant Housing (%)" = "pct_vacant"),
                                       selected = "pct_pre1980"),
                    h5("Address Variables (First Child):"),
                    checkboxGroupInput("addr_covars", NULL,
                                       choices = c("First Child's BLL" = "first_bll",
                                                   "First Test Year" = "first_year",
                                                   "First Child Age (months)" = "first_age_mo",
                                                   "First Child Sex" = "first_sex",
                                                   "Sample Type" = "first_sample_type"),
                                       selected = c("first_bll", "first_year")),
                    conditionalPanel("input.addr_covars.indexOf('first_bll') > -1",
                      checkboxInput("use_log_bll", "Use log(BLL) instead of raw BLL", FALSE)),
                    h5("Race/Ethnicity (Individual):"),
                    p(em("Each group as 0/1 indicator. ~25% missing."),
                      style = "font-size: 10px; color: #888;"),
                    checkboxGroupInput("race_eth_covars", NULL,
                                       choices = c("Hispanic" = "re_hispanic",
                                                   "White (Non-Hispanic)" = "re_white_nh",
                                                   "Black (Non-Hispanic)" = "re_black_nh",
                                                   "AIAN (Non-Hispanic)" = "re_aian_nh",
                                                   "Asian (Non-Hispanic)" = "re_asian_nh",
                                                   "NHOPI (Non-Hispanic)" = "re_nhopi_nh",
                                                   "Two or More (Non-Hispanic)" = "re_multi_nh",
                                                   "Other Race (Non-Hispanic)" = "re_other_nh")),
                    p(em("XGBoost (gradient boosting) is always run alongside for comparison.
                         It is a machine learning algorithm that captures non-linear patterns."),
                      style = "font-size: 11px; color: #666;"),
                    hr(),
                    actionButton("fit_model", "Run Model", class = "btn-primary btn-block")),
                box(width = 8, title = "Results",
                    tabsetPanel(
                      tabPanel("Summary", br(), uiOutput("interpretation")),
                      tabPanel("Sample", br(),
                               checkboxInput("sample_bll_log", "Log scale for BLL plot", FALSE),
                               fluidRow(
                                 column(6, plotlyOutput("sample_plot_addr", height = 240)),
                                 column(6, plotlyOutput("sample_plot_tested", height = 240))),
                               fluidRow(
                                 column(6, plotlyOutput("sample_plot_bll",  height = 240)),
                                 column(6, plotlyOutput("sample_plot_elev", height = 240))),
                               br(), downloadButton("dl_sample_table", "Download Year-by-Year Summary")),
                      tabPanel("Odds Ratios", br(), uiOutput("or_note"),
                               plotlyOutput("coef_plot", height = 350)),
                      tabPanel("Correlations", br(),
                               p("High correlations between predictors can cause unstable estimates. Consider
                   removing one variable if correlation > 0.70.", style = "font-size: 12px;"),
                               radioButtons("cor_method", "Correlation Method:",
                                            c("Pearson" = "pearson", "Spearman" = "spearman"), inline = TRUE),
                               p(em("Pearson measures linear association; Spearman measures monotonic (rank-based) association."),
                                 style = "font-size: 10px; color: #888;"),
                               p(span("Yellow = borderline (0.50-0.70)", style = "background: #fff3cd; padding: 2px 5px;"),
                                 span("Red = high (> 0.70)", style = "background: #f8d7da; padding: 2px 5px; margin-left: 10px;")),
                               DTOutput("cor_table"), br(),
                               h5("Variance Inflation Factors (VIF)"),
                               p("VIF > 5 suggests problematic multicollinearity. VIF > 10 is severe.", style = "font-size: 12px;"),
                               DTOutput("vif_table")),
                      tabPanel("Model Comparison", br(),
                               uiOutput("model_comparison_ui"),
                               hr(),
                               h5("Cumulative Capture Curves"),
                               p("Fraction of elevated addresses identified as more addresses are screened.",
                                 style = "font-size: 12px; color: #666;"),
                               plotlyOutput("capture_overlay_plot", height = 780),
                               hr(),
                               fluidRow(
                                 column(6, h5("AUC by Test Year"),
                                        plotlyOutput("auc_year_overlay_plot", height = 320)),
                                 column(6, h5("Calibration Over Time"),
                                        plotlyOutput("cal_time_overlay_plot", height = 320)))),
                      tabPanel("Technical", br(), verbatimTextOutput("model_summary"))))
              ),
              fluidRow(
                box(width = 12, title = "Predicted Risk Map", status = "primary", solidHeader = TRUE,
                    p("Average predicted risk for addresses in each area. Darker color = higher predicted risk of elevated BLL.",
                      style = "font-size: 13px; color: #555;"),
                    fluidRow(
                      column(2, radioButtons("pred_type", "Data:", c("Test Set" = "test", "All Addresses" = "full"), inline = FALSE)),
                      column(3, selectInput("pred_geo", "Map Level:",
                                            c("County" = "county", "LHD" = "lhd",
                                              "Tract" = "tract", "Block Group" = "bg"),
                                            selected = "tract")),
                      column(3, selectInput("pred_class", "Classification:",
                                            c("Continuous" = "continuous", "Quintile" = "quintile", "Decile" = "decile"),
                                            selected = "quintile")),
                      column(2, selectInput("pred_palette", "Color Scheme:",
                                            c("Red-Yellow" = "risk", "Viridis (colorblind)" = "viridis",
                                              "Cividis (colorblind)" = "cividis"),
                                            selected = "risk"))),
                    leafletOutput("pred_map", height = 760)))
      )
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    fit = NULL, fit_full = NULL, test_preds = NULL, full_preds = NULL,
    roc_obj = NULL, comparison = NULL, cor_matrix = NULL,
    clicked_id = NULL, census_id = NULL, ranked_clicked_id = NULL,
    model_run = FALSE, is_scaled = FALSE,
    train_sample = NULL, test_sample = NULL,
    ranked_display = NULL,
    xgb_fit = NULL, xgb_roc = NULL, xgb_comparison = NULL,
    xgb_test_preds = NULL)

  set_dl <- function(id, data_func, prefix, n_head = NULL) {
    output[[id]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".csv"),
      content  = function(file) {
        d <- data_func()
        if (!is.null(n_head)) d <- head(d, n_head())
        write.csv(d, file, row.names = FALSE)
      })
  }

  use_geo <- reactive(input$use_geo_mean == "geo")

  active_data <- reactive({
    if (!is.null(input$test_selection) && input$test_selection == "max" && !is.null(max_test))
      max_test
    else
      first_test
  })

  get_address_stats <- reactive({
    if (input$addr_pop == "multi") address_stats else address_stats_all
  })
  get_addr_elev_summary <- reactive({
    if (input$addr_pop == "multi") addr_elevated_summary else addr_elevated_summary_all
  })

  observeEvent(input$addr_pop, {
    lbl <- if (input$addr_pop == "multi") "2+ Children" else "All Tested"
    showNotification(paste0("Address population: ", lbl,
                            ". Affects maps (Addresses mode), Priority Addresses, and Data Tables."),
                     type = "message", duration = 4)
  }, ignoreInit = TRUE)

  get_census_covars <- reactive({ input$neighborhood_covars %||% character(0) })

  get_all_covars <- reactive({
    c(get_census_covars(), input$addr_covars, input$race_eth_covars %||% character(0))
  })

  get_state_trend <- reactive({
    geo <- use_geo(); unit <- input$unit_type; dat <- active_data()
    if (unit == "Addresses" && input$addr_pop == "multi") {
      multi_addrs <- dat %>% group_by(address_id) %>%
        summarize(n = n(), .groups = "drop") %>% filter(n >= 2) %>% pull(address_id)
      dat <- dat %>% filter(address_id %in% multi_addrs)
    }
    state_summary(dat, by_year = TRUE, by_address = (unit == "Addresses"), use_geo = geo)
  })
  get_state_all <- reactive({
    geo <- use_geo(); unit <- input$unit_type; dat <- active_data()
    if (unit == "Addresses" && input$addr_pop == "multi") {
      multi_addrs <- dat %>% group_by(address_id) %>%
        summarize(n = n(), .groups = "drop") %>% filter(n >= 2) %>% pull(address_id)
      dat <- dat %>% filter(address_id %in% multi_addrs)
    }
    state_summary(dat, by_address = (unit == "Addresses"), use_geo = geo)
  })

  get_stats <- reactive({
    yr  <- input$year_filter; sup <- input$suppress_n; geo <- use_geo()
    dat <- active_data()
    ft  <- dat %>% filter(sample_year >= yr[1], sample_year <= yr[2])
    if (input$addr_pop == "multi") {
      multi_addrs <- ft %>% group_by(address_id) %>%
        summarize(n = n(), .groups = "drop") %>% filter(n >= 2) %>% pull(address_id)
      ft_addr <- ft %>% filter(address_id %in% multi_addrs)
    } else {
      ft_addr <- ft
    }
    make_stats <- function(g_col) {
      child <- summarize_lead(ft, g_col, sup, geo)
      addr_data <- ft_addr %>% filter(!is.na(.data[[g_col]])) %>%
        group_by(geo_id = .data[[g_col]], address_id) %>%
        summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
        group_by(geo_id) %>%
        summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")
      child %>% select(-any_of(c("n_addresses", "addr_elevated"))) %>%
        left_join(addr_data, by = "geo_id")
    }
    list(lhd = make_stats("lhd"), county = make_stats("county"),
         tract = make_stats("tract_geoid"), bg = make_stats("bg_geoid"))
  })

  # ===== SURVEILLANCE MAP =====

  output$surv_map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>% setView(-99.5, 41.5, zoom = 7)
  })

  output$state_avg_text <- renderUI({
    unit <- input$unit_type; var <- input$fill_var; yr <- input$year_filter
    dat <- active_data()
    ft_sel <- dat %>% filter(sample_year >= yr[1], sample_year <= yr[2])
    avg_sel <- state_summary(ft_sel, by_address = (unit == "Addresses"), use_geo = use_geo())
    avg_all <- get_state_all()

    fmt <- function(a, v) {
      switch(v,
        mean_bll = sprintf("%.2f \u00b5g/dL", a$mean_bll),
        n_elevated = format(a$n_elevated, big.mark = ","),
        n_children = format(a$n_children, big.mark = ","),
        pct_elevated = sprintf("%.2f%%", a$pct_elevated),
        "N/A")
    }

    is_subset <- yr[1] != min(years) || yr[2] != max(years)
    if (is_subset) {
      ft_non_sel <- dat %>% filter(sample_year < yr[1] | sample_year > yr[2])
      avg_non_sel <- state_summary(ft_non_sel, by_address = (unit == "Addresses"), use_geo = use_geo())
      non_sel_yrs <- sort(unique(ft_non_sel$sample_year))
      non_sel_label <- if (length(non_sel_yrs) > 0) {
        paste0(min(non_sel_yrs), "\u2013", max(non_sel_yrs), " (excl. selected)")
      } else "Other Years"
      txt <- paste0("<b>", yr[1], "\u2013", yr[2], ":</b> ", fmt(avg_sel, var),
                    "<br><em style='font-size:12px;'>", non_sel_label, ": ", fmt(avg_non_sel, var), "</em>",
                    "<br><em style='font-size:12px;'>All Years: ", fmt(avg_all, var), "</em>")
    } else {
      txt <- paste0("<b>All Years:</b> ", fmt(avg_all, var))
    }
    mean_label <- if (use_geo()) "Geometric Mean" else "Arithmetic Mean"
    HTML(paste0(txt, "<br><em style='font-size:12px;'>BLL: ", mean_label, "</em>"))
  })

  # State trend plot
  output$state_avg_plot <- renderPlotly({
    var <- input$fill_var; yr_sel <- input$year_filter; unit <- input$unit_type
    show_ci_val <- input$show_trend_ci; test_sel <- input$test_selection
    df <- get_state_trend(); dat <- active_data()
    if (!var %in% names(df)) return(NULL)
    df$y_val <- df[[var]]
    if (!"sample_year" %in% names(df)) return(NULL)

    y_lab <- get_label(var)
    df <- df %>% filter(sample_year >= yr_sel[1], sample_year <= yr_sel[2])
    if (nrow(df) == 0) return(NULL)

    # Compute CI
    show_ci <- isTRUE(input$show_trend_ci)
    if (show_ci) {
      df$ci_lo <- NA_real_; df$ci_hi <- NA_real_
      for (i in seq_len(nrow(df))) {
        ft_yr <- dat %>% filter(sample_year == df$sample_year[i])
        n_t <- nrow(ft_yr); sv <- df$y_val[i]
        if (n_t > 0 && !is.na(sv)) {
          if (var == "pct_elevated") {
            p_hat <- sv / 100; se <- sqrt(p_hat * (1 - p_hat) / n_t) * 100
            df$ci_lo[i] <- max(0, sv - 1.96 * se); df$ci_hi[i] <- sv + 1.96 * se
          } else if (var == "mean_bll") {
            se <- sd(ft_yr$result, na.rm = TRUE) / sqrt(n_t)
            df$ci_lo[i] <- sv - 1.96 * se; df$ci_hi[i] <- sv + 1.96 * se
          } else if (var %in% c("n_elevated", "n_children")) {
            p_hat <- mean(ft_yr$elevated, na.rm = TRUE)
            se <- sqrt(n_t * p_hat * (1 - p_hat))
            df$ci_lo[i] <- max(0, sv - 1.96 * se); df$ci_hi[i] <- sv + 1.96 * se
          }
        }
      }
    }

    fmt_val <- function(v) {
      if (var %in% c("n_elevated", "n_children")) format(round(v), big.mark = ",")
      else round(v, 2)
    }

    df$series <- "Nebraska"
    df$hover_text <- paste0("Nebraska\nYear: ", df$sample_year, "\n", y_lab, ": ", sapply(df$y_val, fmt_val))
    if (show_ci && all(c("ci_lo", "ci_hi") %in% names(df)) && any(!is.na(df$ci_lo)))
      df$hover_text <- paste0(df$hover_text, "\n95% CI: [", sapply(df$ci_lo, fmt_val),
                              ", ", sapply(df$ci_hi, fmt_val), "]")

    # Area overlay
    area_df <- NULL
    if (!is.null(rv$clicked_id)) {
      g_col <- switch(input$geo_level, LHD = "lhd", County = "county",
                      Tract = "tract_geoid", `Block Group` = "bg_geoid")
      ft_area <- dat %>% filter(.data[[g_col]] == rv$clicked_id,
                                sample_year >= yr_sel[1], sample_year <= yr_sel[2])
      if (nrow(ft_area) > 0) {
        area_trend <- state_summary(ft_area, by_year = TRUE,
                                    by_address = (unit == "Addresses"), use_geo = use_geo())
        if (var %in% names(area_trend)) area_trend$y_val <- area_trend[[var]]
        else area_trend$y_val <- NA_real_
        area_df <- area_trend %>% filter(!is.na(y_val))
        if (nrow(area_df) > 0) {
          area_df$series <- rv$clicked_id
          if (show_ci) {
            area_df$ci_lo <- NA_real_; area_df$ci_hi <- NA_real_
            for (i in seq_len(nrow(area_df))) {
              ft_yr <- ft_area %>% filter(sample_year == area_df$sample_year[i])
              n_t <- nrow(ft_yr); sv <- area_df$y_val[i]
              if (n_t > 0 && !is.na(sv)) {
                if (var == "pct_elevated") {
                  p_hat <- sv / 100; se <- sqrt(p_hat * (1 - p_hat) / n_t) * 100
                  area_df$ci_lo[i] <- max(0, sv - 1.96 * se); area_df$ci_hi[i] <- sv + 1.96 * se
                } else if (var == "mean_bll") {
                  se <- sd(ft_yr$result, na.rm = TRUE) / sqrt(n_t)
                  area_df$ci_lo[i] <- sv - 1.96 * se; area_df$ci_hi[i] <- sv + 1.96 * se
                } else if (var %in% c("n_elevated", "n_children")) {
                  p_hat <- mean(ft_yr$elevated, na.rm = TRUE)
                  se <- sqrt(n_t * p_hat * (1 - p_hat))
                  area_df$ci_lo[i] <- max(0, sv - 1.96 * se); area_df$ci_hi[i] <- sv + 1.96 * se
                }
              }
            }
          }
          area_df$hover_text <- paste0(rv$clicked_id, "\nYear: ", area_df$sample_year,
                                       "\n", y_lab, ": ", sapply(area_df$y_val, fmt_val))
          if (show_ci && all(c("ci_lo", "ci_hi") %in% names(area_df)) && any(!is.na(area_df$ci_lo)))
            area_df$hover_text <- paste0(area_df$hover_text, "\n95% CI: [",
                                         sapply(area_df$ci_lo, fmt_val), ", ", sapply(area_df$ci_hi, fmt_val), "]")
        }
      }
    }

    plot_df <- df %>% select(sample_year, y_val, series, hover_text, any_of(c("ci_lo", "ci_hi")))
    if (!is.null(area_df) && nrow(area_df) > 0)
      plot_df <- bind_rows(plot_df, area_df %>% select(sample_year, y_val, series, hover_text))
    plot_df <- plot_df %>% arrange(series, sample_year) %>% filter(!is.na(y_val))

    p <- plot_ly()

    # CI ribbon
    if (show_ci && all(c("ci_lo", "ci_hi") %in% names(df)) && any(!is.na(df$ci_lo))) {
      ci_df <- df %>% filter(!is.na(ci_lo)) %>% arrange(sample_year)
      p <- p %>% add_trace(x = ci_df$sample_year, y = ci_df$ci_hi, type = "scatter", mode = "lines",
                            line = list(width = 0), showlegend = FALSE, hoverinfo = "skip") %>%
        add_trace(x = ci_df$sample_year, y = ci_df$ci_lo, type = "scatter", mode = "lines",
                  fill = "tonexty", fillcolor = "rgba(44,127,184,0.15)",
                  line = list(width = 0), showlegend = FALSE, hoverinfo = "skip")
    }

    ne_df <- plot_df %>% filter(series == "Nebraska")
    if (nrow(ne_df) > 0) {
      p <- p %>% add_trace(x = ne_df$sample_year, y = ne_df$y_val, type = "scatter",
                            mode = "lines+markers", name = "Nebraska",
                            line = list(color = "#2c7fb8", width = 2.5),
                            marker = list(color = "#2c7fb8", size = 6),
                            text = ne_df$hover_text, hoverinfo = "text")
    }

    area_plot <- plot_df %>% filter(series != "Nebraska")
    if (nrow(area_plot) > 0) {
      area_name <- area_plot$series[1]
      if (show_ci && !is.null(area_df) && all(c("ci_lo", "ci_hi") %in% names(area_df)) && any(!is.na(area_df$ci_lo))) {
        aci_df <- area_df %>% filter(!is.na(ci_lo)) %>% arrange(sample_year)
        p <- p %>% add_trace(x = aci_df$sample_year, y = aci_df$ci_hi, type = "scatter", mode = "lines",
                              line = list(width = 0), showlegend = FALSE, hoverinfo = "skip") %>%
          add_trace(x = aci_df$sample_year, y = aci_df$ci_lo, type = "scatter", mode = "lines",
                    fill = "tonexty", fillcolor = "rgba(215,48,39,0.12)",
                    line = list(width = 0), showlegend = FALSE, hoverinfo = "skip")
      }
      p <- p %>% add_trace(x = area_plot$sample_year, y = area_plot$y_val, type = "scatter",
                            mode = "lines+markers", name = area_name,
                            line = list(color = "#d73027", width = 2.5),
                            marker = list(color = "#d73027", size = 6),
                            text = area_plot$hover_text, hoverinfo = "text")
    }

    p %>% layout(title = list(text = paste("Trend:", y_lab), font = list(size = 14, weight = "bold")),
                 xaxis = list(title = "Year"), yaxis = list(title = y_lab),
                 legend = list(orientation = "h", y = -0.15), margin = list(t = 45, b = 50))
  })

  # Trend summary table
  output$trend_summary_table <- renderDT({
    var <- input$fill_var; unit <- input$unit_type; sup <- input$suppress_n
    dat <- active_data()
    g_col <- switch(input$geo_level, LHD = "lhd", County = "county",
                    Tract = "tract_geoid", `Block Group` = "bg_geoid")
    trend_df <- get_state_trend()
    display_var <- if (unit == "Addresses") switch(var, n_children = "n_addresses",
                                                    n_elevated = "addr_elevated", var) else var
    yr_rows <- lapply(years, function(yr) {
      state_row <- trend_df %>% filter(sample_year == yr)
      if (nrow(state_row) == 0) return(NULL)
      geo_stats <- summarize_lead(dat %>% filter(sample_year == yr), g_col, sup, use_geo())
      vals <- geo_stats[[display_var]][!is.na(geo_stats[[display_var]]) & !geo_stats$suppressed]
      ft_yr <- dat %>% filter(sample_year == yr); n_t <- nrow(ft_yr)
      state_val <- state_row[[var]]
      ci_lo <- NA_real_; ci_hi <- NA_real_
      if (n_t > 0 && !is.na(state_val)) {
        sv <- state_val
        if (var == "pct_elevated") {
          p_hat <- sv / 100; se <- sqrt(p_hat * (1 - p_hat) / n_t) * 100
          ci_lo <- round(max(0, sv - 1.96 * se), 2); ci_hi <- round(sv + 1.96 * se, 2)
        } else if (var == "mean_bll") {
          se <- sd(ft_yr$result, na.rm = TRUE) / sqrt(n_t)
          ci_lo <- round(sv - 1.96 * se, 2); ci_hi <- round(sv + 1.96 * se, 2)
        } else if (var %in% c("n_elevated", "n_children")) {
          p_hat <- mean(ft_yr$elevated, na.rm = TRUE); se <- sqrt(n_t * p_hat * (1 - p_hat))
          ci_lo <- round(max(0, sv - 1.96 * se), 0); ci_hi <- round(sv + 1.96 * se, 0)
        }
      }
      n_suppressed <- sum(!is.na(geo_stats[[display_var]]) & geo_stats$suppressed, na.rm = TRUE)
      data.frame(Year = yr, `State Value` = round(state_val, 2),
                 `95% CI Low` = ci_lo, `95% CI High` = ci_hi,
                 `Std Dev` = if (length(vals) > 1) round(sd(vals, na.rm = TRUE), 2) else NA_real_,
                 `Suppressed/NA` = n_suppressed, check.names = FALSE)
    })
    datatable(do.call(rbind, yr_rows),
              options = list(dom = "t", pageLength = 20, scrollY = "150px"), rownames = FALSE)
  })

  # Surveillance map observer
  observeEvent(c(input$geo_level, input$fill_var, input$unit_type, input$year_filter,
                 input$suppress_n, input$use_quintile, input$use_geo_mean, input$addr_pop,
                 input$test_selection, input$show_map_ci), {
    stats <- get_stats(); var <- input$fill_var; unit <- input$unit_type
    display_var <- if (unit == "Addresses") {
      switch(var, n_children = "n_addresses", n_elevated = "addr_elevated", var)
    } else var
    stat_data <- switch(input$geo_level,
                        LHD = stats$lhd, County = stats$county, Tract = stats$tract, `Block Group` = stats$bg)
    if (is.null(stat_data) || nrow(stat_data) == 0) return()
    geo_bounds_list <- list(LHD = lhd_bounds, County = ne_counties, Tract = ne_tracts, `Block Group` = ne_bgs)
    geo_key_list    <- list(LHD = "lhd", County = "NAME", Tract = "GEOID", `Block Group` = "GEOID")
    geo_bounds <- geo_bounds_list[[input$geo_level]]
    if (is.null(geo_bounds)) return()
    id_col <- geo_key_list[[input$geo_level]]
    geo <- geo_bounds %>% left_join(stat_data, by = setNames("geo_id", id_col))
    geo_label <- switch(input$geo_level, LHD = "LHD", County = "County",
                        Tract = "Census Tract", `Block Group` = "Block Group")
    if (!display_var %in% names(geo)) return()
    vals <- geo[[display_var]][!is.na(geo[[display_var]]) & is.finite(geo[[display_var]]) &
                                 (is.null(geo$suppressed) | !geo$suppressed)]
    if (length(vals) == 0) vals <- c(0, 1)
    var_label <- get_label(var)
    if (input$use_quintile) {
      qp <- build_quintile_palette(vals)
      geo$fill_color <- assign_fill_colors(geo, display_var, function(v) {
        qp$colors[max(1, min(findInterval(v, qp$breaks, rightmost.closed = TRUE), 5))]
      })
      legend_labels <- c(qp$labels, "Suppressed", "No Data")
      legend_colors <- c(qp$colors, "#666666", "#cccccc")
    } else {
      domain <- switch(var, mean_bll = c(0, 5), pct_elevated = c(0, 30),
                       c(0, max(vals, na.rm = TRUE)))
      pal <- colorNumeric("YlOrRd", domain = domain, na.color = "#cccccc")
      geo$fill_color <- assign_fill_colors(geo, display_var, function(v) pal(min(v, domain[2])))
      breaks <- seq(domain[1], domain[2], length.out = 6)
      fmt_brk <- function(x) {
        if (var %in% c("n_elevated", "n_children", "n_addresses", "addr_elevated")) format(round(x), big.mark = ",")
        else if (var == "pct_elevated") sprintf("%.1f%%", x)
        else format(round(x, 1), big.mark = ",")
      }
      legend_labels <- c(sapply(1:5, function(i) paste0(fmt_brk(breaks[i]), "-", fmt_brk(breaks[i+1]))),
                         "Suppressed", "No Data")
      legend_colors <- c(pal(breaks[1:5] + diff(breaks[1:2])/2), "#666666", "#cccccc")
    }
    sup_flag <- ifelse(!is.null(geo$suppressed) & !is.na(geo$suppressed) & geo$suppressed, " [Suppressed]", "")
    ml <- bll_mean_label(use_geo())
    na_bll <- ifelse(!is.na(geo$mean_bll), paste0(sprintf("%.2f", geo$mean_bll), " \u00b5g/dL ", ml), "N/A")

    # Compute CIs for map tooltips if enabled
    ci_html_pct <- ""; ci_html_bll <- ""
    if (isTRUE(input$show_map_ci)) {
      ci_pct_lo <- ci_pct_hi <- ci_bll_lo <- ci_bll_hi <- rep(NA_real_, nrow(geo))
      for (i in seq_len(nrow(geo))) {
        n_t <- geo$n_children[i]; n_e <- geo$n_elevated[i]; bll <- geo$mean_bll[i]
        if (!is.na(n_t) && n_t > 0) {
          p_hat <- n_e / n_t; z <- 1.96
          denom <- 1 + z^2 / n_t
          center <- (p_hat + z^2 / (2 * n_t)) / denom
          margin <- z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * n_t)) / n_t) / denom
          ci_pct_lo[i] <- max(0, (center - margin)) * 100
          ci_pct_hi[i] <- min(100, (center + margin)) * 100
          if (!is.na(bll) && n_t >= 2) {
            se_bll <- bll * 0.5 / sqrt(n_t)
            ci_bll_lo[i] <- max(0, bll - 1.96 * se_bll)
            ci_bll_hi[i] <- bll + 1.96 * se_bll
          }
        }
      }
      geo$ci_pct <- ifelse(!is.na(ci_pct_lo), paste0("<em style='font-size:11px; color:#555;'>95% CI: [",
                           sprintf("%.1f", ci_pct_lo), "%, ", sprintf("%.1f", ci_pct_hi), "%]</em>"), "")
      geo$ci_bll <- ifelse(!is.na(ci_bll_lo), paste0("<em style='font-size:11px; color:#555;'>95% CI: [",
                           sprintf("%.2f", ci_bll_lo), ", ", sprintf("%.2f", ci_bll_hi), "]</em>"), "")
    } else {
      geo$ci_pct <- ""; geo$ci_bll <- ""
    }
    if (unit == "Children") {
      geo$tip <- paste0("<b>", geo_label, ": ", geo[[id_col]], "</b>", sup_flag, "<br>",
                        "Elevated Count: ", format(geo$n_elevated, big.mark = ","), "<br>",
                        "Tested Count: ", format(geo$n_children, big.mark = ","), "<br>",
                        "% Elevated: ", ifelse(!is.na(geo$pct_elevated), sprintf("%.1f%%", geo$pct_elevated), "N/A"),
                        ifelse(geo$ci_pct != "", paste0("<br>", geo$ci_pct), ""), "<br>",
                        "Average BLL: ", na_bll,
                        ifelse(geo$ci_bll != "", paste0("<br>", geo$ci_bll), ""))
    } else {
      pct <- ifelse(!is.na(geo$n_addresses) & geo$n_addresses > 0,
                    sprintf("%.1f%%", 100 * geo$addr_elevated / geo$n_addresses), "N/A")
      geo$tip <- paste0("<b>", geo_label, ": ", geo[[id_col]], "</b>", sup_flag, "<br>",
                        "Addresses Elevated: ", format(geo$addr_elevated, big.mark = ","), "<br>",
                        "Addresses Tested: ", format(geo$n_addresses, big.mark = ","), "<br>",
                        "% Elevated: ", pct,
                        ifelse(geo$ci_pct != "", paste0("<br>", geo$ci_pct), ""), "<br>",
                        "Average BLL: ", na_bll,
                        ifelse(geo$ci_bll != "", paste0("<br>", geo$ci_bll), ""))
    }
    leafletProxy("surv_map") %>% clearShapes() %>% clearControls() %>%
      addPolygons(data = geo, fillColor = ~fill_color, fillOpacity = 0.7,
                  color = "white", weight = 1, layerId = as.character(geo[[id_col]]),
                  label = lapply(geo$tip, HTML)) %>%
      addLegend("bottomright", colors = legend_colors, labels = legend_labels,
                title = var_label, opacity = 0.8)
  }, ignoreNULL = FALSE)

  observeEvent(input$surv_map_shape_click, { rv$clicked_id <- input$surv_map_shape_click$id })

  output$click_info <- renderPrint({
    if (is.null(rv$clicked_id)) return(cat("Click a region for details"))
    geo_label <- switch(input$geo_level, LHD = "LHD", County = "County",
                        Tract = "Tract", `Block Group` = "Block Group")
    yr <- input$year_filter
    d <- switch(input$geo_level, LHD = get_stats()$lhd, County = get_stats()$county,
                Tract = get_stats()$tract, `Block Group` = get_stats()$bg) %>%
      filter(geo_id == rv$clicked_id)
    if (nrow(d) == 0) return(cat("No data"))
    ft_sel <- active_data() %>% filter(sample_year >= yr[1], sample_year <= yr[2])
    state <- state_summary(ft_sel, by_address = (input$unit_type == "Addresses"), use_geo = use_geo())
    ml <- bll_mean_label(use_geo())
    cat(geo_label, ":", rv$clicked_id, "\n")
    if (d$suppressed) cat("[Suppressed]\n")
    cat("\n")
    fmt_n <- function(x) format(x, big.mark = ",")
    fmt_p <- function(x) if (!is.na(x)) sprintf("%.2f%%", x) else "N/A"
    fmt_b <- function(x) if (!is.na(x)) sprintf("%.2f", x) else "N/A"
    cat(sprintf("%-22s %-18s\n", "--- Selected Area ---", "--- Nebraska ---"))
    if (input$unit_type == "Children") {
      cat(sprintf("Elevated:   %-12s Elevated:   %s\n", fmt_n(d$n_elevated), fmt_n(state$n_elevated)))
      cat(sprintf("Tested:     %-12s Tested:     %s\n", fmt_n(d$n_children), fmt_n(state$n_children)))
      cat(sprintf("%% Elevated: %-12s %% Elevated: %s\n", fmt_p(d$pct_elevated), fmt_p(state$pct_elevated)))
      cat(sprintf("Avg BLL:    %-12s Avg BLL:    %s\n", fmt_b(d$mean_bll), fmt_b(state$mean_bll)))
    } else {
      pct_d <- if (d$n_addresses > 0) 100 * d$addr_elevated / d$n_addresses else NA_real_
      cat(sprintf("Addr Elev:  %-12s Elevated:   %s\n", fmt_n(d$addr_elevated), fmt_n(state$n_elevated)))
      cat(sprintf("Addr Tested:%-12s Tested:     %s\n", fmt_n(d$n_addresses), fmt_n(state$n_children)))
      cat(sprintf("%% Elevated: %-12s %% Elevated: %s\n", fmt_p(pct_d), fmt_p(state$pct_elevated)))
      cat(sprintf("Avg BLL:    %-12s Avg BLL:    %s\n", fmt_b(d$mean_bll), fmt_b(state$mean_bll)))
    }
    cat(sprintf("BLL Average: %s\n", ml))
    cat(sprintf("Years: %d-%d\n", yr[1], yr[2]))
  })

  # ===== PRIORITY ADDRESSES =====

  output$ranked_map <- renderLeaflet({
    tile <- if (input$map_tiles == "sat") providers$Esri.WorldImagery else providers$CartoDB.Positron
    leaflet() %>% addProviderTiles(tile) %>% setView(-99.5, 41.5, zoom = 7) %>%
      addPolygons(data = lhd_bounds, fillColor = "transparent", color = "#2c7fb8", weight = 2, label = ~lhd) %>%
      addLegend("bottomright",
                colors = c("#67000d", "#a50f15", "#d73027", "#fc8d59", "#c9a84c"),
                labels = c("\u226520 \u00b5g/dL", "10\u201320", "5\u201310", "3.5\u20135", "<3.5"),
                title = "BLL Color Key", opacity = 0.9)
  })
  observeEvent(input$map_tiles, {
    tile <- if (input$map_tiles == "sat") providers$Esri.WorldImagery else providers$CartoDB.Positron
    leafletProxy("ranked_map") %>% clearTiles() %>% addProviderTiles(tile)
  })

  get_ranked_addresses <- reactive({
    addr <- get_address_stats() %>% filter(!is.na(lhd))
    yr <- input$ranked_years
    addr <- addr %>% filter(last_year >= yr[1], first_year <= yr[2])
    if (use_geo()) {
      geo_means <- active_data() %>%
        filter(sample_year >= yr[1], sample_year <= yr[2]) %>%
        group_by(address_id) %>%
        summarize(mean_bll = calc_mean_bll(result, TRUE), .groups = "drop")
      addr <- addr %>% select(-mean_bll) %>% left_join(geo_means, by = "address_id")
    }
    if (!is.null(rv$full_preds)) {
      addr <- addr %>% left_join(rv$full_preds %>% select(address_id, pred_prob = pred) %>% distinct(address_id, .keep_all = TRUE), by = "address_id")
    } else {
      addr$pred_prob <- NA_real_
    }
    if (input$rank_by == "pred_prob" && all(is.na(addr$pred_prob))) return(addr[0, ])
    if (input$ranked_lhd != "All") addr <- addr %>% filter(lhd == input$ranked_lhd)
    addr <- addr %>% mutate(risk_decile = ifelse(!is.na(pred_prob), ntile(pred_prob, 10), NA_integer_))
    sort_var <- input$rank_by
    addr %>%
      group_by(lhd) %>%
      arrange(desc(.data[[sort_var]]), .by_group = TRUE) %>%
      mutate(rank = row_number()) %>%
      ungroup() %>%
      filter(rank <= input$top_n) %>%
      mutate(color = get_bll_color(mean_bll))
  })

  observeEvent(c(input$ranked_lhd, input$top_n, input$rank_by, input$ranked_years,
                 input$use_geo_mean, input$addr_pop), {
    addr <- get_ranked_addresses()
    if (is.null(addr) || nrow(addr) == 0) {
      if (input$rank_by == "pred_prob")
        showNotification("Run the Risk Model first to rank by Predicted Risk.", type = "warning")
      return()
    }
    rv$ranked_display <- addr
    ml <- bll_mean_label(use_geo())
    addr$tip <- paste0(
      "<b>#", addr$rank, " - ", addr$street_address, "</b><br>",
      "LHD: ", addr$lhd, "<br>County: ", addr$county, "<br>",
      "Children Tested: ", format(addr$n_children, big.mark = ","), "<br>",
      "Children Elevated: ", format(addr$n_elevated, big.mark = ","), "<br>",
      "Average BLL ", ml, ": ", round(addr$mean_bll, 2), " \u00b5g/dL<br>",
      "Max BLL: ", round(addr$max_bll, 2), " \u00b5g/dL<br>",
      "First Year: ", addr$first_year, " | Last Year: ", addr$last_year,
      ifelse(!is.na(addr$pred_prob),
             paste0("<br>Risk: ", round(addr$pred_prob * 100, 1), "%"), ""))
    leafletProxy("ranked_map") %>% clearMarkers() %>%
      addCircleMarkers(data = addr, lng = ~lng, lat = ~lat, radius = 8,
                       fillColor = ~color, color = "#333", weight = 1, fillOpacity = 0.85,
                       layerId = ~address_id, label = lapply(addr$tip, HTML),
                       options = pathOptions(pane = "markerPane")) %>%
      addLabelOnlyMarkers(data = addr, lng = ~lng, lat = ~lat,
                          label = as.character(addr$rank),
                          labelOptions = labelOptions(noHide = TRUE, direction = "center",
                                                      textOnly = TRUE, textsize = "10px",
                                                      style = list("font-weight" = "bold", "color" = "white",
                                                                   "pointer-events" = "none")))
  }, ignoreNULL = FALSE)

  observeEvent(input$ranked_map_marker_click, { rv$ranked_clicked_id <- input$ranked_map_marker_click$id })

  output$ranked_click_info <- renderPrint({
    if (is.null(rv$ranked_clicked_id)) return(cat("Click a marker for details"))
    d <- get_ranked_addresses() %>% filter(address_id == rv$ranked_clicked_id)
    if (nrow(d) == 0) return(cat("No data for this address"))
    ml <- bll_mean_label(use_geo())
    cat("Address:", d$street_address[1], "\n")
    cat("LHD:", d$lhd[1], "\nCounty:", d$county[1], "\n\n")
    cat("Children Tested:", d$n_children[1], "\n")
    cat("Elevated:", d$n_elevated[1], "\n")
    cat("Avg BLL", ml, ":", round(d$mean_bll[1], 2), "\n")
    cat("Max BLL:", round(d$max_bll[1], 2), "\n")
    cat("First Test Year:", d$first_year[1], "\nLast Test Year:", d$last_year[1], "\n")
    if (!is.na(d$pred_prob[1])) {
      cat("Risk Decile:", d$risk_decile[1], "\n")
      cat("Predicted Risk:", sprintf("%.1f%%", d$pred_prob[1] * 100), "\n")
    }
  })

  # View Tests modal
  show_address_tests <- function(addr_id) {
    src <- if (!is.null(all_tests)) all_tests else first_test
    tests <- src %>% filter(address_id == addr_id)
    if (nrow(tests) == 0) { showNotification("No test records found", type = "warning"); return() }
    addr_label <- tests$street_address[1]
    ensure_col <- function(df, target, patterns) {
      if (target %in% names(df)) return(df)
      found <- intersect(patterns, names(df))
      if (length(found) == 0) found <- names(df)[tolower(names(df)) %in% tolower(patterns)]
      if (length(found) > 0) { df[[target]] <- df[[found[1]]]; return(df) }
      df[[target]] <- NA_character_; df
    }
    tests <- ensure_col(tests, "SEX", c("Sex", "sex", "SEX", "GENDER", "Gender", "gender",
                                           "SEX_CD", "sex_cd", "first_sex", "PatientSex", "PATIENT_SEX"))
    tests <- ensure_col(tests, "RACE_ETH", c("Race_Eth", "race_eth", "RACE_ETH",
                                               "race_ethnicity_category", "RACE", "Race", "race",
                                               "RaceEthnicity", "RACE_ETHNICITY", "PatientRace", "PATIENT_RACE"))
    tests <- ensure_col(tests, "SAMPLE_TYPE", c("SampleType", "sample_type", "SAMPLE_TYPE",
                                                  "first_sample_type", "SpecimenType", "SPECIMEN_TYPE"))
    tests <- tests %>% mutate(
      SEX = if_else(is.na(SEX) | SEX == "", "Not recorded", as.character(SEX)),
      RACE_ETH = if_else(is.na(RACE_ETH) | RACE_ETH == "", "Not recorded", as.character(RACE_ETH)),
      SAMPLE_TYPE = if_else(is.na(SAMPLE_TYPE) | SAMPLE_TYPE == "", "Not recorded", as.character(SAMPLE_TYPE)))
    cols <- intersect(c("PATIENT_LOCAL_ID", "SAMPLE_DATE", "result", "elevated",
                         "AGE_YR", "AGE_MO", "SEX", "RACE_ETH", "SAMPLE_TYPE", "sample_year"), names(tests))
    display <- tests %>% select(all_of(cols)) %>% arrange(SAMPLE_DATE)
    sex_summary <- if (any(display$SEX != "Not recorded")) {
      tbl <- table(display$SEX[display$SEX != "Not recorded"])
      paste("Sex:", paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", "))
    } else "Sex: Not recorded in source data"
    race_summary <- if (any(display$RACE_ETH != "Not recorded")) {
      tbl <- table(display$RACE_ETH[display$RACE_ETH != "Not recorded"])
      paste("Race/Ethnicity:", paste(paste0(names(tbl), " (", tbl, ")"), collapse = ", "))
    } else "Race/Ethnicity: Not recorded in source data"
    showModal(modalDialog(
      title = paste("All Test Records:", addr_label), size = "l", easyClose = TRUE,
      p(em(paste(nrow(display), "total test records at this address"))),
      p(em(sex_summary), style = "font-size: 11px; color: #666;"),
      p(em(race_summary), style = "font-size: 11px; color: #666;"),
      DTOutput("modal_test_table"),
      footer = tagList(downloadButton("dl_modal_tests", "Download"), modalButton("Close"))))
    output$modal_test_table <- renderDT({
      datatable(display, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    })
    output$dl_modal_tests <- downloadHandler(
      filename = function() paste0("tests_", addr_id, "_", Sys.Date(), ".csv"),
      content = function(file) write.csv(display, file, row.names = FALSE))
  }

  output$ranked_table <- renderDT({
    if (is.null(rv$ranked_display)) return(NULL)
    d <- rv$ranked_display
    ml <- bll_mean_label(use_geo())
    bll_col <- paste("Avg BLL", ml)
    if (!is.null(rv$full_preds) && "single_child_flag" %in% names(rv$full_preds)) {
      sc_flags <- rv$full_preds %>% select(address_id, single_child_flag) %>% distinct()
      d <- d %>% left_join(sc_flags, by = "address_id")
    } else { d$single_child_flag <- FALSE }
    d$single_child_flag[is.na(d$single_child_flag)] <- FALSE
    tbl <- d %>%
      transmute(Rank = rank, Address = street_address, LHD = lhd, County = county,
                `Children` = n_children, `Elevated` = n_elevated,
                !!bll_col := round(mean_bll, 2), `Max BLL` = round(max_bll, 2),
                `Risk Decile` = if_else(!is.na(risk_decile), as.character(risk_decile), "\u2014"),
                `Pred Risk` = if_else(!is.na(pred_prob),
                  paste0(sprintf("%.1f%%", pred_prob * 100), if_else(single_child_flag, "*", "")), "\u2014"),
                `View Tests` = paste0('<button class="btn btn-xs btn-info" ',
                                      'onclick="Shiny.setInputValue(\'view_tests_click\', \'',
                                      address_id, '\', {priority: \'event\'})">View</button>'),
                `Risk` = paste0('<button class="btn btn-xs btn-warning" ',
                                'onclick="Shiny.setInputValue(\'risk_breakdown_click\', \'',
                                address_id, '\', {priority: \'event\'})"',
                                if_else(!is.na(pred_prob), '>', ' disabled title="Run model first">'),
                                if_else(!is.na(pred_prob), 'Explain', '\u2014'), '</button>'))
    datatable(tbl, escape = FALSE, selection = "single",
              options = list(pageLength = 10, scrollX = TRUE, dom = "tp"), rownames = FALSE,
              caption = if (any(d$single_child_flag, na.rm = TRUE))
                htmltools::tags$caption(style = "font-size:11px; color:#888; font-style:italic;",
                  "* Single-child address: predicted risk shown but not included in model validation.") else NULL)
  })
  observeEvent(input$view_tests_click, { show_address_tests(input$view_tests_click) })

  set_dl("dl_ranked", function() {
    if (is.null(rv$ranked_display)) return(data.frame())
    rv$ranked_display %>% select(-any_of(c("tip", "color")))
  }, "priority_addresses")

  # ===== CENSUS DATA =====

  get_census_geo_data <- function() {
    geo <- input$census_geo; var <- input$census_var; sup <- input$suppress_n
    census_df <- if (geo == "tract") tract_census else bg_census
    geo_bounds <- if (geo == "tract") ne_tracts else ne_bgs
    geo_col    <- if (geo == "tract") "tract_geoid" else "bg_geoid"
    id_col <- "GEOID"
    dat <- active_data()
    lead_df <- dat %>%
      filter(!is.na(.data[[geo_col]]), nchar(.data[[geo_col]]) >= 11) %>%
      group_by(geo_id = .data[[geo_col]]) %>%
      summarize(n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
                mean_bll = calc_mean_bll(result, use_geo()), .groups = "drop") %>%
      mutate(pct_elevated = 100 * n_elevated / n_children, suppressed = n_children < sup)
    census_join_col <- if (geo == "tract") "tract_geoid" else "bg_geoid"
    var_missing <- !is.null(census_df) && (!var %in% names(census_df) || all(is.na(census_df[[var]])))
    if (var_missing && !is.null(tract_census) && var %in% names(tract_census) && geo == "bg") {
      tract_vals <- tract_census %>% select(tract_geoid, all_of(var))
      census_df <- census_df %>%
        mutate(.tract_id = substr(bg_geoid, 1, 11)) %>%
        left_join(tract_vals, by = c(".tract_id" = "tract_geoid"), suffix = c("", ".tract")) %>%
        select(-.tract_id)
      if (paste0(var, ".tract") %in% names(census_df)) {
        census_df[[var]] <- coalesce(census_df[[var]], census_df[[paste0(var, ".tract")]])
        census_df[[paste0(var, ".tract")]] <- NULL
      }
    }
    combined <- lead_df %>% left_join(census_df, by = setNames(census_join_col, "geo_id"))
    geo_sf <- geo_bounds %>% left_join(combined, by = setNames("geo_id", id_col))
    list(geo_sf = geo_sf, combined = combined, id_col = id_col)
  }

  output$census_map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>% setView(-99.5, 41.5, zoom = 7)
  })

  observeEvent(c(input$census_var, input$census_geo, input$suppress_n, input$census_quintile, input$use_geo_mean), {
    cd <- get_census_geo_data()
    geo <- cd$geo_sf; var <- input$census_var
    if (is.null(geo) || !var %in% names(geo)) return()
    vals <- geo[[var]][!is.na(geo[[var]]) & is.finite(geo[[var]])]
    if (length(vals) == 0) return()
    var_label <- get_label(var)
    if (input$census_quintile) {
      qp <- build_quintile_palette(vals, c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"))
      geo$fill_color <- assign_fill_colors(geo, var, function(v) {
        qp$colors[max(1, min(findInterval(v, qp$breaks, rightmost.closed = TRUE), 5))]
      })
      legend_labels <- c(qp$labels, "No Data")
      legend_colors <- c(qp$colors, "#cccccc")
    } else {
      pal <- colorNumeric("Blues", domain = range(vals, na.rm = TRUE), na.color = "#cccccc")
      geo$fill_color <- ifelse(is.na(geo[[var]]), "#cccccc", sapply(geo[[var]], pal))
      legend_labels <- NULL
      legend_colors <- NULL
    }
    ml <- bll_mean_label(use_geo())
    geo$tip <- paste0("<b>", geo$GEOID, "</b><br>", var_label, ": ",
                      ifelse(!is.na(geo[[var]]), round(geo[[var]], 2), "N/A"), "<br>",
                      "Mean BLL ", ml, ": ", ifelse(!is.na(geo$mean_bll), round(geo$mean_bll, 2), "N/A"), "<br>",
                      "% Elevated: ", ifelse(!is.na(geo$pct_elevated), round(geo$pct_elevated, 1), "N/A"))
    proxy <- leafletProxy("census_map") %>% clearShapes() %>% clearControls() %>%
      addPolygons(data = geo, fillColor = ~fill_color, fillOpacity = 0.7,
                  color = "white", weight = 0.5, layerId = ~GEOID,
                  label = lapply(geo$tip, HTML))
    if (input$census_quintile) {
      proxy %>% addLegend("bottomright", title = var_label, opacity = 0.8,
                          colors = legend_colors, labels = legend_labels)
    } else {
      proxy %>% addLegend("bottomright", title = var_label, opacity = 0.8,
                          pal = colorNumeric("Blues", range(vals, na.rm = TRUE)),
                          values = seq(min(vals), max(vals), length.out = 5))
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$census_map_shape_click, { rv$census_id <- input$census_map_shape_click$id })

  output$census_info_ui <- renderUI({
    if (is.null(rv$census_id)) return(p("Click an area for details", style = "color: #888;"))
    cd <- get_census_geo_data()
    row <- cd$combined %>% filter(geo_id == rv$census_id)
    if (nrow(row) == 0) return(p("No data"))
    ml <- bll_mean_label(use_geo())
    var <- input$census_var
    census_val <- if (var %in% names(row) && !is.na(row[[var]])) round(row[[var]], 2) else "N/A"
    tagList(
      div(style = "background: #f8f9fa; padding: 8px; border-radius: 4px; border-left: 3px solid #2c7fb8;",
        p(strong("Area:"), rv$census_id, style = "margin: 0 0 4px 0;"),
        p(strong(get_label(var), ":"), census_val, style = "margin: 0;")),
      div(style = "margin-top: 10px; padding: 8px; background: #fafafa; border-radius: 4px;",
        p(strong("Surveillance Metrics"), style = "margin: 0 0 4px 0; font-size: 12px; color: #2c7fb8;"),
        p("Tested:", format(row$n_children, big.mark = ","), " | Elevated:", format(row$n_elevated, big.mark = ","),
          style = "margin: 2px 0; font-size: 12px;"),
        p("% Elevated:", sprintf("%.2f%%", row$pct_elevated), style = "margin: 2px 0; font-size: 12px;"),
        p("Mean BLL", ml, ":", sprintf("%.2f", row$mean_bll), "\u00b5g/dL", style = "margin: 2px 0; font-size: 12px;")))
  })

  output$census_state_overall <- renderUI({
    dat <- active_data()
    s <- state_summary(dat, use_geo = use_geo())
    ml <- bll_mean_label(use_geo())
    var <- input$census_var
    census_df <- if (input$census_geo == "tract") tract_census else bg_census
    state_census_val <- if (!is.null(census_df) && var %in% names(census_df))
      round(mean(census_df[[var]], na.rm = TRUE), 1) else "N/A"
    div(class = "state-avg-box",
        p(strong("Tested:"), format(s$n_children, big.mark = ","),
          " | ", strong("Elevated:"), format(s$n_elevated, big.mark = ",")),
        p(strong("% Elevated:"), sprintf("%.2f%%", s$pct_elevated),
          " | ", strong("Avg BLL"), ml, ":", sprintf("%.2f", s$mean_bll)),
        p(strong(get_label(var), "(state avg):"), state_census_val))
  })

  get_quintile_data <- reactive({
    cd <- get_census_geo_data()
    var <- input$census_var; d <- cd$combined
    if (!var %in% names(d)) return(NULL)
    d <- d %>% filter(!is.na(.data[[var]]), !is.na(pct_elevated))
    if (nrow(d) < 10) return(NULL)
    d$quintile <- ntile(d[[var]], 5)
    d %>% group_by(quintile) %>%
      summarize(pct_elev = mean(pct_elevated, na.rm = TRUE),
                avg_bll = calc_mean_bll(mean_bll, use_geo()),
                n_areas = n(), tested = sum(n_children, na.rm = TRUE),
                elevated = sum(n_elevated, na.rm = TRUE),
                range_lo = round(min(.data[[var]], na.rm = TRUE), 1),
                range_hi = round(max(.data[[var]], na.rm = TRUE), 1),
                .groups = "drop")
  })

  output$quintile_pct_plot <- renderPlotly({
    qd <- get_quintile_data(); if (is.null(qd)) return(NULL)
    var_label <- get_label(input$census_var)
    max_y <- max(qd$pct_elev, na.rm = TRUE) * 1.25
    plot_ly(qd, x = ~paste0("Q", quintile), y = ~pct_elev, type = "bar",
            marker = list(color = "#2c7fb8"),
            text = ~paste0(round(pct_elev, 1), "%"),
            textposition = "outside", textfont = list(size = 11),
            hovertext = ~paste0("Q", quintile, " (", range_lo, "-", range_hi, ")",
                               "\n% Elevated: ", round(pct_elev, 2),
                               "\nAreas: ", n_areas, " | Tested: ", tested),
            hoverinfo = "text") %>%
      layout(title = list(text = paste("% Elevated by", var_label), font = list(size = 12)),
             xaxis = list(title = paste(var_label, "Quintile")),
             yaxis = list(visible = FALSE, range = c(0, max_y)),
             margin = list(l = 20, b = 70, t = 40))
  })

  output$quintile_bll_plot <- renderPlotly({
    qd <- get_quintile_data(); if (is.null(qd)) return(NULL)
    var_label <- get_label(input$census_var)
    ml <- bll_mean_label(use_geo())
    max_y <- max(qd$avg_bll, na.rm = TRUE) * 1.25
    plot_ly(qd, x = ~paste0("Q", quintile), y = ~avg_bll, type = "bar",
            marker = list(color = "#d73027"),
            text = ~paste0(round(avg_bll, 2), " \u00b5g/dL"),
            textposition = "outside", textfont = list(size = 10),
            hovertext = ~paste0("Q", quintile, " (", range_lo, "-", range_hi, ")",
                               "\nAvg BLL ", ml, ": ", round(avg_bll, 2), " \u00b5g/dL",
                               "\nTested: ", format(tested, big.mark = ","),
                               "\nAreas: ", format(n_areas, big.mark = ",")),
            hoverinfo = "text") %>%
      layout(title = list(text = paste("Avg BLL by", var_label), font = list(size = 12)),
             xaxis = list(title = paste(var_label, "Quintile")),
             yaxis = list(visible = FALSE, range = c(0, max_y)),
             margin = list(l = 20, b = 70, t = 40))
  })

  output$census_quintile_table <- renderDT({
    qd <- get_quintile_data(); if (is.null(qd)) return(NULL)
    tbl <- qd %>% transmute(Quintile = paste0("Q", quintile),
                             Range = paste0(range_lo, "-", range_hi),
                             Areas = n_areas, Tested = tested, Elevated = elevated,
                             `% Elevated (elevated/tested)` = round(pct_elev, 2),
                             `Avg BLL` = round(avg_bll, 2))
    datatable(tbl, options = list(dom = "t", pageLength = 5), rownames = FALSE)
  })

  set_dl("dl_quintile", function() {
    qd <- get_quintile_data(); if (is.null(qd)) return(data.frame())
    qd %>% transmute(Quintile = paste0("Q", quintile), Range = paste0(range_lo, "-", range_hi),
                     Areas = n_areas, Tested = tested, Elevated = elevated,
                     `% Elevated (elevated/tested)` = round(pct_elev, 2), `Avg BLL` = round(avg_bll, 2))
  }, "quintile_summary")

  output$census_area_table <- renderDT({
    cd <- get_census_geo_data()
    d <- cd$combined %>% filter(!is.na(n_children))
    ml <- bll_mean_label(use_geo())
    if (!is.null(geo_county_lookup))
      d <- d %>% left_join(geo_county_lookup %>% select(GEOID, County), by = c("geo_id" = "GEOID"))
    else d$County <- NA_character_
    d %>% transmute(County, Area = geo_id, Tested = n_children, Elevated = n_elevated,
                    `% Elevated` = round(pct_elevated, 2),
                    !!paste("Avg BLL", ml) := round(mean_bll, 2)) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
  })

  set_dl("dl_area_detail", function() {
    cd <- get_census_geo_data()
    d <- cd$combined %>% filter(!is.na(n_children))
    if (!is.null(geo_county_lookup))
      d <- d %>% left_join(geo_county_lookup %>% select(GEOID, County), by = c("geo_id" = "GEOID"))
    else d$County <- NA
    d %>% select(County, geo_id, n_children, n_elevated, pct_elevated, mean_bll) %>% as.data.frame()
  }, "area_details")

  # ===== DATA TABLES =====

  make_geo_dt <- function(geo_col, id_name, show_county = FALSE, show_lhd = FALSE) {
    yr <- input$table_years; unit <- input$table_unit
    dat <- active_data()
    ft <- dat %>% filter(sample_year >= yr[1], sample_year <= yr[2])
    d <- summarize_lead(ft, geo_col, input$suppress_n, use_geo())
    ml <- bll_mean_label(use_geo())
    if (show_county && !is.null(geo_county_lookup)) {
      d <- d %>% left_join(geo_county_lookup %>% select(GEOID, County), by = c("geo_id" = "GEOID"))
      d$county <- coalesce(d$County, d$county); d$County <- NULL
    }
    if (show_lhd) {
      if (!is.null(geo_lhd_lookup)) {
        lhd_lkp <- geo_lhd_lookup %>% select(GEOID, lhd = LHD) %>% distinct(GEOID, .keep_all = TRUE)
        d <- d %>% left_join(lhd_lkp, by = c("geo_id" = "GEOID"))
      } else if (!is.null(lhd_map_df) && "county" %in% names(d)) {
        d <- d %>% left_join(lhd_map_df, by = "county") %>% rename(lhd = lhd)
      }
    }
    tbl <- format_table(d, id_name, unit, show_county = show_county, show_lhd = show_lhd, mean_label = ml)
    list(table = tbl, data = d)
  }

  output$lhd_table <- renderDT({
    res <- make_geo_dt("lhd", "LHD")
    datatable(res$table, options = list(pageLength = as.integer(input$table_rows), scrollX = TRUE), rownames = FALSE)
  })
  output$county_table <- renderDT({
    res <- make_geo_dt("county", "County", show_lhd = TRUE)
    datatable(res$table, options = list(pageLength = as.integer(input$table_rows), scrollX = TRUE), rownames = FALSE)
  })
  output$geo_table <- renderDT({
    g_col <- if (input$table_geo == "tract") "tract_geoid" else "bg_geoid"
    id_name <- if (input$table_geo == "tract") "Tract" else "Block Group"
    res <- make_geo_dt(g_col, id_name, show_county = TRUE, show_lhd = TRUE)
    datatable(res$table, options = list(pageLength = as.integer(input$table_rows), scrollX = TRUE), rownames = FALSE)
  })

  output$addr_table <- renderDT({
    yr <- input$table_years; ml <- bll_mean_label(use_geo())
    astats <- get_address_stats() %>% filter(last_year >= yr[1], first_year <= yr[2])
    if (use_geo()) {
      dat <- active_data()
      geo_means <- dat %>% filter(sample_year >= yr[1], sample_year <= yr[2]) %>%
        group_by(address_id) %>%
        summarize(mean_bll_geo = calc_mean_bll(result, TRUE), .groups = "drop")
      astats <- astats %>% left_join(geo_means, by = "address_id") %>%
        mutate(mean_bll = coalesce(mean_bll_geo, mean_bll)) %>% select(-mean_bll_geo)
    }
    if (!is.null(rv$full_preds)) {
      astats <- astats %>% left_join(
        rv$full_preds %>% select(address_id, pred_prob = pred,
                                  single_child_flag) %>% distinct(address_id, .keep_all = TRUE),
        by = "address_id")
    } else { astats$pred_prob <- NA_real_; astats$single_child_flag <- FALSE }
    astats$single_child_flag[is.na(astats$single_child_flag)] <- FALSE
    bll_col <- paste("Avg BLL", ml)
    tbl <- astats %>%
      transmute(Address = street_address, LHD = lhd, County = county,
                Tested = n_children, Elevated = n_elevated,
                `% Elevated (elevated/tested)` = round(100 * n_elevated / pmax(n_children, 1), 2),
                !!bll_col := round(mean_bll, 2), `Max BLL` = round(max_bll, 2),
                `First Year` = first_year, `Last Year` = last_year,
                `Pred Risk` = if_else(!is.na(pred_prob),
                  paste0(sprintf("%.1f%%", pred_prob * 100), if_else(single_child_flag, "*", "")), "\u2014"),
                `View Tests` = paste0('<button class="btn btn-xs btn-info" ',
                                      'onclick="Shiny.setInputValue(\'view_tests_addr_click\', \'',
                                      address_id, '\', {priority: \'event\'})">View</button>'),
                `Risk Breakdown` = paste0('<button class="btn btn-xs btn-warning" ',
                                          'onclick="Shiny.setInputValue(\'risk_breakdown_click\', \'',
                                          address_id, '\', {priority: \'event\'})"',
                                          if_else(!is.na(pred_prob), '>', ' disabled title="Run model first">'),
                                          if_else(!is.na(pred_prob), 'Explain', '\u2014'), '</button>'))
    datatable(tbl, escape = FALSE, options = list(pageLength = as.integer(input$table_rows), scrollX = TRUE), rownames = FALSE,
              caption = if (any(astats$single_child_flag, na.rm = TRUE))
                htmltools::tags$caption(style = "font-size:11px; color:#888; font-style:italic;",
                  "* Single-child address: predicted risk shown but not included in model validation.") else NULL)
  })
  observeEvent(input$view_tests_addr_click, { show_address_tests(input$view_tests_addr_click) })

  # Risk Breakdown modal
  show_risk_breakdown <- function(addr_id) {
    fit <- rv$fit; preds <- rv$full_preds
    if (is.null(fit) || is.null(preds)) {
      showNotification("Please run the Risk Model first.", type = "warning"); return()
    }
    row <- preds %>% filter(address_id == addr_id)
    if (nrow(row) == 0) { showNotification("Address not found in model data.", type = "warning"); return() }
    row <- row[1, ]
    is_single_child <- isTRUE(row$single_child_flag)
    addr_label <- tryCatch({
      a <- get_address_stats() %>% filter(address_id == addr_id)
      if (nrow(a) > 0) a$street_address[1] else addr_id
    }, error = function(e) addr_id)
    use_re <- inherits(fit, "glmerMod")
    coefs <- if (use_re) fixef(fit) else coef(fit)
    intercept <- coefs["(Intercept)"]
    beta <- coefs[names(coefs) != "(Intercept)"]
    friendly_name <- function(v) {
      lkp <- c(pct_pre1980 = "% Pre-1980 Housing", pct_pre1950 = "% Pre-1950 Housing",
        pct_poverty = "Poverty Rate (%)", pct_rent_burden = "Rent Burden (%)",
        pct_renter = "% Renter-Occupied", median_income = "Median Income ($)",
        pct_nonwhite = "% Non-White", pct_vacant = "% Vacant Housing",
        log_housing_density = "Housing Density (log)", pct_children_u6 = "% Children Under 6",
        n_children = "Children Tested at Address", n_elevated = "Prior Elevated Count",
        first_bll = "First BLL (\u00b5g/dL)", mean_bll = "Average BLL", max_bll = "Max BLL",
        first_year = "First Test Year", last_year = "Last Test Year",
        income_cat = "Income Category", first_sex = "Child Sex", first_sample_type = "Sample Type",
        re_hispanic = "% Hispanic", re_white_nh = "% White (Non-Hispanic)",
        re_black_nh = "% Black (Non-Hispanic)", re_aian_nh = "% Am. Indian/AK Native",
        re_asian_nh = "% Asian (Non-Hispanic)", re_nhopi_nh = "% Native HI/Pac. Islander",
        re_multi_nh = "% Two+ Races", re_other_nh = "% Other Race")
      if (v %in% names(lkp)) return(lkp[v])
      gsub("_", " ", tools::toTitleCase(v))
    }
    contrib_rows <- list(); running_logodds <- intercept
    for (nm in names(beta)) {
      b <- beta[nm]; base_var <- nm; val <- NA_real_; val_display <- ""; factor_match <- FALSE
      for (fv in c("income_cat", "first_sex", "first_sample_type")) {
        if (startsWith(nm, fv) && nchar(nm) > nchar(fv)) {
          base_var <- fv; level_name <- substring(nm, nchar(fv) + 1)
          actual_val <- as.character(row[[fv]]); is_active <- (actual_val == level_name)
          val <- if (is_active) 1 else 0
          val_display <- paste0(actual_val, if (is_active) paste0(" (= ", level_name, ")")
                                else paste0(" (\u2260 ", level_name, ")"))
          factor_match <- TRUE; break
        }
      }
      if (!factor_match) {
        if (nm %in% names(row) && !is.na(row[[nm]])) {
          val <- as.numeric(row[[nm]]); val_display <- round(val, 3)
        } else { val <- 0; val_display <- "N/A" }
      }
      contribution <- b * val; running_logodds <- running_logodds + contribution
      contrib_rows[[length(contrib_rows) + 1]] <- list(
        name = friendly_name(base_var), value = val_display,
        contribution = round(contribution, 4),
        direction = if (contribution > 0.005) "\u2191 increases risk"
                    else if (contribution < -0.005) "\u2193 decreases risk"
                    else "\u2194 minimal effect")
    }
    contrib_rows <- contrib_rows[order(-sapply(contrib_rows, function(x) abs(x$contribution)))]
    re_text <- NULL
    if (use_re) {
      re_var <- names(ranef(fit))[1]; re_val <- row[[re_var]]
      if (!is.null(re_val) && !is.na(re_val)) {
        re_vals <- ranef(fit)[[re_var]]
        re_adj <- tryCatch(re_vals[as.character(re_val), 1], error = function(e) 0)
        if (!is.na(re_adj) && re_adj != 0) {
          running_logodds <- running_logodds + re_adj
          re_text <- div(style = "background: #f0f7ff; padding: 12px; border-radius: 6px; margin: 10px 0;",
            p(strong("\u2699 Neighborhood adjustment:"),
              sprintf("This address is in %s %s, which adjusts the score by %+.3f based on area-level patterns.",
                      re_var, re_val, re_adj), style = "margin: 0; font-size: 13px;"))
        }
      }
    }
    final_prob <- 1 / (1 + exp(-running_logodds)); final_pct <- round(final_prob * 100, 1)
    risk_color <- if (final_pct >= 30) "#d73027" else if (final_pct >= 15) "#fc8d59"
                  else if (final_pct >= 5) "#fee08b" else "#91cf60"
    risk_label <- if (final_pct >= 30) "High" else if (final_pct >= 15) "Moderate"
                  else if (final_pct >= 5) "Low-Moderate" else "Low"
    tbl_header <- '<table style="width:100%; border-collapse:collapse; font-size:13px; margin: 10px 0;">
      <tr style="background:#f5f5f5; font-weight:bold; border-bottom:2px solid #ddd;">
        <td style="padding:8px;">Factor</td><td style="padding:8px;">This Address</td>
        <td style="padding:8px; text-align:right;">Score Impact</td>
        <td style="padding:8px;">Direction</td></tr>'
    tbl_rows <- ""
    for (cr in contrib_rows) {
      bg <- if (cr$contribution > 0.005) "#fff5f5" else if (cr$contribution < -0.005) "#f0fff0" else "#fff"
      tbl_rows <- paste0(tbl_rows,
        '<tr style="border-bottom:1px solid #eee; background:', bg, ';">',
        '<td style="padding:7px;">', cr$name, '</td>',
        '<td style="padding:7px;">', cr$value, '</td>',
        '<td style="padding:7px; text-align:right; font-weight:bold;">', sprintf("%+.4f", cr$contribution), '</td>',
        '<td style="padding:7px; font-size:12px;">', cr$direction, '</td></tr>')
    }
    tbl_html <- paste0(tbl_header, tbl_rows, '</table>')
    showModal(modalDialog(
      title = paste0("Risk Breakdown: ", addr_label), size = "l", easyClose = TRUE,
      div(style = "background:#eef6ff; padding:14px; border-radius:8px; margin-bottom:14px;",
        p(strong("What is this?"), style = "margin:0 0 5px 0;"),
        p("This breakdown shows how the model estimated this address's risk of a future elevated blood lead level.",
          "The model looks at address characteristics and neighborhood data, weighs each factor, and combines",
          "them into a single risk score.", style = "margin:0; font-size:13px;")),
      if (is_single_child) {
        div(style = "background:#fff8e1; border:1px solid #ffc107; padding:12px; border-radius:6px; margin-bottom:14px;",
          p(HTML("<b>\u26a0 Single-child address:</b> This address has only one recorded child. Predicted risk is shown but was not included in model validation."),
            style = "margin:0; font-size:12px; color:#856404;"))
      } else NULL,
      div(style = paste0("text-align:center; padding:18px; border-radius:8px; border:2px solid ", risk_color, "; margin-bottom:16px;"),
        h3(paste0(final_pct, "% estimated risk"), style = paste0("color:", risk_color, "; margin:0 0 4px 0;")),
        p(HTML(paste0("Risk Level: <b>", risk_label, "</b>")), style = "margin:0; font-size:14px;"),
        p(em(paste0("Based on the model, roughly ", final_pct, " out of 100 similar addresses",
             " would be expected to have a future elevated blood lead test.")),
          style = "margin:6px 0 0 0; font-size:12px; color:#666;")),
      h5("How did the model get this number?", style = "color:#2c7fb8; margin-top:16px;"),
      p("The model starts with a baseline score of", strong(sprintf("%.3f", intercept)),
        "then adjusts based on each factor. Red rows push risk up; green rows push risk down. Sorted by impact size.",
        style = "font-size:13px; margin-bottom:8px;"),
      HTML(tbl_html),
      if (!is.null(re_text)) re_text else NULL,
      div(style = "background:#f9f9f9; padding:12px; border-radius:6px; margin:12px 0;",
        p(HTML(paste0("Combined score: <b>", sprintf("%.3f", running_logodds),
                      "</b> \u2192 Risk = 1 / (1 + e<sup>-(", sprintf("%.3f", running_logodds),
                      ")</sup>) = <b>", final_pct, "%</b>")),
          style = "font-size:13px; margin:0; text-align:center;")),
      br(),
      p(em("Note: This is a statistical estimate, not a certainty. The model uses historical data and",
           " neighborhood patterns to estimate risk. Always consider local knowledge and on-the-ground conditions."),
        style = "font-size:11px; color:#888;"),
      footer = modalButton("Close")))
  }
  observeEvent(input$risk_breakdown_click, { show_risk_breakdown(input$risk_breakdown_click) })

  # Download handlers
  set_dl("dl_lhd_disp", function() make_geo_dt("lhd", "LHD")$table, "lhd_data")
  set_dl("dl_lhd_all", function() make_geo_dt("lhd", "LHD")$table, "lhd_all")
  set_dl("dl_county_disp", function() make_geo_dt("county", "County", show_lhd = TRUE)$table, "county_data")
  set_dl("dl_county_all", function() make_geo_dt("county", "County", show_lhd = TRUE)$table, "county_all")
  set_dl("dl_geo_disp", function() {
    g <- if (input$table_geo == "tract") "tract_geoid" else "bg_geoid"
    n <- if (input$table_geo == "tract") "Tract" else "Block Group"
    make_geo_dt(g, n, show_county = TRUE, show_lhd = TRUE)$table
  }, "geo_data")
  set_dl("dl_geo_all", function() {
    g <- if (input$table_geo == "tract") "tract_geoid" else "bg_geoid"
    n <- if (input$table_geo == "tract") "Tract" else "Block Group"
    make_geo_dt(g, n, show_county = TRUE, show_lhd = TRUE)$table
  }, "geo_all")
  set_dl("dl_addr_disp", function() {
    get_address_stats() %>%
      filter(last_year >= input$table_years[1], first_year <= input$table_years[2]) %>%
      select(street_address, lhd, county, n_children, n_elevated, mean_bll, max_bll, first_year, last_year)
  }, "addresses")
  set_dl("dl_addr_all", function() {
    get_address_stats() %>%
      select(street_address, lhd, county, n_children, n_elevated, mean_bll, max_bll, first_year, last_year)
  }, "addresses_all")

  # ===== RISK MODEL =====

  do_fit_model <- function() {
    showNotification("Fitting model...", id = "model_note", type = "message", duration = NULL)
    on.exit(removeNotification("model_note"))

    re <- input$random_effect
    use_re <- re != "none"
    md <- if (re == "bg_geoid") model_data_bg else model_data_tract
    if (is.null(md)) { showNotification("Model data not loaded", type = "error"); return(NULL) }

    all_covars <- get_all_covars()
    if (length(all_covars) == 0) { showNotification("Select at least one variable", type = "error"); return(NULL) }

    factor_covars  <- intersect(all_covars, c("income_cat", "first_sex", "first_sample_type"))
    re_covars      <- intersect(all_covars, paste0("re_", c("hispanic", "white_nh", "black_nh", "aian_nh",
                                                             "asian_nh", "nhopi_nh", "multi_nh", "other_nh")))
    numeric_covars <- setdiff(all_covars, c(factor_covars, re_covars))

    d <- md
    if (length(numeric_covars) > 0) d <- d %>% filter(if_all(all_of(numeric_covars), ~ !is.na(.) & is.finite(.)))
    if (length(factor_covars) > 0)  d <- d %>% filter(if_all(all_of(factor_covars), ~ !is.na(.)))
    if (length(re_covars) > 0)
      d <- d %>% filter(rowSums(!is.na(across(all_of(re_covars)))) > 0)

    train_yrs <- input$train_years
    if (input$train_all) {
      train <- d; test <- NULL
    } else {
      train <- d %>% filter(first_year >= train_yrs[1], first_year <= train_yrs[2])
      test  <- d %>% filter(first_year > train_yrs[2])
    }
    if (nrow(train) < 50) { showNotification("Too few training records", type = "error"); return(NULL) }

    rv$train_sample <- train; rv$test_sample <- test

    if (isTRUE(input$use_log_bll) && "first_bll" %in% all_covars) {
      log_transform <- function(df) { df$first_bll <- log(pmax(df$first_bll, 0.1)); df }
      train <- log_transform(train)
      if (!is.null(test)) test <- log_transform(test)
      d <- log_transform(d)
    }

    scaled_vars <- character(0); scale_params <- list()
    never_scale <- c(factor_covars, re_covars)
    if (input$scale_vars && length(numeric_covars) > 0) {
      to_scale <- setdiff(numeric_covars, never_scale)
      for (v in to_scale) {
        m <- mean(train[[v]], na.rm = TRUE); s <- sd(train[[v]], na.rm = TRUE)
        if (!is.na(s) && s > 0) {
          scale_params[[v]] <- list(mean = m, sd = s)
          train[[v]] <- (train[[v]] - m) / s
          if (!is.null(test)) test[[v]] <- (test[[v]] - m) / s
          d[[v]] <- (d[[v]] - m) / s
          scaled_vars <- c(scaled_vars, v)
        }
      }
    }
    rv$is_scaled <- length(scaled_vars) > 0

    covar_str <- paste(all_covars, collapse = " + ")
    if (use_re) {
      f <- as.formula(paste0("outcome ~ ", covar_str, " + (1 | ", re, ")"))
    } else {
      f <- as.formula(paste0("outcome ~ ", covar_str))
    }

    fit <- suppressWarnings(tryCatch({
      if (use_re) {
        glmer(f, data = train, family = binomial, nAGQ = 1,
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))
      } else {
        glm(f, data = train, family = binomial)
      }
    }, error = function(e) { showNotification(paste("Model error:", e$message), type = "error"); NULL }))
    if (is.null(fit)) return(NULL)

    if (use_re && inherits(fit, "glmerMod") && isSingular(fit))
      showNotification("Warning: Model is singular (random effect variance near zero).", type = "warning", duration = 8)

    rv$fit <- fit; rv$model_run <- TRUE

    if (!is.null(test) && nrow(test) > 0) {
      test$pred <- predict(fit, newdata = test, type = "response",
                           allow.new.levels = if (use_re) TRUE else NULL)
      rv$test_preds <- test
      rv$roc_obj <- compute_roc(test$outcome, test$pred)
    }
    d$pred <- predict(fit, newdata = d, type = "response",
                      allow.new.levels = if (use_re) TRUE else NULL)
    d$single_child_flag <- FALSE
    rv$full_preds <- d

    # Predict for single-child addresses
    tryCatch({
      all_addr <- address_stats_all
      if (!is.null(all_addr) && nrow(all_addr) > 0) {
        existing_ids <- unique(d$address_id)
        single_addr <- all_addr %>% filter(!address_id %in% existing_ids, n_children == 1)
        if (nrow(single_addr) > 0) {
          census_join_col <- if (re == "bg_geoid") "bg_geoid" else "tract_geoid"
          census_src <- if (re == "bg_geoid") bg_census else tract_census
          single_df <- if (!is.null(census_src)) {
            single_addr %>% left_join(census_src, by = setNames(census_join_col, census_join_col))
          } else single_addr
          missing_covars <- setdiff(all_covars, names(single_df))
          for (mc in missing_covars) single_df[[mc]] <- NA_real_
          if (length(numeric_covars) > 0)
            single_df <- single_df %>% filter(if_all(all_of(numeric_covars), ~ !is.na(.) & is.finite(.)))
          if (length(factor_covars) > 0)
            single_df <- single_df %>% filter(if_all(all_of(factor_covars), ~ !is.na(.)))
          if (nrow(single_df) > 0) {
            if (isTRUE(input$use_log_bll) && "first_bll" %in% all_covars)
              single_df$first_bll <- log(pmax(single_df$first_bll, 0.1))
            if (input$scale_vars && length(scaled_vars) > 0) {
              for (v in scaled_vars) {
                if (v %in% names(scale_params) && v %in% names(single_df)) {
                  m <- scale_params[[v]]$mean; s <- scale_params[[v]]$sd
                  if (!is.na(s) && s > 0) single_df[[v]] <- (single_df[[v]] - m) / s
                }
              }
            }
            if (use_re && !re %in% names(single_df)) single_df[[re]] <- NA_character_
            single_df$outcome <- NA_integer_
            single_df$pred <- predict(fit, newdata = single_df, type = "response", allow.new.levels = TRUE)
            single_df$single_child_flag <- TRUE
            shared_cols <- intersect(names(rv$full_preds), names(single_df))
            rv$full_preds <- bind_rows(rv$full_preds[shared_cols], single_df[shared_cols])
          }
        }
      }
    }, error = function(e) { message("[Single-child prediction] ", e$message) })

    all_numeric <- c(numeric_covars, re_covars)
    if (length(all_numeric) >= 2) {
      rv$cor_matrix_pearson <- cor(train[all_numeric], use = "pairwise.complete.obs", method = "pearson")
      rv$cor_matrix_spearman <- cor(train[all_numeric], use = "pairwise.complete.obs", method = "spearman")
      rv$cor_matrix <- rv$cor_matrix_pearson
    } else { rv$cor_matrix <- NULL; rv$cor_matrix_pearson <- NULL; rv$cor_matrix_spearman <- NULL }

    # Build comparison stats
    if (!is.null(test) && nrow(test) > 0) {
      n_test <- nrow(test); n_elev <- sum(test$outcome == 1)
      screen_n <- round(n_test * 0.5)
      test_ranked <- test %>% arrange(desc(pred))
      model_det <- sum(test_ranked$outcome[1:screen_n] == 1)
      random_det <- round(n_elev * 0.5)
      rv$comparison <- list(
        n_test = n_test, n_elev = n_elev, screen_n = screen_n,
        model_det = model_det, random_det = random_det,
        model_pct = round(100 * model_det / max(n_elev, 1), 1), random_pct = 50,
        improve = round(100 * (model_det - random_det) / max(random_det, 1), 1),
        auc = round(rv$roc_obj$auc, 3),
        tests_per_pos_model = round(screen_n / max(model_det, 1), 1),
        tests_per_pos_random = round(screen_n / max(random_det, 1), 1))
    }

    # XGBoost comparison
    showNotification("Fitting XGBoost comparison...", id = "xgb_note", type = "message", duration = NULL)
    rv$xgb_fit <- NULL; rv$xgb_roc <- NULL; rv$xgb_comparison <- NULL; rv$xgb_test_preds <- NULL
    tryCatch({
      x_formula <- as.formula(paste0("~ ", covar_str))
      x_train <- model.matrix(x_formula, data = train)[, -1, drop = FALSE]
      y_train <- train$outcome
      dtrain <- xgb.DMatrix(data = x_train, label = y_train)
      if (!is.null(test) && nrow(test) > 0) {
        x_test <- model.matrix(x_formula, data = test)[, -1, drop = FALSE]
        dtest <- xgb.DMatrix(data = x_test, label = test$outcome)
        watchlist <- list(train = dtrain, test = dtest)
      } else { watchlist <- list(train = dtrain); dtest <- NULL }
      params <- list(objective = "binary:logistic", eval_metric = "auc",
                     max_depth = 4, eta = 0.1, subsample = 0.8)
      xgb_fit <- xgb.train(params = params, data = dtrain, nrounds = 200,
                           watchlist = watchlist, early_stopping_rounds = 20, verbose = 0)
      rv$xgb_fit <- xgb_fit
      if (!is.null(test) && nrow(test) > 0) {
        xgb_preds <- predict(xgb_fit, dtest)
        rv$xgb_test_preds <- xgb_preds
        rv$xgb_roc <- compute_roc(test$outcome, xgb_preds)
        xgb_ord <- order(xgb_preds, decreasing = TRUE)
        n_half <- ceiling(nrow(test) / 2)
        xgb_capture <- sum(test$outcome[xgb_ord[1:n_half]]) / sum(test$outcome)
        rv$xgb_comparison <- list(
          glmer_auc = rv$roc_obj$auc, xgb_auc = rv$xgb_roc$auc,
          glmer_capture = rv$comparison$model_pct / 100, xgb_capture = xgb_capture,
          n_test = nrow(test),
          importance = xgb.importance(model = xgb_fit) %>% head(10))
      }
    }, error = function(e) { showNotification(paste("XGBoost warning:", e$message), type = "warning") })
    removeNotification("xgb_note")
    showNotification("Model fitting complete!", type = "message", duration = 3)
    return(fit)
  }

  observeEvent(input$fit_model, {
    tryCatch(do_fit_model(), error = function(e) {
      showNotification(paste("Model error:", e$message), type = "error", duration = 10)
    })
  }, ignoreInit = TRUE)

  observe({
    invalidateLater(2000, session)
    isolate({
      if (!rv$model_run) {
        rv$model_run <- TRUE
        tryCatch(do_fit_model(), error = function(e) {
          showNotification(paste("Auto-run error:", e$message), type = "warning")
        })
      }
    })
  })

  output$interpretation <- renderUI({
    fit <- rv$fit
    if (is.null(fit)) return(div(class = "metric-box", p("Click 'Run Model' to fit.")))
    is_glmer <- inherits(fit, "glmerMod")
    coefs <- tryCatch({
      if (is_glmer) {
        tidy(fit, effects = "fixed", conf.int = TRUE) %>%
          filter(term != "(Intercept)") %>%
          mutate(or = exp(estimate), var_name = gsub("scale\\(|\\)", "", term))
      } else {
        tidy(fit, conf.int = TRUE) %>%
          filter(term != "(Intercept)") %>%
          mutate(or = exp(estimate), var_name = gsub("scale\\(|\\)", "", term))
      }
    }, error = function(e) data.frame())
    is_scaled <- rv$is_scaled
    unit_text <- if (is_scaled) "1 standard deviation" else "1 unit"
    findings <- if (nrow(coefs) > 0) {
      lapply(1:nrow(coefs), function(i) {
        r <- coefs[i, ]; pct_change <- round((r$or - 1) * 100, 0); abs_pct <- abs(pct_change)
        tags$li(strong(get_label(r$var_name)), ": ",
          "A ", unit_text, " ", if (pct_change > 0) "increase" else "decrease", " is associated with ",
          abs_pct, "% ", if (r$or > 1) "higher" else "lower", " odds of elevated BLL ",
          "(OR = ", round(r$or, 2), ")")
      })
    } else { list(tags$li("Model uses geographic grouping only.")) }
    auc_val <- rv$comparison$auc
    auc_text <- if (!is.null(auc_val)) {
      tagList(p(strong("GLMER AUC = ", auc_val)),
        p("The AUC measures how well the model discriminates between addresses that will vs. won't have elevated children.
           Values range from 0.5 (no better than chance) to 1.0 (perfect). AUC > 0.7 is acceptable; > 0.8 is good."),
        p(if (auc_val >= 0.8) "This model shows good discrimination."
          else if (auc_val >= 0.7) "This model shows acceptable discrimination."
          else "This model shows modest discrimination. Consider adding more predictors."))
    }
    xgb_text <- if (!is.null(rv$xgb_roc)) {
      xgb_auc <- round(rv$xgb_roc$auc, 3)
      diff_auc <- round(xgb_auc - (auc_val %||% 0.5), 3)
      comp_word <- if (diff_auc > 0.01) "outperforms" else if (diff_auc < -0.01) "underperforms" else "matches"
      tagList(p(strong("XGBoost AUC = ", xgb_auc)),
        p(paste0("XGBoost ", comp_word, " GLMER by ", abs(diff_auc), " AUC points.")))
    }
    div(class = "metric-box",
        h4("What the Model Found:"), tags$ul(findings),
        if (is_scaled) p(em("Note: Variables were standardized. ORs represent 1 SD change, not 1-unit change."),
                         style = "font-size: 11px; color: #666;"),
        hr(), h4("Model Performance:"), auc_text, xgb_text,
        hr(), h4("About XGBoost"),
        p("XGBoost is a machine learning algorithm that builds an ensemble of decision trees sequentially. ",
          "Unlike the mixed-effects logistic regression (GLMER), XGBoost can automatically capture non-linear ",
          "relationships and complex interactions between variables."),
        p("XGBoost is always run alongside the GLMER for benchmarking. If XGBoost substantially outperforms ",
          "GLMER, it suggests the data contains non-linear patterns that a logistic model misses."),
        p(em("XGBoost results are shown in the Model Comparison tab."), style = "font-size: 11px; color: #666;"))
  })

  output$or_note <- renderUI({
    if (is.null(rv$fit)) return(NULL)
    div(class = "warning-box",
        if (rv$is_scaled) {
          "These odds ratios are STANDARDIZED. Each OR represents the change in odds for a 1 standard deviation
           increase in the predictor."
        } else {
          "These odds ratios are UNSTANDARDIZED. Each OR represents the change in odds for a 1-unit increase
           in the predictor. Magnitudes are not directly comparable across variables with different scales."
        })
  })

  output$coef_plot <- renderPlotly({
    fit <- rv$fit; if (is.null(fit)) return(NULL)
    is_glmer <- inherits(fit, "glmerMod")
    coefs <- tryCatch({
      raw <- if (is_glmer) tidy(fit, effects = "fixed", conf.int = TRUE) else tidy(fit, conf.int = TRUE)
      raw %>% filter(term != "(Intercept)") %>%
        mutate(or = exp(estimate), or_lo = exp(conf.low), or_hi = exp(conf.high),
               var_name = gsub("scale\\(|\\)", "", term),
               label = sapply(var_name, get_label),
               Effect = ifelse(or > 1.5, "Strong positive",
                               ifelse(or > 1.1, "Moderate positive",
                                      ifelse(or < 0.67, "Strong negative",
                                             ifelse(or < 0.9, "Moderate negative", "Weak/null")))))
    }, error = function(e) data.frame())
    if (nrow(coefs) == 0) return(NULL)
    p <- ggplot(coefs, aes(or, reorder(label, or),
                           text = paste0(label, "\nOR: ", round(or, 2),
                                         "\n95% CI: ", round(or_lo, 2), "-", round(or_hi, 2)))) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_errorbar(aes(xmin = or_lo, xmax = or_hi), width = 0.2, color = "gray40", orientation = "y") +
      geom_point(aes(color = Effect), size = 4) +
      scale_color_manual(values = c("Strong positive" = "#d73027", "Moderate positive" = "#fc8d59",
                                    "Weak/null" = "#999999",
                                    "Moderate negative" = "#91bfdb", "Strong negative" = "#4575b4")) +
      labs(title = "Odds Ratios (with 95% CI)", x = "Odds Ratio", y = NULL) +
      theme_minimal() + theme(legend.position = "bottom")
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.2))
  })

  output$cor_table <- renderDT({
    cm <- if (isTRUE(input$cor_method == "spearman")) rv$cor_matrix_spearman else rv$cor_matrix_pearson
    if (is.null(cm)) cm <- rv$cor_matrix
    if (is.null(cm)) return(NULL)
    cor_df <- as.data.frame(round(cm, 3))
    cor_df$Variable <- sapply(rownames(cor_df), get_label)
    cor_df <- cor_df %>% select(Variable, everything())
    names(cor_df)[-1] <- sapply(names(cor_df)[-1], get_label)
    datatable(cor_df, options = list(dom = "t", pageLength = 15, scrollX = TRUE), rownames = FALSE) %>%
      formatStyle(columns = 2:ncol(cor_df),
                  backgroundColor = styleInterval(c(-0.7, -0.5, 0.5, 0.7, 0.99),
                                                  c('#f8d7da', '#fff3cd', 'white', '#fff3cd', '#f8d7da', 'white')))
  })

  output$vif_table <- renderDT({
    cm <- if (isTRUE(input$cor_method == "spearman")) rv$cor_matrix_spearman else rv$cor_matrix_pearson
    if (is.null(cm)) cm <- rv$cor_matrix
    if (is.null(cm) || ncol(cm) < 2) return(NULL)
    tryCatch({
      inv_cm <- solve(cm); vifs <- diag(inv_cm)
      vif_df <- data.frame(Variable = sapply(names(vifs), get_label),
                           VIF = round(vifs, 2), stringsAsFactors = FALSE) %>% arrange(desc(VIF))
      datatable(vif_df, options = list(dom = "t", pageLength = 15), rownames = FALSE) %>%
        formatStyle("VIF", backgroundColor = styleInterval(c(5, 10), c("white", "#fff3cd", "#f8d7da")))
    }, error = function(e) NULL)
  })

  # Model Comparison
  output$model_comparison_ui <- renderUI({
    comp <- rv$xgb_comparison; rcomp <- rv$comparison
    if (is.null(comp)) return(p("Model comparison will appear after fitting."))
    glmer_auc <- round(comp$glmer_auc, 3); xgb_auc <- round(comp$xgb_auc, 3)
    auc_diff <- round(xgb_auc - glmer_auc, 3)
    glmer_cap <- comp$glmer_capture * 100; xgb_cap <- comp$xgb_capture * 100
    xgb_screen_n <- rcomp$screen_n; xgb_det <- NA_integer_; xgb_tpp <- NA_real_
    if (!is.null(rv$xgb_test_preds) && !is.null(rv$test_preds)) {
      xgb_ord <- order(rv$xgb_test_preds, decreasing = TRUE)
      xgb_det <- sum(rv$test_preds$outcome[xgb_ord[1:xgb_screen_n]] == 1)
      xgb_tpp <- round(xgb_screen_n / max(xgb_det, 1), 1)
    }
    tagList(
      if (!is.null(rcomp)) {
        div(class = "compare-box",
            h4("Screening Efficiency: All Models at 50% Screen Rate", style = "margin-top: 0;"),
            p("If you could only test 50% of addresses (",
              format(rcomp$screen_n, big.mark = ","), " of ", format(rcomp$n_test, big.mark = ","), "):"),
            hr(),
            fluidRow(
              column(4, div(style = "text-align: center; padding: 15px; background: #e8f4f8; border-radius: 6px;",
                  h5("GLMER", style = "color: #2c7fb8; margin: 0;"),
                  h2(format(rcomp$model_det, big.mark = ","), style = "margin: 10px 0; color: #2c7fb8;"),
                  p("elevated addresses found"), p(strong(rcomp$model_pct, "%"), " of all elevated"),
                  hr(), p(strong(rcomp$tests_per_pos_model), " tests per positive"))),
              column(4, div(style = "text-align: center; padding: 15px; background: #fff5f5; border-radius: 6px;",
                  h5("XGBoost", style = "color: #d73027; margin: 0;"),
                  h2(if (!is.na(xgb_det)) format(xgb_det, big.mark = ",") else "N/A",
                     style = "margin: 10px 0; color: #d73027;"),
                  p("elevated addresses found"), p(strong(sprintf("%.1f%%", xgb_cap)), " of all elevated"),
                  hr(), p(strong(if (!is.na(xgb_tpp)) xgb_tpp else "N/A"), " tests per positive"))),
              column(4, div(style = "text-align: center; padding: 15px; background: #f5f5f5; border-radius: 6px;",
                  h5("Random Selection", style = "color: #666; margin: 0;"),
                  h2(format(rcomp$random_det, big.mark = ","), style = "margin: 10px 0; color: #666;"),
                  p("elevated addresses found"), p(strong("50%"), " of all elevated"),
                  hr(), p(strong(rcomp$tests_per_pos_random), " tests per positive")))),
            hr(),
            div(style = "text-align: center;",
                h4(paste0("GLMER: +", rcomp$improve, "% over random | XGBoost: +",
                           round((xgb_cap - 50) / 50 * 100, 1), "% over random"),
                   style = "color: #4a9; margin: 0;"),
                p(em("Lower tests per positive = more efficient screening"))))
      },
      hr(),
      fluidRow(
        column(4, div(style = "background:#f0f7fb; padding:14px; border-radius:6px; border-left:4px solid #2c7fb8;",
          h4("GLMER", style = "margin-top:0; color:#2c7fb8;"),
          p(strong("AUC:"), glmer_auc), p(strong("Capture @ 50%:"), sprintf("%.1f%%", glmer_cap)))),
        column(4, div(style = "background:#fff5f5; padding:14px; border-radius:6px; border-left:4px solid #d73027;",
          h4("XGBoost", style = "margin-top:0; color:#d73027;"),
          p(strong("AUC:"), xgb_auc), p(strong("Capture @ 50%:"), sprintf("%.1f%%", xgb_cap)),
          p(em(if (abs(auc_diff) < 0.02) "Similar to GLMER"
               else if (auc_diff > 0) paste0("+", auc_diff, " vs GLMER")
               else paste0(auc_diff, " vs GLMER")), style = "font-size:11px; color:#555;"))),
        column(4, div(style = "background:#f5f5f5; padding:14px; border-radius:6px; border-left:4px solid #999;",
          h4("Random Baseline", style = "margin-top:0; color:#999;"),
          p(strong("AUC:"), "0.500"), p(strong("Capture @ 50%:"), "50.0%")))),
      if (!is.null(comp$importance) && nrow(comp$importance) > 0) {
        tagList(hr(), h5("XGBoost Top Features (by Gain)"), DTOutput("xgb_importance_table"))
      })
  })
  output$xgb_importance_table <- renderDT({
    comp <- rv$xgb_comparison
    if (is.null(comp) || is.null(comp$importance)) return(NULL)
    datatable(comp$importance %>% mutate(across(where(is.numeric), ~ round(., 4))),
              options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })

  output$capture_overlay_plot <- renderPlotly({
    test <- rv$test_preds
    if (is.null(test) || nrow(test) == 0) return(NULL)
    n <- nrow(test); n_elev <- sum(test$outcome)
    if (n_elev == 0) return(NULL)
    glmer_ord <- order(test$pred, decreasing = TRUE)
    glmer_cum <- cumsum(test$outcome[glmer_ord]) / n_elev
    idx <- unique(c(1, round(seq(1, n, length.out = 500)), n))
    pct_scr <- idx / n * 100
    p <- plot_ly() %>%
      add_trace(x = pct_scr, y = glmer_cum[idx] * 100, type = "scatter", mode = "lines",
                name = "GLMER", line = list(color = "#2c7fb8", width = 3),
                hovertemplate = "GLMER<br>Screen: %{x:.1f}%<br>Capture: %{y:.1f}%<extra></extra>")
    if (!is.null(rv$xgb_test_preds)) {
      xgb_ord <- order(rv$xgb_test_preds, decreasing = TRUE)
      xgb_cum <- cumsum(test$outcome[xgb_ord]) / n_elev
      p <- p %>% add_trace(x = pct_scr, y = xgb_cum[idx] * 100, type = "scatter", mode = "lines",
                            name = "XGBoost", line = list(color = "#d73027", width = 3),
                            hovertemplate = "XGBoost<br>Screen: %{x:.1f}%<br>Capture: %{y:.1f}%<extra></extra>")
    }
    p %>% add_trace(x = c(0, 100), y = c(0, 100), type = "scatter", mode = "lines",
                      name = "Random", line = list(color = "gray60", width = 2, dash = "dash")) %>%
      layout(title = list(text = "Cumulative Capture Curves", font = list(size = 15)),
             xaxis = list(title = "% Addresses Screened", ticksuffix = "%", range = c(0, 100)),
             yaxis = list(title = "% Elevated Captured", ticksuffix = "%", range = c(0, 100)),
             legend = list(orientation = "h", y = -0.15), margin = list(t = 50, b = 60))
  })

  output$auc_year_overlay_plot <- renderPlotly({
    test <- rv$test_preds
    if (is.null(test) || nrow(test) == 0) return(NULL)
    test_years <- sort(unique(test$first_year))
    glmer_aucs <- sapply(test_years, function(yr) {
      sub <- test %>% filter(first_year == yr)
      if (sum(sub$outcome) < 3 || sum(!sub$outcome) < 3) return(NA)
      compute_roc(sub$outcome, sub$pred)$auc
    })
    valid <- !is.na(glmer_aucs)
    if (sum(valid) < 2) return(NULL)
    p <- plot_ly() %>%
      add_trace(x = test_years[valid], y = glmer_aucs[valid], type = "scatter", mode = "lines+markers",
                name = "GLMER", line = list(color = "#2c7fb8", width = 2.5),
                marker = list(color = "#2c7fb8", size = 7))
    if (!is.null(rv$xgb_test_preds)) {
      xgb_aucs <- sapply(test_years, function(yr) {
        idx <- which(test$first_year == yr)
        if (sum(test$outcome[idx]) < 3) return(NA)
        compute_roc(test$outcome[idx], rv$xgb_test_preds[idx])$auc
      })
      xvalid <- !is.na(xgb_aucs)
      if (sum(xvalid) >= 2)
        p <- p %>% add_trace(x = test_years[xvalid], y = xgb_aucs[xvalid], type = "scatter", mode = "lines+markers",
                              name = "XGBoost", line = list(color = "#d73027", width = 2.5),
                              marker = list(color = "#d73027", size = 7))
    }
    all_aucs <- c(glmer_aucs[valid])
    p %>% add_trace(x = range(test_years), y = c(0.5, 0.5), type = "scatter", mode = "lines",
                     line = list(color = "gray60", width = 1, dash = "dash"), showlegend = FALSE, hoverinfo = "skip") %>%
      layout(title = list(text = "AUC by Test Year", font = list(size = 14)),
             xaxis = list(title = "Year"),
             yaxis = list(title = "AUC", range = c(max(0.4, min(all_aucs) - 0.05), 1)),
             legend = list(orientation = "h", y = -0.2), margin = list(t = 40))
  })

  output$cal_time_overlay_plot <- renderPlotly({
    test <- rv$test_preds
    if (is.null(test) || nrow(test) == 0) return(NULL)
    test_years <- sort(unique(test$first_year))
    glmer_pred <- sapply(test_years, function(yr) mean(test$pred[test$first_year == yr]))
    observed <- sapply(test_years, function(yr) mean(test$outcome[test$first_year == yr]))
    p <- plot_ly() %>%
      add_trace(x = test_years, y = observed * 100, type = "scatter", mode = "lines+markers",
                name = "Observed", line = list(color = "#333", width = 2.5),
                marker = list(color = "#333", size = 7)) %>%
      add_trace(x = test_years, y = glmer_pred * 100, type = "scatter", mode = "lines+markers",
                name = "GLMER predicted", line = list(color = "#2c7fb8", width = 2.5, dash = "dash"),
                marker = list(color = "#2c7fb8", size = 7))
    if (!is.null(rv$xgb_test_preds)) {
      xgb_pred <- sapply(test_years, function(yr) {
        idx <- which(test$first_year == yr); mean(rv$xgb_test_preds[idx])
      })
      p <- p %>% add_trace(x = test_years, y = xgb_pred * 100, type = "scatter", mode = "lines+markers",
                            name = "XGBoost predicted", line = list(color = "#d73027", width = 2.5, dash = "dash"),
                            marker = list(color = "#d73027", size = 7))
    }
    p %>% layout(title = list(text = "Calibration Over Time", font = list(size = 14)),
                 xaxis = list(title = "Year"),
                 yaxis = list(title = "Elevated Rate (%)", ticksuffix = "%"),
                 legend = list(orientation = "h", y = -0.2), margin = list(t = 40))
  })

  # Sample tab plots
  get_model_data_for_plot <- reactive({
    re <- input$random_effect
    if (re == "bg_geoid") model_data_bg else model_data_tract
  })

  output$sample_plot_addr <- renderPlotly({
    if (is.null(addr_by_year) || nrow(addr_by_year) == 0) return(NULL)
    p <- ggplot(addr_by_year, aes(sample_year, n_addr, fill = type,
                                  text = paste0("Year: ", sample_year, "\n", type, ": ",
                                               format(n_addr, big.mark = ",")))) +
      geom_col(position = "stack") +
      scale_fill_manual(values = c("2+ Children" = "#2c7fb8", "1 Child" = "#b0d5e8")) +
      scale_y_continuous(labels = comma) +
      labs(title = "Unique Addresses Over Time (by Child Count)", x = "Year", y = "Addresses", fill = "") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggplotly(p, tooltip = "text") %>%
      layout(autosize = TRUE, xaxis = list(autorange = TRUE), yaxis = list(autorange = TRUE))
  })

  output$sample_plot_tested <- renderPlotly({
    md <- get_model_data_for_plot()
    if (is.null(md)) return(NULL)
    d <- md %>% group_by(sample_year = first_year) %>%
      summarize(Elevated = sum(outcome == 1, na.rm = TRUE),
                `Not Elevated` = sum(outcome == 0, na.rm = TRUE), .groups = "drop")
    plot_ly(d, x = ~sample_year) %>%
      add_bars(y = ~Elevated, name = "Elevated", marker = list(color = "#d73027"),
               hovertext = ~paste0("Year: ", sample_year, "\nElevated: ", format(Elevated, big.mark = ",")),
               hoverinfo = "text", textposition = "none") %>%
      add_bars(y = ~`Not Elevated`, name = "Not Elevated", marker = list(color = "#b0d5e8"),
               hovertext = ~paste0("Year: ", sample_year, "\nNot Elevated: ", format(`Not Elevated`, big.mark = ",")),
               hoverinfo = "text", textposition = "none") %>%
      layout(barmode = "stack",
             title = list(text = "Addresses: Elevated vs Not", font = list(size = 13)),
             xaxis = list(title = "First Test Year", autorange = TRUE),
             yaxis = list(title = "Count", autorange = TRUE),
             legend = list(orientation = "h", y = -0.2), autosize = TRUE)
  })

  output$sample_plot_bll <- renderPlotly({
    md <- get_model_data_for_plot()
    if (is.null(md)) return(NULL)
    d <- md %>% group_by(sample_year = first_year) %>%
      summarize(bll = mean(first_bll, na.rm = TRUE), .groups = "drop")
    train_end <- if (!input$train_all) input$train_years[2] else max(d$sample_year)
    d$split <- if_else(d$sample_year <= train_end, "Train", "Test")
    use_log <- isTRUE(input$sample_bll_log)
    p <- ggplot(d, aes(sample_year, bll, group = 1,
                       text = paste0("Year: ", sample_year, "\nMean BLL: ", round(bll, 2), "\n", split))) +
      geom_line(linewidth = 0.5, color = "#555") + geom_point(aes(color = split), size = 2) +
      scale_color_manual(values = c("Train" = "#2c7fb8", "Test" = "#d73027")) +
      labs(title = "Mean First Child BLL by Year", x = "First Test Year", y = "\u00b5g/dL", color = "") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    if (use_log) p <- p + scale_y_log10()
    ggplotly(p, tooltip = "text") %>%
      layout(autosize = TRUE, xaxis = list(autorange = TRUE), yaxis = list(autorange = TRUE))
  })

  output$sample_plot_elev <- renderPlotly({
    md <- get_model_data_for_plot()
    if (is.null(md)) return(NULL)
    d <- md %>% group_by(sample_year = first_year) %>%
      summarize(pct = 100 * mean(outcome, na.rm = TRUE), .groups = "drop")
    train_end <- if (!input$train_all) input$train_years[2] else max(d$sample_year)
    d$split <- if_else(d$sample_year <= train_end, "Train", "Test")
    p <- ggplot(d, aes(sample_year, pct, group = 1,
                       text = paste0("Year: ", sample_year, "\n% Elevated: ", round(pct, 1), "%\n", split))) +
      geom_line(linewidth = 0.5, color = "#555") + geom_point(aes(color = split), size = 2) +
      scale_color_manual(values = c("Train" = "#2c7fb8", "Test" = "#d73027")) +
      labs(title = "% Elevated (Subsequent Child)", x = "First Test Year", y = "%", color = "") +
      theme_minimal(base_size = 11) + theme(legend.position = "bottom")
    ggplotly(p, tooltip = "text") %>%
      layout(autosize = TRUE, xaxis = list(autorange = TRUE), yaxis = list(autorange = TRUE))
  })

  output$dl_sample_table <- downloadHandler(
    filename = function() paste0("sample_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      md <- if (input$random_effect == "bg_geoid") model_data_bg else model_data_tract
      if (is.null(md)) { write.csv(data.frame(), file); return() }
      d <- md %>% group_by(Year = first_year) %>%
        summarize(Addresses = n(), `% Elevated` = round(100 * mean(outcome, na.rm = TRUE), 2),
                  `Mean First BLL` = round(mean(first_bll, na.rm = TRUE), 2), .groups = "drop")
      write.csv(d, file, row.names = FALSE)
    })

  output$model_summary <- renderPrint({
    if (is.null(rv$fit)) return(cat("No model fitted"))
    summary(rv$fit)
  })

  # ===== PREDICTION MAP OBSERVER =====
  observeEvent(c(input$pred_type, input$pred_geo, input$pred_class, input$pred_palette, rv$model_run), {
    preds <- if (input$pred_type == "test") rv$test_preds else rv$full_preds
    if (is.null(preds) || nrow(preds) == 0) return()
    geo_col <- switch(input$pred_geo, county = "county", lhd = "lhd",
                      tract = "tract_geoid", bg = "bg_geoid", "tract_geoid")
    geo_sf <- switch(input$pred_geo, county = ne_counties, lhd = lhd_bounds,
                     tract = ne_tracts, bg = ne_bgs, ne_tracts)
    geo_id_col <- switch(input$pred_geo, county = "NAME", lhd = "lhd",
                         tract = "GEOID", bg = "GEOID", "GEOID")
    if (!geo_col %in% names(preds) || is.null(geo_sf)) return()
    pred_col <- if ("pred" %in% names(preds)) "pred" else if ("pred_prob" %in% names(preds)) "pred_prob" else return()
    agg <- preds %>% filter(!is.na(.data[[geo_col]])) %>%
      group_by(geo_id = .data[[geo_col]]) %>%
      summarize(mean_pred = mean(.data[[pred_col]], na.rm = TRUE), n_addr = n(),
                pct_elevated = 100 * mean(outcome, na.rm = TRUE), .groups = "drop")
    geo <- geo_sf %>% left_join(agg, by = setNames("geo_id", geo_id_col))
    vals <- geo$mean_pred[!is.na(geo$mean_pred)]
    if (length(vals) == 0) return()

    class_mode <- input$pred_class %||% "quintile"
    use_viridis <- isTRUE(input$pred_palette == "viridis")
    use_cividis <- isTRUE(input$pred_palette == "cividis")

    if (class_mode %in% c("quintile", "decile")) {
      n_bins <- if (class_mode == "quintile") 5 else 10
      qbreaks <- quantile(vals, probs = seq(0, 1, 1/n_bins), na.rm = TRUE)
      risk_colors <- if (use_viridis) viridis::viridis(n_bins)
                     else if (use_cividis) viridis::cividis(n_bins)
                     else colorRampPalette(c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026"))(n_bins)
      geo$bin <- findInterval(geo$mean_pred, qbreaks, rightmost.closed = TRUE)
      geo$bin <- pmin(pmax(geo$bin, 1), n_bins)
      geo$fill <- ifelse(is.na(geo$mean_pred), "#e0e0e0", risk_colors[geo$bin])
      prefix <- if (class_mode == "quintile") "Q" else "D"
      leg_labels <- paste0(prefix, 1:n_bins, ": ", sprintf("%.1f%%", qbreaks[1:n_bins] * 100), " \u2013 ",
                           sprintf("%.1f%%", qbreaks[2:(n_bins+1)] * 100))
    } else {
      max_val <- max(0.20, max(vals, na.rm = TRUE))
      risk_pal <- if (use_viridis) {
        colorNumeric(viridis::viridis(100), domain = c(0, max_val), na.color = "#e0e0e0")
      } else if (use_cividis) {
        colorNumeric(viridis::cividis(100), domain = c(0, max_val), na.color = "#e0e0e0")
      } else {
        colorNumeric(c("#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c",
                        "#fc4e2a", "#e31a1c", "#bd0026", "#800026"),
                      domain = c(0, max_val), na.color = "#e0e0e0")
      }
      geo$fill <- ifelse(is.na(geo$mean_pred), "#e0e0e0", sapply(geo$mean_pred, risk_pal))
    }

    geo_label_col <- geo[[geo_id_col]]
    geo$risk_cat <- case_when(
      is.na(geo$mean_pred) ~ "No Data",
      geo$mean_pred >= 0.15 ~ "\u26a0 High Risk",
      geo$mean_pred >= 0.08 ~ "Moderate Risk",
      geo$mean_pred >= 0.04 ~ "Low-Moderate",
      TRUE ~ "Lower Risk")
    geo$tip <- paste0(
      "<div style='font-size:13px; line-height:1.6;'>",
      "<b style='font-size:14px;'>", geo_label_col, "</b><br>",
      "<span style='font-size:18px; font-weight:bold; color:",
      ifelse(is.na(geo$mean_pred), "#999",
             ifelse(geo$mean_pred >= 0.15, "#bd0026",
                    ifelse(geo$mean_pred >= 0.08, "#fc4e2a", "#333"))), ";'>",
      ifelse(!is.na(geo$mean_pred), sprintf("%.1f%%", geo$mean_pred * 100), "N/A"),
      "</span> predicted risk<br>",
      "<b>", geo$risk_cat, "</b><br>",
      "Observed: ", ifelse(!is.na(geo$pct_elevated), sprintf("%.1f%%", geo$pct_elevated), "N/A"), " elevated<br>",
      "Addresses: ", ifelse(!is.na(geo$n_addr), format(geo$n_addr, big.mark = ","), "0"), "</div>")

    proxy <- leafletProxy("pred_map") %>% clearShapes() %>% clearControls() %>%
      addPolygons(data = geo, fillColor = ~fill, fillOpacity = 0.85,
                  color = "#444", weight = 0.3, opacity = 0.5,
                  highlightOptions = highlightOptions(weight = 2.5, color = "#000", fillOpacity = 0.95),
                  label = lapply(geo$tip, HTML))

    if (class_mode %in% c("quintile", "decile")) {
      proxy %>% addLegend("bottomright", colors = c(risk_colors, "#e0e0e0"),
                           labels = c(leg_labels, "No Data"),
                           title = HTML("<b>Predicted Risk</b>"), opacity = 0.9)
    } else {
      proxy %>% addLegend("bottomright", pal = risk_pal,
                           values = seq(0, max_val, length.out = 6),
                           title = HTML("<b>Predicted Risk</b>"),
                           labFormat = labelFormat(suffix = "%", transform = function(x) round(x * 100, 1)),
                           opacity = 0.9)
    }
  })

  output$pred_map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron, group = "Light") %>%
      addProviderTiles(providers$CartoDB.DarkMatter, group = "Dark") %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Satellite") %>%
      addLayersControl(baseGroups = c("Light", "Dark", "Satellite"), position = "topright") %>%
      setView(-99.5, 41.5, zoom = 7)
  })
}

# ==============================================================================
# RUN APP
# ==============================================================================

shinyApp(ui, server)
