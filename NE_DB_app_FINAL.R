###############################################################################
# NE_DB_app_FINAL.R — Nebraska Childhood Lead Surveillance Dashboard
#
# FINAL consolidated version combining all features from prior iterations:
#   - Mixed-effects logistic regression (GLMER), logistic (GLM), and XGBoost
#   - ACS modes: year-matched, non-overlapping fixed (2013/2018/2023), endpoint
#   - First test vs max BLL child selection for surveillance and modeling
#   - All addresses vs 2+ children population toggle
#   - Arithmetic vs geometric mean BLL
#   - Temporal outcome recomputation to prevent data leakage
#   - Pre-aggregated geo stats for fast map rendering
#   - Lazy-loaded model data to conserve memory
#   - Race/ethnicity covariates for modeling
#   - ACS under-6 population denominators for rate calculations
#   - Cumulative capture curves, calibration, confusion matrices
#   - Single-child address predictions in 2+ children mode
#
# Requires: NE_DB_prep_optimized.R output in K:/CLPPP/TrentonS/shiny_data
#
# Data leakage prevention:
#   1. Outcomes are recomputed at runtime with temporal boundary — training
#      labels never include information from future test-period children
#   2. The feature child (first_id) is excluded from their own outcome
#   3. Max-BLL features are rebuilt within temporal windows when using
#      max BLL child mode with temporal splits
#   4. train_all mode warns that evaluation metrics are not valid
###############################################################################
# 1. CTRL + A → 2. CTRL + Enter → 3. Wait ~60s for default model

packages <- c("shiny", "shinydashboard", "leaflet", "sf", "tidyverse",
              "DT", "lme4", "broom.mixed", "plotly", "xgboost", "scales")
invisible(lapply(packages, library, character.only = TRUE))

'%||%' <- function(a, b) if (!is.null(a)) a else b
DATA_PATH <- "K:/CLPPP/TrentonS/shiny_data"; BLL_CUTOFF <- 3.5

# ==============================================================================
# HELPERS
# ==============================================================================

re_flag <- function(re, pat, ex = TRUE) {
  h <- grepl(pat, re) & !is.na(re); if (ex) h <- h & !grepl("hispanic|latino", re); as.integer(h)
}

rebuild_addr_features <- function(ft_subset, mode = "first") {
  sliced <- if (mode == "max") {
    ft_subset |> slice_max(result, by = AddressIdentifier, with_ties = FALSE)
  } else {
    ft_subset |> slice_min(sample_year, by = AddressIdentifier, with_ties = FALSE)
  }
  sliced |>
    mutate(re = tolower(trimws(RACE)),
           re_hispanic = re_flag(re, "hispanic|latino", ex = FALSE),
           re_white_nh = re_flag(re, "white"), re_black_nh = re_flag(re, "black|african"),
           re_aian_nh = re_flag(re, "indian|alaska|aian|native am"),
           re_asian_nh = re_flag(re, "asian"), re_nhopi_nh = re_flag(re, "pacific|hawaiian|nhopi"),
           re_multi_nh = re_flag(re, "multi|two|more"), re_other_nh = re_flag(re, "other"),
           across(starts_with("re_"), \(x) if_else(is.na(RACE), NA_integer_, x))) |>
    select(AddressIdentifier, first_id = PATIENT_LOCAL_ID, first_year = sample_year,
           first_bll = result, first_age_mo = AGE_MO, first_sex = patient_SEX,
           first_sample_type = sample_type, bg_geoid, tract_geoid, county, lhd,
           lat = Latitude, lng = Longitude, starts_with("re_"))
}

calc_mean_bll <- function(x, use_geo = FALSE) {
  x <- x[!is.na(x)]; if (length(x) == 0) return(NA_real_)
  if (use_geo) exp(mean(log(pmax(x, 0.1)))) else mean(x)
}

bll_mean_label <- function(geo) if (geo) "(Geometric)" else "(Arithmetic)"

state_summary <- function(data, by_year = FALSE, by_address = FALSE, use_geo = FALSE) {
  if (by_address) {
    grp <- if (by_year) group_by(data, sample_year, AddressIdentifier) else group_by(data, AddressIdentifier)
    addr <- grp |> summarize(e = max(elevated, na.rm = TRUE),
                              bll = calc_mean_bll(result, use_geo), .groups = "drop")
    if (by_year) addr <- addr |> group_by(sample_year)
    addr |> summarize(mean_bll = mean(bll, na.rm = TRUE), n_children = n(),
                       n_elevated = sum(e > 0, na.rm = TRUE),
                       pct_elevated = 100 * mean(e > 0, na.rm = TRUE), .groups = "drop")
  } else {
    if (by_year) data <- data |> group_by(sample_year)
    data |> summarize(mean_bll = calc_mean_bll(result, use_geo), n_children = n(),
                       n_elevated = sum(elevated, na.rm = TRUE),
                       pct_elevated = 100 * mean(elevated, na.rm = TRUE), .groups = "drop")
  }
}

format_table <- function(df, id_name, unit = "Children", show_county = FALSE,
                         show_lhd = FALSE, mean_label = "") {
  tested   <- if (unit == "Addresses") df$n_addresses else df$n_children
  elevated <- if (unit == "Addresses") df$addr_elevated else df$n_elevated
  bll_col <- if (nchar(mean_label) > 0) paste("Average BLL", mean_label) else "Average BLL"
  out <- tibble(!!id_name := df$geo_id, Tested = tested, Elevated = elevated,
                `% Elevated (elevated/tested)` = round(100 * elevated / pmax(tested, 1), 2),
                !!bll_col := round(df$mean_bll, 2), Suppressed = df$suppressed)
  if ("n_acs_u6" %in% names(df)) {
    out$`ACS Under-6 Pop` <- round(df$n_acs_u6, 0)
    out$`Elevated/1k` <- round(df$rate_per_1k, 2)
    out$`Tested/1k` <- round(df$tested_rate_per_1k, 2)
  }
  if (show_lhd && "lhd" %in% names(df)) out <- bind_cols(tibble(LHD = df$lhd), out)
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

compute_trend_ci <- function(ft_yr, var, sv) {
  n_t <- nrow(ft_yr); if (n_t == 0 || is.na(sv)) return(c(NA_real_, NA_real_))
  switch(var,
    pct_elevated = { p <- sv/100; se <- sqrt(p*(1-p)/n_t)*100; c(max(0, sv-1.96*se), sv+1.96*se) },
    mean_bll = { se <- sd(ft_yr$result, na.rm=TRUE)/sqrt(n_t); c(sv-1.96*se, sv+1.96*se) },
    n_elevated = { p <- mean(ft_yr$elevated, na.rm=TRUE); se <- sqrt(n_t*p*(1-p)); c(max(0,sv-1.96*se),sv+1.96*se) },
    c(NA_real_, NA_real_))
}

build_quintile_palette <- function(vals, colors = c("#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026")) {
  q <- quantile(vals, probs = seq(0, 1, 0.2), na.rm = TRUE)
  list(breaks = q, colors = colors, labels = paste0("Q", 1:5, ": ", round(q[1:5], 1), "-", round(q[2:6], 1)))
}

assign_fill_colors <- function(geo, display_var, color_func) {
  v <- geo[[display_var]]
  sup <- if (!is.null(geo$suppressed)) !is.na(geo$suppressed) & geo$suppressed else rep(FALSE, length(v))
  ifelse(sup, "#666666", ifelse(is.na(v) | !is.finite(v), "#cccccc", sapply(v, color_func)))
}

base_map <- function() leaflet() |> addProviderTiles(providers$CartoDB.Positron) |> setView(-99.5, 41.5, 7)

make_dt <- function(df, pg = 10, ...) datatable(df, options = list(pageLength = pg, scrollX = TRUE, ...), rownames = FALSE)

get_risk_palette <- function(pal_name, n) {
  switch(pal_name, viridis = viridis::viridis(n), cividis = viridis::cividis(n),
         colorRampPalette(c("#ffffb2","#fecc5c","#fd8d3c","#f03b20","#bd0026"))(n))
}

add_model_traces <- function(p, rv, x_vals, y_multi, y_xgb = NULL, y_glm = NULL,
                              mode = "lines+markers", marker_size = 7) {
  cfgs <- list(list(y=y_multi, name="Multilevel", color="#2c7fb8", dash="solid"),
               list(y=y_xgb, name="XGBoost", color="#d73027", dash="solid"),
               list(y=y_glm, name="Logistic", color="#1a9850", dash="dot"))
  for (cfg in cfgs) {
    if (is.null(cfg$y)) next; valid <- !is.na(cfg$y); if (sum(valid) < 2) next
    p <- p |> add_trace(x=x_vals[valid], y=cfg$y[valid], type="scatter", mode=mode, name=cfg$name,
                         line=list(color=cfg$color, width=2.5, dash=cfg$dash),
                         marker=list(color=cfg$color, size=marker_size))
  }; p
}

get_bll_color <- function(bll) {
  case_when(bll >= 20 ~ "#67000d", bll >= 10 ~ "#a50f15", bll >= 5 ~ "#d73027",
            bll >= 3.5 ~ "#fc8d59", TRUE ~ "#c9a84c")
}

LABELS <- c(
  pct_pre1980="Pre-1980 Housing (%)", pct_pre1950="Pre-1950 Housing (%)",
  pct_poverty="Poverty Rate (%)", pct_rent_burden="Rent Burden \u226530% (%)",
  pct_renter="Renter-Occupied (%)", median_income="Median Income ($)",
  median_house_value="Median House Value ($)", income_cat="Income Category",
  pct_nonwhite="Non-White (%)", race_cat="Race Category",
  log_housing_density="Housing Density (log)", pct_vacant="Vacant Housing (%)",
  first_bll="First Child's BLL", first_year="First Test Year",
  first_year_cat="First Test Year (cat)", first_age_yr_cat="Child Age (years, cat)",
  first_sex="First Child Sex", first_age_mo="First Child Age (months)",
  first_sample_type="Sample Type (Capillary/Venous)",
  re_hispanic="Hispanic", re_white_nh="White (Non-Hispanic)",
  re_black_nh="Black (Non-Hispanic)", re_aian_nh="AIAN (Non-Hispanic)",
  re_asian_nh="Asian (Non-Hispanic)", re_nhopi_nh="NHOPI (Non-Hispanic)",
  re_multi_nh="Two or More Races (Non-Hispanic)", re_other_nh="Other Race (Non-Hispanic)",
  n_elevated="Elevated Count", n_children="Tested Count",
  pct_elevated="% Elevated (elevated/tested)", mean_bll="Average BLL (\u00b5g/dL)",
  n_addresses="Addresses Tested", addr_elevated="Addresses Elevated",
  rate_per_1k="Elevated per 1,000 Under-6", tested_rate_per_1k="Tested per 1,000 Under-6")
get_label <- function(v) if (v %in% names(LABELS)) unname(LABELS[v]) else gsub("_"," ",tools::toTitleCase(v))

# Geo-level lookup tables
GEO_COL <- c(LHD="lhd", County="county", Tract="tract_geoid", `Block Group`="bg_geoid")
GEO_ID  <- c(LHD="lhd", County="NAME", Tract="GEOID", `Block Group`="GEOID")
GEO_LBL <- c(LHD="LHD", County="County", Tract="Census Tract", `Block Group`="Block Group")
GEO_KEY <- c(LHD="lhd", County="county", Tract="tract", `Block Group`="bg")

# ==============================================================================
# LOAD DATA (core only — model data lazy-loaded)
# ==============================================================================

load_rds <- function(name) {
  path <- file.path(DATA_PATH, paste0(name, ".rds"))
  if (file.exists(path)) readRDS(path) else NULL
}

first_test     <- load_rds("first_test")
max_test       <- load_rds("max_test")
address_stats  <- load_rds("address_stats")
address_stats_all <- load_rds("address_stats_all")
ne_counties    <- load_rds("ne_counties")
lhd_bounds     <- load_rds("lhd_bounds")
ne_bgs         <- load_rds("ne_bgs")
ne_tracts      <- load_rds("ne_tracts")
bg_census      <- load_rds("bg_census")
tract_census   <- load_rds("tract_census")
lhd_map_df     <- load_rds("lhd_map")
all_tests      <- load_rds("all_tests")
tract_census_ts <- load_rds("tract_census_ts")
bg_census_ts   <- load_rds("bg_census_ts")
acs_u6_ts      <- load_rds("acs_u6_ts")
addr_by_year   <- load_rds("addr_by_year")
geo_county_lookup <- load_rds("geo_county_lookup")
geo_lhd_lookup <- load_rds("geo_lhd_lookup")

# Pre-aggregated geo stats (new from optimized prep)
geo_stats_first <- load_rds("geo_stats_first")
geo_stats_max   <- load_rds("geo_stats_max")

required_data <- list(first_test=first_test, address_stats_all=address_stats_all, ne_counties=ne_counties)
missing_req <- names(required_data)[sapply(required_data, is.null)]
if (length(missing_req) > 0) stop("Missing: ", paste(missing_req, collapse=", "), "\nRun NE_DB_prep first.")

sums <- load_rds("startup_summaries")
state_avgs <- sums$state_avgs; state_avg_all <- sums$state_avg_all
state_avgs_addr <- sums$state_avgs_addr; state_avg_all_addr <- sums$state_avg_all_addr
addr_elevated_summary <- sums$addr_elevated_summary
addr_elevated_summary_all <- sums$addr_elevated_summary_all %||% sums$addr_elevated_summary
years <- sums$years %||% sort(unique(first_test$sample_year))
lhd_names <- sums$lhd_names %||% character(0)

acs_year_map <- function(yr) pmax(2010L, pmin(as.integer(yr), 2023L))

# Lazy-load cache for model data
.md_cache <- new.env(parent = emptyenv())
get_model_data <- function(name) {
  if (!exists(name, envir = .md_cache)) .md_cache[[name]] <- load_rds(name)
  .md_cache[[name]]
}

resolve_model_data <- function(acs_mode, geo_key, use_max, use_all_addr) {
  base <- if (use_all_addr && use_max) "model_data_all_max"
          else if (use_all_addr) "model_data_all"
          else if (use_max) "model_data_max"
          else "model_data"
  if (acs_mode == "endpoint") return(NULL)
  if (acs_mode == "base") return(get_model_data(base))
  suffix <- paste0("_", geo_key, "_", acs_mode)
  get_model_data(paste0(base, suffix))
}

join_census_for_mode <- function(df, acs_mode, geo_col, census_ts, endpoint_acs_year = NULL) {
  if (is.null(census_ts)) return(NULL)
  if (acs_mode == "ym") {
    df |> mutate(acs_year = acs_year_map(first_year)) |>
      left_join(census_ts, by = c(geo_col, "acs_year")) |> select(-acs_year)
  } else if (acs_mode == "nof") {
    nof_yrs <- c(2013L, 2018L, 2023L)
    df |> mutate(acs_year = sapply(first_year, \(y) nof_yrs[which.min(abs(as.integer(y) - nof_yrs))])) |>
      left_join(census_ts, by = c(geo_col, "acs_year")) |> select(-acs_year)
  } else if (acs_mode == "endpoint") {
    ep_acs <- endpoint_acs_year %||% acs_year_map(max(years))
    df |> left_join(census_ts |> filter(acs_year == ep_acs) |> select(-acs_year), by = geo_col)
  } else {
    df |> mutate(acs_year = acs_year_map(first_year)) |>
      left_join(census_ts, by = c(geo_col, "acs_year")) |> select(-acs_year)
  }
}

# Shared: attach ACS under-6 rates
attach_rates <- function(stat_data, geo_level, yr_range) {
  if (is.null(acs_u6_ts)) return(stat_data)
  level_key <- GEO_KEY[geo_level] %||% "tract"
  u6_src <- acs_u6_ts[[level_key]]
  if (is.null(u6_src) || nrow(u6_src) == 0) return(stat_data)
  acs_years <- acs_year_map(yr_range[1]):acs_year_map(yr_range[2])
  avail <- unique(u6_src$acs_year)
  matched <- intersect(acs_years, avail)
  if (length(matched) == 0) matched <- avail[which.min(abs(avail - mean(acs_years)))]
  u6_agg <- u6_src |> filter(acs_year %in% matched) |>
    summarize(n_acs_u6 = sum(n_children_u6, na.rm = TRUE),
              n_acs_u6_avg = mean(n_children_u6, na.rm = TRUE), .by = geo_id)
  stat_data |> left_join(u6_agg, by = "geo_id") |>
    mutate(rate_per_1k = if_else(!is.na(n_acs_u6) & n_acs_u6 > 0, n_elevated/n_acs_u6*1000, NA_real_),
           tested_rate_per_1k = if_else(!is.na(n_acs_u6) & n_acs_u6 > 0, n_children/n_acs_u6*1000, NA_real_))
}

# Shared: summarize_lead for tables (still needed for filtered address-level recomputation)
summarize_lead <- function(data, geo_col, sup = 0, use_geo = FALSE) {
  data <- data |> filter(!is.na(.data[[geo_col]]))
  child <- data |> group_by(geo_id = .data[[geo_col]]) |>
    summarize(n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
              mean_bll = calc_mean_bll(result, use_geo), county = first(county), .groups = "drop")
  addr <- data |> group_by(geo_id = .data[[geo_col]], AddressIdentifier) |>
    summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") |>
    group_by(geo_id) |> summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")
  left_join(child, addr, by = "geo_id") |>
    mutate(pct_elevated = 100 * n_elevated / n_children, suppressed = n_children < sup)
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
      menuItem("Risk Model",         tabName = "model",    icon = icon("chart-line"))),
    hr(),
    numericInput("suppress_n", "Suppress areas with fewer than:", 10, min = 0, max = 50, step = 5),
    radioButtons("addr_pop", "Address Population:",
                 c("All Tested (1+ child)" = "all", "2+ Children Tested" = "multi"), inline = TRUE),
    p(em("Affects address counts. Switch 'Count By' to 'Addresses' on maps."),
      style = "padding: 0 15px; font-size: 10px; color: #aaa;"),
    radioButtons("test_selection", "Child Test Selection:",
                 c("First Test per Child" = "first", "Maximum BLL per Child" = "max"), inline = FALSE),
    p(em("First = earliest test. Maximum = highest BLL recorded."),
      style = "padding: 0 15px; font-size: 10px; color: #aaa;"),
    radioButtons("use_geo_mean", "BLL Average Type:",
                 c("Arithmetic Mean" = "arith", "Geometric Mean" = "geo"), inline = TRUE),
    hr(),
    sliderInput("map_opacity", "Map Fill Opacity:", min = 0, max = 100, value = 70, step = 5, post = "%")),

  dashboardBody(
    tags$head(tags$style(HTML("
      .skin-blue .main-header .logo { background-color: #00607F; }
      .skin-blue .main-header .navbar { background-color: #00607F; }
      .skin-blue .main-header .logo:hover { background-color: #004d66; }
      .skin-blue .main-sidebar { background-color: #333; }
      .section-title { font-size: 20px; font-weight: 600; margin-bottom: 6px; color: #00607F; }
      .section-desc { color: #444; margin-bottom: 15px; font-size: 13.5px; line-height: 1.6; }
      .metric-box { background: #f8f9fa; border-left: 3px solid #00607F; padding: 14px; margin: 8px 0; border-radius: 4px; }
      .state-avg-box { background: #eef5f8; border: 1px solid #00607F; padding: 10px; border-radius: 4px; margin-top: 10px; font-size: 13px; }
      .compare-box { background: #f8f9f5; border: 1px solid #BABF33; padding: 15px; border-radius: 6px; margin: 10px 0; }
      .warning-box { background: #FFF8E1; border: 1px solid #FFC843; padding: 10px; border-radius: 4px; margin: 5px 0; font-size: 12px; }
      .cor-high { background-color: #f8d7da !important; }
      .cor-borderline { background-color: #FFF8E1 !important; }
      .content-wrapper { font-size: 13.5px; }
      .box-title { font-size: 15px; color: #00607F; }
      .box.box-primary > .box-header { background-color: #00607F; }
      .box.box-primary { border-top-color: #00607F; }
      .btn-primary { background-color: #00607F; border-color: #004d66; }
      .btn-primary:hover { background-color: #004d66; border-color: #003d52; }
      .btn-info { background-color: #B9C8D3; border-color: #a0b4c2; color: #333; }
      .btn-warning { background-color: #FFC843; border-color: #e6b23a; color: #333; }
      .nav-tabs > li.active > a { border-top: 2px solid #00607F; }
      .tab-content .row > .col-sm-3 { width: 29.5% !important; }
      .tab-content .row > .col-sm-9 { width: 70.5% !important; }
    ")),
    tags$script(HTML("
      $(document).on('click', '.act-view-tests', function(e) {
        e.stopPropagation(); Shiny.setInputValue('view_tests_click', $(this).attr('data-addr'), {priority: 'event'});
      });
      $(document).on('click', '.act-view-tests-addr', function(e) {
        e.stopPropagation(); Shiny.setInputValue('view_tests_addr_click', $(this).attr('data-addr'), {priority: 'event'});
      });
      $(document).on('click', '.act-risk-breakdown', function(e) {
        e.stopPropagation(); Shiny.setInputValue('risk_breakdown_click', $(this).attr('data-addr'), {priority: 'event'});
      });
    "))),

    tabItems(
      # ===== OVERVIEW =====
      tabItem(tabName = "overview",
        h3(class = "section-title", "Nebraska Childhood Lead Surveillance Dashboard"),
        p(class = "section-desc",
          "Interactive tools for the Nebraska Lead Team to analyze childhood BLL testing data (2010-2023)."),
        p(class = "section-desc",
          strong("Surveillance Map:"), " View metrics by LHD, county, tract, or block group.", br(),
          strong("Priority Addresses:"), " Rank addresses with elevated children.", br(),
          strong("Census Data:"), " Explore neighborhood characteristics associated with lead risk.", br(),
          strong("Data Tables:"), " Browse and download complete datasets.", br(),
          strong("Risk Model:"), " Predictive modeling to identify high-risk addresses."),
        p(em("Elevated BLL defined as \u22653.5 \u00b5g/dL per 2021 CDC reference value. Data includes first test per child, ages 0-6.")),
        hr(),
        fluidRow(
          valueBox(format(nrow(first_test), big.mark=","), "Children Tested", icon=icon("child"), color="light-blue", width=6),
          valueBox(paste0(format(sum(first_test$elevated, na.rm=TRUE), big.mark=","), " (",
                          round(100*mean(first_test$elevated, na.rm=TRUE),2), "%)"),
                   "Children Elevated", icon=icon("exclamation-triangle"), color="blue", width=6)),
        fluidRow(
          valueBox(format(nrow(address_stats_all), big.mark=","), "Addresses Tested", icon=icon("home"), color="green", width=6),
          valueBox(paste0(format(sum(addr_elevated_summary_all$e > 0, na.rm=TRUE), big.mark=","),
                          " (", round(100*mean(addr_elevated_summary_all$e > 0, na.rm=TRUE),2), "%)"),
                   "Addresses with Elevated Child", icon=icon("home"), color="olive", width=6)),
        fluidRow(
          valueBox(format(nrow(address_stats), big.mark=","), "Addresses with 2+ Children Tested",
                   icon=icon("users"), color="yellow", width=6),
          valueBox(paste0(format(sum(address_stats$n_elevated > 0, na.rm=TRUE), big.mark=","), " (",
                          round(100*mean(address_stats$n_elevated > 0, na.rm=TRUE),2), "%)"),
                   "2+ Child Addresses with Elevated", icon=icon("users"), color="orange", width=6)),
        fluidRow(
          valueBox({
            if (!is.null(acs_u6_ts$state)) {
              pop <- acs_u6_ts$state$n_children_u6[acs_u6_ts$state$acs_year == 2023]
              if (length(pop) > 0 && pop > 0) paste0(round(sum(first_test$elevated, na.rm=TRUE)/pop*1000,1)," per 1,000")
              else "N/A"
            } else "N/A"
          }, "Cumulative Elevated Rate (per 1,000 Under-6, 2023 ACS)", icon=icon("chart-line"), color="purple", width=6),
          valueBox({
            if (!is.null(acs_u6_ts$state)) {
              pop <- acs_u6_ts$state$n_children_u6[acs_u6_ts$state$acs_year == 2023]
              if (length(pop) > 0 && pop > 0) format(pop, big.mark=",") else "N/A"
            } else "N/A"
          }, "ACS Under-6 Population (2023)", icon=icon("users"), color="teal", width=6)),
        p(em("Cumulative Elevated Rate = total elevated children (all years) / 2023 ACS Under-6 pop \u00d7 1,000."),
          style = "font-size: 11px; color: #666; padding: 0 15px;")),

      # ===== SURVEILLANCE MAP =====
      tabItem(tabName = "map",
        h3(class = "section-title", "Geographic Surveillance"),
        p(class = "section-desc", "Lead testing metrics aggregated to geographic levels. Click any region for details."),
        fluidRow(
          box(width = 9,
              leafletOutput("surv_map", height = 450), hr(),
              fluidRow(column(12, plotlyOutput("state_avg_plot", height = 220))),
              hr(), h5("Annual Summary Statistics"), DTOutput("trend_summary_table")),
          box(width = 3, title = "Map Options",
              sliderInput("year_filter", "Years:", min(years), max(years), c(min(years), max(years)), sep = ""),
              selectInput("geo_level", "Geographic Level:",
                          c("LHD"="LHD","County"="County","Census Tract"="Tract","Block Group"="Block Group")),
              selectInput("fill_var", "Measure:",
                          c("Average BLL (\u00b5g/dL)"="mean_bll","Elevated Count"="n_elevated",
                            "% Elevated (elevated/tested)"="pct_elevated","Tested Count"="n_children",
                            "Elevated per 1,000 Under-6"="rate_per_1k","Tested per 1,000 Under-6"="tested_rate_per_1k")),
              radioButtons("unit_type", "Count By:", c("Children","Addresses"), inline = TRUE),
              hr(),
              checkboxInput("use_quintile", "Display by Quintile", FALSE),
              checkboxInput("show_map_ci", "Show 95% CIs in Tooltips", FALSE),
              checkboxInput("show_trend_ci", "Show 95% CI on Trend Plot", FALSE),
              hr(), h5("Selected Area:"), verbatimTextOutput("click_info"),
              hr(), div(class = "state-avg-box", h5("Nebraska", style="margin:0 0 5px 0;"), uiOutput("state_avg_text"))))),

      # ===== PRIORITY ADDRESSES =====
      tabItem(tabName = "ranked",
        h3(class = "section-title", "Priority Addresses for Outreach"),
        p(class = "section-desc", "Addresses ranked by elevated count, BLL, or predicted risk. Colors indicate BLL severity."),
        fluidRow(
          box(width = 9,
              leafletOutput("ranked_map", height = 500), hr(),
              downloadButton("dl_ranked", "Download Displayed Addresses"), br(), br(),
              div(style = "height: 300px; overflow-y: auto;", DTOutput("ranked_table"))),
          box(width = 3, title = "Options",
              sliderInput("ranked_years", "Year Range:", min(years), max(years), c(min(years), max(years)), sep = ""),
              selectInput("ranked_lhd", "LHD:", c("Statewide"="All", lhd_names)),
              sliderInput("top_n", "Show Top N per LHD:", 5, 50, 10, 5),
              radioButtons("rank_by", "Rank By:",
                           c("Elevated Count"="n_elevated","Average BLL"="mean_bll","Predicted Risk"="pred_prob")),
              hr(), h5("Selected Address:"), verbatimTextOutput("ranked_click_info"),
              hr(), radioButtons("map_tiles", "Base Map:", c("Street"="street","Satellite"="sat"), inline=TRUE)))),

      # ===== CENSUS DATA =====
      tabItem(tabName = "census",
        h3(class = "section-title", "Census Neighborhood Characteristics"),
        p(class = "section-desc", "ACS data linked to lead testing outcomes."),
        fluidRow(
          box(width = 9,
              leafletOutput("census_map", height = 380), hr(),
              h5("Lead Outcomes by Quintile"),
              fluidRow(column(6, plotlyOutput("quintile_pct_plot", height=320)),
                       column(6, plotlyOutput("quintile_bll_plot", height=320))),
              hr(), h5("Quintile Summary"), downloadButton("dl_quintile", "Download Quintile Summary"),
              div(style="max-height:350px;overflow-y:auto;", DTOutput("census_quintile_table")),
              hr(), h5("Area Details"), downloadButton("dl_area_detail", "Download Area Details")),
          box(width = 3, title = "Options",
              selectInput("census_var", "Map Variable:",
                          c("Pre-1980 Housing (%)"="pct_pre1980","Pre-1950 Housing (%)"="pct_pre1950",
                            "Poverty Rate (%)"="pct_poverty","Renter-Occupied (%)"="pct_renter",
                            "Median Income ($)"="median_income","Median House Value ($)"="median_house_value",
                            "Non-White (%)"="pct_nonwhite","Housing Density (log)"="log_housing_density",
                            "Children Under 6 (%)"="pct_children_u6")),
              radioButtons("census_geo", "Level:", c("Tract"="tract","Block Group"="bg"), inline=TRUE),
              sliderInput("census_year", "ACS Vintage:", 2010, 2023, 2023, sep="", animate=animationOptions(interval=1200)),
              p(em("ACS 5-year estimate. Pre-2020 BG uses tract-level values."), style="font-size:10px;color:#888;"),
              checkboxInput("census_quintile", "Display by Quintile", FALSE),
              hr(), h5("Selected Area:"), verbatimTextOutput("census_click_info"),
              hr(), h5("Nebraska Overall:"), uiOutput("census_state_overall"),
              hr(), p(em("Data from CLPPP BLL records joined with ACS 5-year estimates (2013\u20132023)."),
                       style="font-size:10px;color:#888;")))),

      # ===== DATA TABLES =====
      tabItem(tabName = "tables",
        h3(class = "section-title", "Data Tables"),
        p(class = "section-desc", "Browse and download datasets."),
        fluidRow(
          column(3, sliderInput("table_years", "Years:", min(years), max(years), c(min(years), max(years)), sep="")),
          column(3, radioButtons("table_unit", "Unit:", c("Children","Addresses"), inline=TRUE)),
          column(3, radioButtons("table_geo", "Small Area:", c("Tract"="tract","Block Group"="bg"), inline=TRUE)),
          column(3, selectInput("table_rows", "Display Rows:", c(10,25,50,100), selected=25))),
        tabsetPanel(
          tabPanel("LHDs", br(),
                   fluidRow(column(6, downloadButton("dl_lhd_disp","Download Displayed")),
                            column(6, downloadButton("dl_lhd_all","Download All"))),
                   br(), DTOutput("lhd_table")),
          tabPanel("Counties", br(),
                   fluidRow(column(6, downloadButton("dl_county_disp","Download Displayed")),
                            column(6, downloadButton("dl_county_all","Download All"))),
                   br(), DTOutput("county_table")),
          tabPanel("Small Area", br(),
                   fluidRow(column(6, downloadButton("dl_geo_disp","Download Displayed")),
                            column(6, downloadButton("dl_geo_all","Download All"))),
                   br(), DTOutput("geo_table")),
          tabPanel("Addresses", br(),
                   fluidRow(column(6, downloadButton("dl_addr_disp","Download Displayed")),
                            column(6, downloadButton("dl_addr_all","Download All"))),
                   br(), DTOutput("addr_table")))),

      # ===== RISK MODEL =====
      tabItem(tabName = "model",
        h3(class = "section-title", "Predictive Risk Model"),
        p(class = "section-desc",
          "Mixed-effects logistic regression with geographic random effects. XGBoost and Logistic comparisons run alongside."),
        fluidRow(
          box(width = 4, title = "Model Settings",
              sliderInput("train_years", "Training Years:", min(years), max(years), c(2010, 2016), sep=""),
              checkboxInput("train_all", "Train on all years", FALSE),
              div(class="warning-box", em("Train on all years: predictions for future addresses but no held-out test set.")),
              hr(),
              radioButtons("model_addr_pop", "Training Population:",
                           c("2+ Children (subsequent child elevated)"="multi",
                             "All Addresses (test-period outcome)"="all"), selected="multi"),
              p(em("2+ Children: first child features predict whether subsequent child is elevated. ",
                   "All Addresses: first child features; outcome = any other child elevated in test period."),
                style="font-size:11px;color:#666;"),
              hr(),
              radioButtons("acs_mode", "ACS Census Linkage:",
                           c("Non-overlapping fixed (2013/2018/2023)"="nof",
                             "Year-matched (each address \u2192 closest ACS)"="ym",
                             "Endpoint split (train-end / test-end ACS)"="endpoint"),
                           selected="nof"),
              p(em("Fixed: uses only 2013, 2018, 2023 non-overlapping windows. ",
                   "Year-matched: closest ACS vintage (may have overlap). ",
                   "Endpoint: train addresses use train-end ACS, test use test-end ACS."),
                style="font-size:11px;color:#666;"),
              hr(),
              radioButtons("random_effect", "Random Intercept Level:",
                           c("Census Tract"="tract_geoid","Block Group"="bg_geoid"), selected="tract_geoid"),
              p(em("Tracts are larger (more stable); block groups are smaller (more local)."),
                style="font-size:11px;color:#666;"),
              hr(),
              checkboxInput("scale_vars", "Standardize variables", TRUE),
              div(class="warning-box", em("Standardization makes ORs comparable. Categorical/time variables never standardized.")),
              hr(),
              h5("Neighborhood Variables:"),
              checkboxGroupInput("neighborhood_covars", NULL,
                                 choices = c("Pre-1980 Housing (%)"="pct_pre1980","Pre-1950 Housing (%)"="pct_pre1950",
                                             "Poverty Rate (%)"="pct_poverty","Renter-Occupied (%)"="pct_renter",
                                             "Median Income ($)"="median_income","Median House Value ($)"="median_house_value",
                                             "Non-White (%)"="pct_nonwhite","Housing Density (log)"="log_housing_density"),
                                 selected = "pct_pre1980"),
              h5("Address Variables (First Child):"),
              checkboxGroupInput("addr_covars", NULL,
                                 choices = c("First Child's BLL"="first_bll","First Test Year"="first_year",
                                             "First Test Year (categorical)"="first_year_cat",
                                             "First Child Age (years, cat)"="first_age_yr_cat",
                                             "First Child Sex"="first_sex","Sample Type"="first_sample_type"),
                                 selected = c("first_bll","first_year")),
              conditionalPanel("input.addr_covars.indexOf('first_bll') > -1",
                checkboxInput("use_log_bll", "Use log(BLL) instead of raw BLL", FALSE)),
              hr(),
              h5("Models:"),
              p(em(strong("Multilevel"), " (GLMER), ", strong("Logistic"), " (GLM), ", strong("XGBoost"), " run simultaneously."),
                style="font-size:11px;color:#666;"),
              hr(),
              actionButton("fit_model", "Run Model", class="btn-primary btn-block")),
          box(width = 8, title = "Results",
              tabsetPanel(
                tabPanel("Summary", br(), uiOutput("interpretation")),
                tabPanel("Sample", br(),
                         checkboxInput("sample_bll_log", "Log scale for BLL plot", FALSE),
                         h5("Cross-Tabulation Summary"), DTOutput("sample_crosstab_table"), hr(),
                         fluidRow(column(6, plotlyOutput("sample_plot_addr", height=240)),
                                  column(6, plotlyOutput("sample_plot_tested", height=240))),
                         fluidRow(column(6, plotlyOutput("sample_plot_bll", height=240)),
                                  column(6, plotlyOutput("sample_plot_elev", height=240))),
                         br(), downloadButton("dl_sample_table", "Download Year-by-Year Summary")),
                tabPanel("Odds Ratios", br(), uiOutput("or_note"),
                         radioButtons("or_model_select", "Show Odds Ratios for:",
                                      c("Multilevel"="multilevel","Logistic (GLM)"="glm"), selected="multilevel", inline=TRUE),
                         plotlyOutput("coef_plot", height=350)),
                tabPanel("Correlations", br(),
                         p("High correlations (>0.70) can cause unstable estimates.", style="font-size:12px;"),
                         radioButtons("cor_method", "Correlation Method:", c("Pearson"="pearson","Spearman"="spearman"), inline=TRUE),
                         p(span("Yellow = borderline (0.50-0.70)", style="background:#fff3cd;padding:2px 5px;"),
                           span("Red = high (>0.70)", style="background:#f8d7da;padding:2px 5px;margin-left:10px;")),
                         DTOutput("cor_table"), br(), h5("VIF"), DTOutput("vif_table")),
                tabPanel("Model Comparison", br(),
                         uiOutput("model_comparison_ui"), hr(),
                         h5("Mean BLL by Screening Quintile"),
                         fluidRow(column(6, plotlyOutput("mean_bll_screened_plot", height=360)),
                                  column(6, plotlyOutput("capture_by_year_plot", height=360))),
                         hr(), h5("Cumulative Capture Curves"),
                         plotlyOutput("capture_overlay_plot", height=780), hr(),
                         fluidRow(column(6, h5("ROC Curve"), plotlyOutput("roc_youden_plot", height=400)),
                                  column(6, h5("AUC by Test Year"), plotlyOutput("auc_year_overlay_plot", height=400))),
                         hr(), h5("Calibration Over Time"),
                         checkboxGroupInput("cal_models", "Models:",
                                             c("Multilevel"="multi","Logistic"="glm","XGBoost"="xgb"),
                                             selected=c("multi","glm","xgb"), inline=TRUE),
                         plotlyOutput("cal_time_overlay_plot", height=400), hr(),
                         h5("Calibration: Risk Decile"), plotlyOutput("cal_decile_plot", height=400), hr(),
                         h5("Calibration: 10%-Width Bins"), plotlyOutput("cal_pctbin_plot", height=400), hr(),
                         h5("Confusion Matrix"), uiOutput("confusion_matrix_ui")),
                tabPanel("Technical", br(), verbatimTextOutput("model_summary"))))),
        fluidRow(
          box(width = 12, title = "Predicted Risk Map", status = "primary", solidHeader = TRUE,
              p("Average predicted risk by area. Darker = higher risk.", style="font-size:13px;color:#555;"),
              fluidRow(
                column(2, radioButtons("pred_type", "Data:", c("Test Set"="test","All Addresses"="full"))),
                column(2, selectInput("pred_model", "Model:", c("Multilevel"="glmer","Logistic (GLM)"="glm","XGBoost"="xgb"), selected="glmer")),
                column(2, selectInput("pred_geo", "Map Level:", c("County"="county","LHD"="lhd","Tract"="tract","Block Group"="bg"), selected="tract")),
                column(3, selectInput("pred_class", "Classification:", c("Continuous"="continuous","Quintile"="quintile","Decile"="decile"), selected="quintile")),
                column(2, selectInput("pred_palette", "Colors:", c("Red-Yellow"="risk","Viridis"="viridis","Cividis"="cividis"), selected="risk"))),
              leafletOutput("pred_map", height = 988))))
    )
  )
)

# ==============================================================================
# SERVER
# ==============================================================================

server <- function(input, output, session) {

  rv <- reactiveValues(
    fit = NULL, test_preds = NULL, full_preds = NULL,
    roc_obj = NULL, comparison = NULL, cor_matrix = NULL,
    clicked_id = NULL, census_id = NULL, ranked_clicked_id = NULL,
    model_run = FALSE, is_scaled = FALSE,
    train_sample = NULL, test_sample = NULL, ranked_display = NULL,
    xgb_fit = NULL, xgb_roc = NULL, xgb_comparison = NULL, xgb_test_preds = NULL,
    glm_fit = NULL, glm_roc = NULL, glm_test_preds = NULL)

  set_dl <- function(id, data_func, prefix) {
    output[[id]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".csv"),
      content = function(file) write.csv(data_func(), file, row.names = FALSE))
  }

  use_geo <- reactive(input$use_geo_mean == "geo")

  active_data <- reactive({
    if (!is.null(input$test_selection) && input$test_selection == "max" && !is.null(max_test)) max_test
    else first_test
  })

  get_address_stats <- reactive({
    if (input$addr_pop == "multi") address_stats else address_stats_all
  })

  observeEvent(input$addr_pop, {
    lbl <- if (input$addr_pop == "multi") "2+ Children" else "All Tested"
    showNotification(paste0("Address population: ", lbl), type = "message", duration = 4)
  }, ignoreInit = TRUE)

  get_all_covars <- reactive(c(input$neighborhood_covars %||% character(0), input$addr_covars))

  # Pre-aggregated geo stats: filter from prep output instead of recomputing
  get_stats <- reactive({
    yr <- input$year_filter; sup <- input$suppress_n
    geo_src <- if (isTRUE(input$test_selection == "max") && !is.null(geo_stats_max)) geo_stats_max
               else if (!is.null(geo_stats_first)) geo_stats_first
               else NULL
    # If no pre-aggregated data, fall back to runtime computation
    if (is.null(geo_src)) {
      dat <- active_data()
      ft <- dat |> filter(sample_year >= yr[1], sample_year <= yr[2])
      result <- list()
      for (nm in c("lhd","county")) result[[nm]] <- summarize_lead(ft, nm, sup, use_geo())
      active <- GEO_KEY[input$geo_level]
      if (active %in% c("tract","bg")) result[[active]] <- summarize_lead(ft, GEO_COL[input$geo_level], sup, use_geo())
      return(result)
    }
    # Use pre-aggregated: filter years, aggregate across years, apply suppression
    bll_col <- if (use_geo()) "mean_bll_geo" else "mean_bll"
    agg_geo <- function(level_key) {
      d <- geo_src |> filter(geo_level == level_key, sample_year >= yr[1], sample_year <= yr[2])
      if (nrow(d) == 0) return(tibble())
      d <- d |> filter(!is.na(.data[[bll_col]]), !is.na(n_children))
      if (nrow(d) == 0) return(tibble())
      d |> summarize(
        n_children = sum(n_children), n_elevated = sum(n_elevated),
        mean_bll = if (use_geo()) exp(weighted.mean(log(pmax(.data[[bll_col]], 0.1)), n_children))
                   else weighted.mean(.data[[bll_col]], n_children),
        n_addresses = sum(n_addresses, na.rm = TRUE),
        addr_elevated = sum(addr_elevated, na.rm = TRUE),
        county = first(county), .by = geo_id) |>
        mutate(pct_elevated = 100 * n_elevated / n_children, suppressed = n_children < sup)
    }
    result <- list(lhd = agg_geo("lhd"), county = agg_geo("county"))
    active <- GEO_KEY[input$geo_level]
    if (active == "tract") result$tract <- agg_geo("tract")
    else if (active == "bg") result$bg <- agg_geo("bg")
    result
  })

  get_filtered_data <- reactive({
    dat <- active_data()
    if (input$unit_type == "Addresses" && input$addr_pop == "multi") {
      multi <- dat |> summarize(n = n(), .by = AddressIdentifier) |> filter(n >= 2) |> pull(AddressIdentifier)
      dat <- dat |> filter(AddressIdentifier %in% multi)
    }
    dat
  })

  get_state_trend <- reactive(state_summary(get_filtered_data(), by_year=TRUE, by_address=(input$unit_type=="Addresses"), use_geo=use_geo()))
  get_state_all   <- reactive(state_summary(get_filtered_data(), by_address=(input$unit_type=="Addresses"), use_geo=use_geo()))

  # ===== SURVEILLANCE MAP =====

  output$surv_map <- renderLeaflet(base_map())

  output$state_avg_text <- renderUI({
    unit <- input$unit_type; var <- input$fill_var; yr <- input$year_filter
    dat <- active_data()
    ft_sel <- dat |> filter(sample_year >= yr[1], sample_year <= yr[2])
    avg_sel <- state_summary(ft_sel, by_address=(unit=="Addresses"), use_geo=use_geo())
    avg_all <- get_state_all()
    fmt <- function(a, v) {
      switch(v, mean_bll = sprintf("%.2f \u00b5g/dL", a$mean_bll),
        n_elevated = format(a$n_elevated, big.mark=","), n_children = format(a$n_children, big.mark=","),
        pct_elevated = sprintf("%.2f%%", a$pct_elevated),
        rate_per_1k = { acs_yrs <- acs_year_map(yr[1]):acs_year_map(yr[2])
          pop <- if (!is.null(acs_u6_ts$state)) sum(acs_u6_ts$state$n_children_u6[acs_u6_ts$state$acs_year %in% acs_yrs], na.rm=TRUE) else NA
          if (!is.na(pop) && pop > 0) sprintf("%.1f per 1k", a$n_elevated/pop*1000) else "N/A" },
        tested_rate_per_1k = { acs_yrs <- acs_year_map(yr[1]):acs_year_map(yr[2])
          pop <- if (!is.null(acs_u6_ts$state)) sum(acs_u6_ts$state$n_children_u6[acs_u6_ts$state$acs_year %in% acs_yrs], na.rm=TRUE) else NA
          if (!is.na(pop) && pop > 0) sprintf("%.1f per 1k", a$n_children/pop*1000) else "N/A" }, "N/A")
    }
    is_subset <- yr[1] != min(years) || yr[2] != max(years)
    if (is_subset) {
      ft_non <- dat |> filter(sample_year < yr[1] | sample_year > yr[2])
      avg_non <- state_summary(ft_non, by_address=(unit=="Addresses"), use_geo=use_geo())
      non_yrs <- sort(unique(ft_non$sample_year))
      non_lbl <- if (length(non_yrs) > 0) paste0(min(non_yrs),"\u2013",max(non_yrs)," (excl.)") else "Other"
      txt <- paste0("<b>", yr[1], "\u2013", yr[2], ":</b> ", fmt(avg_sel, var),
                    "<br><em style='font-size:12px;'>", non_lbl, ": ", fmt(avg_non, var), "</em>",
                    "<br><em style='font-size:12px;'>All Years: ", fmt(avg_all, var), "</em>")
    } else txt <- paste0("<b>All Years:</b> ", fmt(avg_all, var))
    ml <- if (use_geo()) "Geometric Mean" else "Arithmetic Mean"
    HTML(paste0(txt, "<br><em style='font-size:12px;'>BLL: ", ml, "</em>"))
  })

  # State trend plot
  output$state_avg_plot <- renderPlotly({
    var <- input$fill_var; yr_sel <- input$year_filter; unit <- input$unit_type
    df <- get_state_trend(); dat <- active_data()
    if (!var %in% names(df)) return(NULL)
    df$y_val <- df[[var]]
    if (!"sample_year" %in% names(df)) return(NULL)
    if (var %in% c("rate_per_1k","tested_rate_per_1k") && !is.null(acs_u6_ts$state)) {
      df <- df |> mutate(acs_yr = acs_year_map(sample_year)) |>
        left_join(acs_u6_ts$state |> rename(acs_yr = acs_year), by = "acs_yr")
      df$y_val <- if_else(!is.na(df$n_children_u6) & df$n_children_u6 > 0,
        if (var == "rate_per_1k") df$n_elevated/df$n_children_u6*1000 else df$n_children/df$n_children_u6*1000, NA_real_)
    }
    y_lab <- get_label(var)
    df <- df |> filter(sample_year >= yr_sel[1], sample_year <= yr_sel[2])
    if (nrow(df) == 0) return(NULL)
    show_ci <- isTRUE(input$show_trend_ci)
    if (show_ci) {
      df$ci_lo <- NA_real_; df$ci_hi <- NA_real_
      for (i in seq_len(nrow(df))) {
        ci <- compute_trend_ci(dat |> filter(sample_year == df$sample_year[i]), var, df$y_val[i])
        df$ci_lo[i] <- ci[1]; df$ci_hi[i] <- ci[2]
      }
    }
    fmt_val <- function(v) {
      if (var %in% c("n_elevated","n_children")) format(round(v), big.mark=",")
      else if (var %in% c("rate_per_1k","tested_rate_per_1k")) round(v, 1) else round(v, 2)
    }
    df$series <- "Nebraska"
    df$hover_text <- paste0("Nebraska\nYear: ", df$sample_year, "\n", y_lab, ": ", sapply(df$y_val, fmt_val))
    if (show_ci && any(!is.na(df$ci_lo)))
      df$hover_text <- paste0(df$hover_text, "\n95% CI: [", sapply(df$ci_lo, fmt_val), ", ", sapply(df$ci_hi, fmt_val), "]")
    # Area overlay
    area_df <- NULL
    if (!is.null(rv$clicked_id)) {
      g_col <- GEO_COL[input$geo_level]
      ft_area <- dat |> filter(.data[[g_col]] == rv$clicked_id, sample_year >= yr_sel[1], sample_year <= yr_sel[2])
      if (nrow(ft_area) > 0) {
        area_trend <- state_summary(ft_area, by_year=TRUE, by_address=(unit=="Addresses"), use_geo=use_geo())
        area_trend$y_val <- if (var %in% names(area_trend)) area_trend[[var]] else NA_real_
        area_df <- area_trend |> filter(!is.na(y_val))
        if (nrow(area_df) > 0) {
          area_df$series <- rv$clicked_id
          if (show_ci) {
            area_df$ci_lo <- NA_real_; area_df$ci_hi <- NA_real_
            for (i in seq_len(nrow(area_df))) {
              ci <- compute_trend_ci(ft_area |> filter(sample_year == area_df$sample_year[i]), var, area_df$y_val[i])
              area_df$ci_lo[i] <- ci[1]; area_df$ci_hi[i] <- ci[2]
            }
          }
          area_df$hover_text <- paste0(rv$clicked_id, "\nYear: ", area_df$sample_year, "\n", y_lab, ": ", sapply(area_df$y_val, fmt_val))
          if (show_ci && any(!is.na(area_df$ci_lo)))
            area_df$hover_text <- paste0(area_df$hover_text, "\n95% CI: [", sapply(area_df$ci_lo, fmt_val), ", ", sapply(area_df$ci_hi, fmt_val), "]")
        }
      }
    }
    plot_df <- df |> select(sample_year, y_val, series, hover_text, any_of(c("ci_lo","ci_hi")))
    if (!is.null(area_df) && nrow(area_df) > 0)
      plot_df <- bind_rows(plot_df, area_df |> select(sample_year, y_val, series, hover_text))
    plot_df <- plot_df |> arrange(series, sample_year) |> filter(!is.na(y_val))
    p <- plot_ly()
    if (show_ci && any(!is.na(df$ci_lo))) {
      ci_df <- df |> filter(!is.na(ci_lo)) |> arrange(sample_year)
      p <- p |> add_trace(x=ci_df$sample_year, y=ci_df$ci_hi, type="scatter", mode="lines",
                           line=list(width=0), showlegend=FALSE, hoverinfo="skip") |>
        add_trace(x=ci_df$sample_year, y=ci_df$ci_lo, type="scatter", mode="lines",
                  fill="tonexty", fillcolor="rgba(44,127,184,0.15)", line=list(width=0), showlegend=FALSE, hoverinfo="skip")
    }
    ne_df <- plot_df |> filter(series == "Nebraska")
    if (nrow(ne_df) > 0)
      p <- p |> add_trace(x=ne_df$sample_year, y=ne_df$y_val, type="scatter", mode="lines+markers",
                           name="Nebraska", line=list(color="#2c7fb8",width=2.5), marker=list(color="#2c7fb8",size=6),
                           text=ne_df$hover_text, hoverinfo="text")
    area_plot <- plot_df |> filter(series != "Nebraska")
    if (nrow(area_plot) > 0) {
      if (show_ci && !is.null(area_df) && any(!is.na(area_df$ci_lo))) {
        aci <- area_df |> filter(!is.na(ci_lo)) |> arrange(sample_year)
        p <- p |> add_trace(x=aci$sample_year, y=aci$ci_hi, type="scatter", mode="lines",
                             line=list(width=0), showlegend=FALSE, hoverinfo="skip") |>
          add_trace(x=aci$sample_year, y=aci$ci_lo, type="scatter", mode="lines",
                    fill="tonexty", fillcolor="rgba(215,48,39,0.12)", line=list(width=0), showlegend=FALSE, hoverinfo="skip")
      }
      p <- p |> add_trace(x=area_plot$sample_year, y=area_plot$y_val, type="scatter", mode="lines+markers",
                           name=area_plot$series[1], line=list(color="#d73027",width=2.5), marker=list(color="#d73027",size=6),
                           text=area_plot$hover_text, hoverinfo="text")
    }
    p |> layout(title=list(text=paste("Trend:", y_lab), font=list(size=14)),
                xaxis=list(title="Year"), yaxis=list(title=y_lab),
                legend=list(orientation="h", y=-0.15), margin=list(t=45, b=50))
  })

  # Trend summary table
  output$trend_summary_table <- renderDT({
    var <- input$fill_var; unit <- input$unit_type; sup <- input$suppress_n; dat <- active_data()
    g_col <- GEO_COL[input$geo_level]; trend_df <- get_state_trend()
    display_var <- if (unit=="Addresses") switch(var, n_children="n_addresses", n_elevated="addr_elevated", var) else var
    is_rate <- var %in% c("rate_per_1k","tested_rate_per_1k")
    # Pre-compute geo-level stats for all years
    geo_by_yr <- dat |> filter(!is.na(.data[[g_col]])) |>
      group_by(sample_year, geo_id = .data[[g_col]]) |>
      summarize(n_children=n(), n_elevated=sum(elevated, na.rm=TRUE),
                mean_bll=calc_mean_bll(result, use_geo()), .groups="drop") |>
      mutate(pct_elevated = 100*n_elevated/n_children, suppressed = n_children < sup)
    if (unit=="Addresses" && display_var %in% c("n_addresses","addr_elevated")) {
      addr_yr <- dat |> filter(!is.na(.data[[g_col]])) |>
        group_by(sample_year, geo_id=.data[[g_col]], AddressIdentifier) |>
        summarize(e=max(elevated, na.rm=TRUE), .groups="drop") |>
        group_by(sample_year, geo_id) |>
        summarize(n_addresses=n(), addr_elevated=sum(e>0, na.rm=TRUE), .groups="drop")
      geo_by_yr <- geo_by_yr |> left_join(addr_yr, by=c("sample_year","geo_id"))
    }
    if (is_rate) {
      geo_by_yr <- tryCatch({
        u6_src <- acs_u6_ts[[GEO_KEY[input$geo_level]]]
        if (!is.null(u6_src)) {
          geo_by_yr |> mutate(acs_yr=acs_year_map(sample_year)) |>
            left_join(u6_src |> rename(acs_yr=acs_year), by=c("geo_id","acs_yr")) |>
            mutate(rate_per_1k = if_else(!is.na(n_children_u6)&n_children_u6>0, n_elevated/n_children_u6*1000, NA_real_),
                   tested_rate_per_1k = if_else(!is.na(n_children_u6)&n_children_u6>0, n_children/n_children_u6*1000, NA_real_)) |>
            select(-acs_yr, -n_children_u6)
        } else geo_by_yr
      }, error = function(e) geo_by_yr)
    }
    geo_sf_list <- list(LHD=lhd_bounds, County=ne_counties, Tract=ne_tracts, `Block Group`=ne_bgs)
    yr_rows <- lapply(years, function(yr) {
      state_row <- trend_df |> filter(sample_year == yr)
      if (nrow(state_row) == 0) return(NULL)
      geo_yr <- geo_by_yr |> filter(sample_year == yr)
      vals <- geo_yr[[display_var]][!is.na(geo_yr[[display_var]]) & !geo_yr$suppressed]
      ft_yr <- dat |> filter(sample_year == yr); state_val <- state_row[[var]]
      if (is_rate && !is.null(acs_u6_ts$state)) {
        pop <- acs_u6_ts$state$n_children_u6[acs_u6_ts$state$acs_year == acs_year_map(yr)]
        if (length(pop) > 0 && pop > 0) state_val <- if (var=="rate_per_1k") state_row$n_elevated/pop*1000 else state_row$n_children/pop*1000
        else state_val <- NA_real_
      }
      ci <- compute_trend_ci(ft_yr, var, state_val)
      geo_sf <- geo_sf_list[[input$geo_level]]; total_geos <- if (!is.null(geo_sf)) nrow(geo_sf) else nrow(geo_yr)
      data.frame(Year=yr, `State Value`=round(state_val,2),
                 `95% CI Low`=if(!is.na(ci[1])) round(ci[1],2) else NA_real_,
                 `95% CI High`=if(!is.na(ci[2])) round(ci[2],2) else NA_real_,
                 `Std Dev`=if(length(vals)>1) round(sd(vals, na.rm=TRUE),2) else NA_real_,
                 Suppressed=sum(!is.na(geo_yr[[display_var]]) & geo_yr$suppressed, na.rm=TRUE),
                 `No Data`=max(0, total_geos - nrow(geo_yr)), check.names=FALSE)
    })
    datatable(do.call(rbind, yr_rows), options=list(dom="t", pageLength=20, scrollY="150px"), rownames=FALSE)
  })

  # Surveillance map observer
  observeEvent(c(input$geo_level, input$fill_var, input$unit_type, input$year_filter,
                 input$suppress_n, input$use_quintile, input$use_geo_mean, input$addr_pop,
                 input$test_selection, input$show_map_ci, input$map_opacity), {
    stats <- get_stats(); var <- input$fill_var; unit <- input$unit_type
    display_var <- if (unit=="Addresses") switch(var, n_children="n_addresses", n_elevated="addr_elevated", var) else var
    stat_key <- GEO_KEY[input$geo_level]
    stat_data <- stats[[stat_key]]
    if (is.null(stat_data) || nrow(stat_data) == 0) return()
    stat_data <- attach_rates(stat_data, input$geo_level, input$year_filter)
    GEO_SF <- list(LHD=lhd_bounds, County=ne_counties, Tract=ne_tracts, `Block Group`=ne_bgs)
    geo_bounds <- GEO_SF[[input$geo_level]]
    if (is.null(geo_bounds)) return()
    id_col <- GEO_ID[input$geo_level]
    geo <- geo_bounds |> left_join(stat_data, by = setNames("geo_id", id_col))
    geo_label <- GEO_LBL[input$geo_level]
    if (!display_var %in% names(geo)) return()
    vals <- geo[[display_var]][!is.na(geo[[display_var]]) & is.finite(geo[[display_var]]) &
                                 (is.null(geo$suppressed) | !geo$suppressed)]
    if (length(vals) == 0) vals <- c(0, 1)
    var_label <- get_label(var)
    if (input$use_quintile) {
      qp <- build_quintile_palette(vals)
      geo$fill_color <- assign_fill_colors(geo, display_var, function(v)
        qp$colors[max(1, min(findInterval(v, qp$breaks, rightmost.closed=TRUE), 5))])
      legend_labels <- c(qp$labels, "Suppressed", "No Data"); legend_colors <- c(qp$colors, "#666666", "#cccccc")
    } else {
      domain <- switch(var, mean_bll=c(0,5), pct_elevated=c(0,30),
        rate_per_1k=c(0, max(30, max(vals, na.rm=TRUE))),
        tested_rate_per_1k=c(0, max(500, max(vals, na.rm=TRUE))),
        c(0, max(vals, na.rm=TRUE)))
      pal <- colorNumeric("YlOrRd", domain=domain, na.color="#cccccc")
      geo$fill_color <- assign_fill_colors(geo, display_var, function(v) pal(min(v, domain[2])))
      breaks <- seq(domain[1], domain[2], length.out=6)
      fmt_brk <- function(x) {
        if (var %in% c("n_elevated","n_children","n_addresses","addr_elevated")) format(round(x), big.mark=",")
        else if (var == "pct_elevated") sprintf("%.1f%%", x)
        else if (var %in% c("rate_per_1k","tested_rate_per_1k")) sprintf("%.1f", x)
        else format(round(x,1), big.mark=",")
      }
      legend_labels <- c(sapply(1:5, \(i) paste0(fmt_brk(breaks[i]),"-",fmt_brk(breaks[i+1]))), "Suppressed", "No Data")
      legend_colors <- c(pal(breaks[1:5] + diff(breaks[1:2])/2), "#666666", "#cccccc")
    }
    sup_flag <- ifelse(!is.null(geo$suppressed) & !is.na(geo$suppressed) & geo$suppressed, " [Suppressed]", "")
    ml <- bll_mean_label(use_geo())
    na_bll <- ifelse(!is.na(geo$mean_bll), paste0(sprintf("%.2f", geo$mean_bll), " \u00b5g/dL ", ml), "N/A")
    ci_html_pct <- ""; ci_html_bll <- ""
    if (isTRUE(input$show_map_ci)) {
      n_t <- geo$n_children; n_e <- geo$n_elevated; bll <- geo$mean_bll; valid <- !is.na(n_t) & n_t > 0
      ci_pct_lo <- ci_pct_hi <- ci_bll_lo <- ci_bll_hi <- rep(NA_real_, nrow(geo))
      if (any(valid)) {
        a <- 0.05
        ci_pct_lo[valid] <- ifelse(n_e[valid]==0, 0, qbeta(a/2, n_e[valid], n_t[valid]-n_e[valid]+1))*100
        ci_pct_hi[valid] <- ifelse(n_e[valid]==n_t[valid], 100, qbeta(1-a/2, n_e[valid]+1, n_t[valid]-n_e[valid])*100)
        bll_valid <- valid & !is.na(bll) & n_t >= 2
        if (any(bll_valid)) {
          state_cv <- sd(active_data()$result, na.rm=TRUE) / mean(active_data()$result, na.rm=TRUE)
          state_cv <- min(max(state_cv, 0.3), 1.5)
          se <- bll[bll_valid] * state_cv / sqrt(n_t[bll_valid])
          ci_bll_lo[bll_valid] <- pmax(0, bll[bll_valid] - 1.96*se)
          ci_bll_hi[bll_valid] <- bll[bll_valid] + 1.96*se
        }
      }
      geo$ci_pct <- ifelse(!is.na(ci_pct_lo), paste0("<em style='font-size:11px;color:#555;'>95% CI: [",
                           sprintf("%.1f",ci_pct_lo),"%, ",sprintf("%.1f",ci_pct_hi),"%]</em>"), "")
      geo$ci_bll <- ifelse(!is.na(ci_bll_lo), paste0("<em style='font-size:11px;color:#555;'>95% CI: [",
                           sprintf("%.2f",ci_bll_lo),", ",sprintf("%.2f",ci_bll_hi),"]</em>"), "")
    } else { geo$ci_pct <- ""; geo$ci_bll <- "" }
    rate_e <- ifelse("rate_per_1k" %in% names(geo) & !is.na(geo$rate_per_1k), sprintf("%.1f",geo$rate_per_1k), "N/A")
    rate_t <- ifelse("tested_rate_per_1k" %in% names(geo) & !is.na(geo$tested_rate_per_1k), sprintf("%.1f",geo$tested_rate_per_1k), "N/A")
    acs_pop <- ifelse("n_acs_u6_avg" %in% names(geo) & !is.na(geo$n_acs_u6_avg), format(round(geo$n_acs_u6_avg), big.mark=","), "N/A")
    geo$tip <- paste0("<b>", geo_label, ": ", geo[[id_col]], "</b>", sup_flag, "<br>",
      "Elevated: ", format(geo$n_elevated, big.mark=","), "<br>Tested: ", format(geo$n_children, big.mark=","), "<br>",
      "% Elevated: ", ifelse(!is.na(geo$pct_elevated), sprintf("%.1f%%",geo$pct_elevated), "N/A"),
      ifelse(geo$ci_pct != "", paste0("<br>",geo$ci_pct), ""), "<br>Average BLL: ", na_bll,
      ifelse(geo$ci_bll != "", paste0("<br>",geo$ci_bll), ""), "<br>",
      "Elev/1k: ", rate_e, "<br>Tested/1k: ", rate_t, "<br>ACS U6 Pop: ", acs_pop,
      ifelse(rep(unit=="Addresses" & !is.null(geo$n_addresses), nrow(geo)),
        paste0("<br><em>Addr: ", format(geo$n_addresses,big.mark=","), " | Addr Elev: ", format(geo$addr_elevated,big.mark=","), "</em>"), ""))
    leafletProxy("surv_map") |> clearShapes() |> clearControls() |>
      addPolygons(data=geo, fillColor=~fill_color, fillOpacity=input$map_opacity/100,
                  color="white", weight=1, layerId=as.character(geo[[id_col]]), label=lapply(geo$tip, HTML)) |>
      addLegend("bottomright", colors=legend_colors, labels=legend_labels, title=var_label, opacity=0.8)
  }, ignoreNULL = FALSE)

  observeEvent(input$surv_map_shape_click, { rv$clicked_id <- input$surv_map_shape_click$id })

  output$click_info <- renderPrint({
    if (is.null(rv$clicked_id)) return(cat("Click a region for details"))
    yr <- input$year_filter; stats_list <- get_stats()
    d <- stats_list[[GEO_KEY[input$geo_level]]] |> filter(geo_id == rv$clicked_id)
    if (nrow(d) == 0) return(cat("No data"))
    d <- attach_rates(d, input$geo_level, yr)
    ft_sel <- active_data() |> filter(sample_year >= yr[1], sample_year <= yr[2])
    state <- state_summary(ft_sel, by_address=(input$unit_type=="Addresses"), use_geo=use_geo())
    ml <- bll_mean_label(use_geo()); fn <- function(x) format(x, big.mark=",")
    fp <- function(x) if (!is.na(x)) sprintf("%.2f%%",x) else "N/A"
    fb <- function(x) if (!is.na(x)) sprintf("%.2f",x) else "N/A"
    cat(GEO_LBL[input$geo_level], ":", rv$clicked_id, "\n")
    if (d$suppressed) cat("[Suppressed]\n")
    cat(sprintf("\n%-22s %-18s\n", "--- Selected ---", "--- Nebraska ---"))
    if (input$unit_type == "Children") {
      cat(sprintf("Elevated:   %-12s Elevated:   %s\n", fn(d$n_elevated), fn(state$n_elevated)))
      cat(sprintf("Tested:     %-12s Tested:     %s\n", fn(d$n_children), fn(state$n_children)))
    } else {
      cat(sprintf("Addr Elev:  %-12s Elevated:   %s\n", fn(d$addr_elevated), fn(state$n_elevated)))
      cat(sprintf("Addr Tested:%-12s Tested:     %s\n", fn(d$n_addresses), fn(state$n_children)))
    }
    cat(sprintf("%% Elevated: %-12s %% Elevated: %s\n",
        fp(if(input$unit_type=="Addresses") 100*d$addr_elevated/d$n_addresses else d$pct_elevated), fp(state$pct_elevated)))
    cat(sprintf("Avg BLL:    %-12s Avg BLL:    %s\n", fb(d$mean_bll), fb(state$mean_bll)))
    cat(sprintf("BLL Average: %s | Years: %d-%d\n", ml, yr[1], yr[2]))
    if ("rate_per_1k" %in% names(d) && !is.na(d$rate_per_1k))
      cat(sprintf("Elev/1k: %.1f | Tested/1k: %.1f | ACS U6: %s\n", d$rate_per_1k, d$tested_rate_per_1k, fn(d$n_acs_u6)))
  })

  # ===== PRIORITY ADDRESSES =====

  output$ranked_map <- renderLeaflet({
    tile <- if (input$map_tiles == "sat") providers$Esri.WorldImagery else providers$CartoDB.Positron
    leaflet() |> addProviderTiles(tile) |> setView(-99.5, 41.5, 7) |>
      addPolygons(data=lhd_bounds, fillColor="transparent", color="#2c7fb8", weight=2, label=~lhd) |>
      addLegend("bottomright", colors=c("#67000d","#a50f15","#d73027","#fc8d59","#c9a84c"),
                labels=c("\u226520 \u00b5g/dL","10\u201320","5\u201310","3.5\u20135","<3.5"), title="BLL Color Key", opacity=0.9)
  })
  observeEvent(input$map_tiles, {
    tile <- if (input$map_tiles == "sat") providers$Esri.WorldImagery else providers$CartoDB.Positron
    leafletProxy("ranked_map") |> clearTiles() |> addProviderTiles(tile)
  })

  get_ranked_addresses <- reactive({
    addr <- get_address_stats() |> filter(!is.na(lhd))
    yr <- input$ranked_years; addr <- addr |> filter(last_year >= yr[1], first_year <= yr[2])
    if (use_geo()) {
      geo_means <- active_data() |> filter(sample_year >= yr[1], sample_year <= yr[2]) |>
        group_by(AddressIdentifier) |> summarize(mean_bll = calc_mean_bll(result, TRUE), .groups = "drop")
      addr <- addr |> select(-mean_bll) |> left_join(geo_means, by = "AddressIdentifier")
    }
    if (!is.null(rv$full_preds)) {
      addr <- addr |> left_join(rv$full_preds |> select(AddressIdentifier, pred_prob=pred) |>
                                  distinct(AddressIdentifier, .keep_all=TRUE), by="AddressIdentifier")
    } else addr$pred_prob <- NA_real_
    if (input$rank_by == "pred_prob" && all(is.na(addr$pred_prob))) return(addr[0, ])
    if (input$ranked_lhd != "All") addr <- addr |> filter(lhd == input$ranked_lhd)
    addr <- addr |> mutate(risk_decile = ifelse(!is.na(pred_prob), ntile(pred_prob, 10), NA_integer_))
    addr |> group_by(lhd) |> arrange(desc(.data[[input$rank_by]]), .by_group=TRUE) |>
      mutate(rank = row_number()) |> ungroup() |> filter(rank <= input$top_n) |>
      mutate(color = get_bll_color(mean_bll))
  })

  observeEvent(c(input$ranked_lhd, input$top_n, input$rank_by, input$ranked_years,
                 input$use_geo_mean, input$addr_pop, input$map_opacity), {
    addr <- get_ranked_addresses()
    if (is.null(addr) || nrow(addr) == 0) {
      if (input$rank_by == "pred_prob") showNotification("Run the Risk Model first.", type="warning"); return()
    }
    rv$ranked_display <- addr; ml <- bll_mean_label(use_geo())
    addr$tip <- paste0("<b>#",addr$rank," - ",addr$AddressIdentifier,"</b><br>",
      "LHD: ",addr$lhd,"<br>County: ",addr$county,"<br>",
      "Tested: ",format(addr$n_children,big.mark=","),"<br>Elevated: ",format(addr$n_elevated,big.mark=","),"<br>",
      "Avg BLL ",ml,": ",round(addr$mean_bll,2)," \u00b5g/dL<br>Max BLL: ",round(addr$max_bll,2)," \u00b5g/dL<br>",
      "Years: ",addr$first_year," \u2013 ",addr$last_year,
      ifelse(!is.na(addr$pred_prob), paste0("<br>Risk: ",round(addr$pred_prob*100,1),"%"), ""))
    leafletProxy("ranked_map") |> clearMarkers() |>
      addCircleMarkers(data=addr, lng=~lng, lat=~lat, radius=8,
                       fillColor=~color, color="#333", weight=1, fillOpacity=input$map_opacity/100,
                       layerId=~AddressIdentifier, label=lapply(addr$tip, HTML), options=pathOptions(pane="markerPane")) |>
      addLabelOnlyMarkers(data=addr, lng=~lng, lat=~lat, label=as.character(addr$rank),
                          labelOptions=labelOptions(noHide=TRUE, direction="center", textOnly=TRUE, textsize="10px",
                                                    style=list("font-weight"="bold","color"="white","pointer-events"="none")))
  }, ignoreNULL = FALSE)

  observeEvent(input$ranked_map_marker_click, { rv$ranked_clicked_id <- input$ranked_map_marker_click$id })

  output$ranked_click_info <- renderPrint({
    if (is.null(rv$ranked_clicked_id)) return(cat("Click a marker for details"))
    d <- get_ranked_addresses() |> filter(AddressIdentifier == rv$ranked_clicked_id)
    if (nrow(d) == 0) return(cat("No data")); ml <- bll_mean_label(use_geo())
    cat("Address:", d$AddressIdentifier[1], "\nLHD:", d$lhd[1], "| County:", d$county[1], "\n")
    cat("Tested:", d$n_children[1], "| Elevated:", d$n_elevated[1], "\n")
    cat("Avg BLL", ml, ":", round(d$mean_bll[1],2), "| Max:", round(d$max_bll[1],2), "\n")
    cat("Years:", d$first_year[1], "-", d$last_year[1], "\n")
    if (!is.na(d$pred_prob[1])) cat("Risk:", sprintf("%.1f%%",d$pred_prob[1]*100), "| Decile:", d$risk_decile[1], "\n")
  })

  show_address_tests <- function(addr_id) {
    src <- if (!is.null(all_tests)) all_tests else first_test
    tests <- src |> filter(AddressIdentifier == addr_id)
    if (nrow(tests) == 0) { showNotification("No records", type="warning"); return() }
    for (col in c("SEX","RACE_ETH","SAMPLE_TYPE"))
      if (!col %in% names(tests)) tests[[col]] <- NA_character_
    tests <- tests |> mutate(
      SEX = if_else(is.na(SEX)|SEX=="", "Not recorded", as.character(SEX)),
      RACE_ETH = if_else(is.na(RACE_ETH)|RACE_ETH=="", "Not recorded", as.character(RACE_ETH)),
      SAMPLE_TYPE = if_else(is.na(SAMPLE_TYPE)|SAMPLE_TYPE=="", "Not recorded", as.character(SAMPLE_TYPE)))
    cols <- intersect(c("PATIENT_LOCAL_ID","SAMPLE_DATE","result","elevated","AGE_YR","AGE_MO","SEX","RACE_ETH","SAMPLE_TYPE","sample_year"), names(tests))
    display <- tests |> select(all_of(cols)) |> arrange(SAMPLE_DATE)
    showModal(modalDialog(
      title=paste("Test Records:", tests$AddressIdentifier[1]), size="l", easyClose=TRUE,
      p(em(paste(nrow(display), "records at this address"))),
      DTOutput("modal_test_table"), footer=tagList(downloadButton("dl_modal_tests","Download"), modalButton("Close"))))
    output$modal_test_table <- renderDT(make_dt(display, 25))
    output$dl_modal_tests <- downloadHandler(
      filename=function() paste0("tests_",addr_id,"_",Sys.Date(),".csv"),
      content=function(file) write.csv(tests, file, row.names=FALSE))
  }

  output$ranked_table <- renderDT({
    if (is.null(rv$ranked_display)) return(NULL); d <- rv$ranked_display; ml <- bll_mean_label(use_geo())
    if (!is.null(rv$full_preds) && "single_child_flag" %in% names(rv$full_preds)) {
      d <- d |> left_join(rv$full_preds |> select(AddressIdentifier, single_child_flag) |> distinct(), by="AddressIdentifier")
    } else d$single_child_flag <- FALSE
    d$single_child_flag[is.na(d$single_child_flag)] <- FALSE
    tbl <- d |> transmute(Rank=rank, Address=AddressIdentifier, LHD=lhd, County=county,
      Children=n_children, Elevated=n_elevated, !!paste("Avg BLL",ml):=round(mean_bll,2), `Max BLL`=round(max_bll,2),
      `Risk Decile`=if_else(!is.na(risk_decile), as.character(risk_decile), "\u2014"),
      `Pred Risk`=if_else(!is.na(pred_prob), paste0(sprintf("%.1f%%",pred_prob*100), if_else(single_child_flag,"*","")), "\u2014"),
      `View`=paste0('<button class="btn btn-xs btn-info act-view-tests" data-addr="',htmltools::htmlEscape(AddressIdentifier,TRUE),'">View</button>'),
      `Risk`=paste0('<button class="btn btn-xs btn-warning act-risk-breakdown" data-addr="',htmltools::htmlEscape(AddressIdentifier,TRUE),'"',
                    if_else(!is.na(pred_prob),'>',' disabled>'), if_else(!is.na(pred_prob),'Explain','\u2014'),'</button>'))
    datatable(tbl, escape=FALSE, selection="single", options=list(pageLength=10, scrollX=TRUE, dom="tp"), rownames=FALSE,
              caption=if(any(d$single_child_flag,na.rm=TRUE))
                htmltools::tags$caption(style="font-size:11px;color:#888;font-style:italic;","* Single-child address") else NULL)
  })
  observeEvent(input$view_tests_click, show_address_tests(input$view_tests_click))
  set_dl("dl_ranked", function() { if(is.null(rv$ranked_display)) data.frame() else rv$ranked_display |> select(-any_of(c("tip","color"))) }, "priority_addresses")

  # ===== CENSUS DATA =====

  get_census_geo_data <- reactive({
    geo <- input$census_geo; var <- input$census_var; sup <- input$suppress_n
    sel_year <- input$census_year %||% 2023
    census_df <- if (!is.null(tract_census_ts) && !is.null(bg_census_ts)) {
      if (geo=="tract") tract_census_ts |> filter(acs_year==sel_year) |> select(-acs_year)
      else bg_census_ts |> filter(acs_year==sel_year) |> select(-acs_year)
    } else { if (geo=="tract") tract_census else bg_census }
    geo_bounds <- if (geo=="tract") ne_tracts else ne_bgs
    geo_col <- if (geo=="tract") "tract_geoid" else "bg_geoid"
    dat <- active_data()
    lead_df <- dat |> filter(!is.na(.data[[geo_col]]), nchar(.data[[geo_col]]) >= 11) |>
      group_by(geo_id=.data[[geo_col]]) |>
      summarize(n_children=n(), n_elevated=sum(elevated,na.rm=TRUE),
                mean_bll=calc_mean_bll(result,use_geo()), .groups="drop") |>
      mutate(pct_elevated=100*n_elevated/n_children, suppressed=n_children<sup)
    var_missing <- !is.null(census_df) && (!var %in% names(census_df) || all(is.na(census_df[[var]])))
    if (var_missing && !is.null(tract_census_ts) && geo=="bg") {
      tract_yr <- tract_census_ts |> filter(acs_year==sel_year) |> select(-acs_year)
      if (var %in% names(tract_yr)) {
        census_df <- census_df |> mutate(.tid=substr(bg_geoid,1,11)) |>
          left_join(tract_yr |> select(tract_geoid, all_of(var)), by=c(".tid"="tract_geoid"), suffix=c("",".t")) |>
          select(-.tid)
        if (paste0(var,".t") %in% names(census_df)) {
          census_df[[var]] <- coalesce(census_df[[var]], census_df[[paste0(var,".t")]]); census_df[[paste0(var,".t")]] <- NULL
        }
      }
    }
    combined <- lead_df |> left_join(census_df, by=setNames(geo_col, "geo_id"))
    list(geo_sf=geo_bounds |> left_join(combined, by=setNames("geo_id","GEOID")), combined=combined, id_col="GEOID")
  })

  output$census_map <- renderLeaflet(base_map())

  observeEvent(c(input$census_var, input$census_geo, input$suppress_n, input$census_quintile, input$use_geo_mean, input$census_year, input$map_opacity), {
    cd <- get_census_geo_data(); geo <- cd$geo_sf; var <- input$census_var
    if (is.null(geo) || !var %in% names(geo)) return()
    vals <- geo[[var]][!is.na(geo[[var]]) & is.finite(geo[[var]])]
    if (length(vals) == 0) return()
    var_label <- get_label(var)
    if (input$census_quintile) {
      qp <- build_quintile_palette(vals, c("#eff3ff","#bdd7e7","#6baed6","#3182bd","#08519c"))
      geo$fill_color <- assign_fill_colors(geo, var, function(v)
        qp$colors[max(1, min(findInterval(v, qp$breaks, rightmost.closed=TRUE), 5))])
      leg_l <- c(qp$labels, "No Data"); leg_c <- c(qp$colors, "#cccccc")
    } else {
      pal <- colorNumeric("Blues", domain=range(vals, na.rm=TRUE), na.color="#cccccc")
      geo$fill_color <- ifelse(is.na(geo[[var]]), "#cccccc", sapply(geo[[var]], pal))
      leg_l <- NULL; leg_c <- NULL
    }
    ml <- bll_mean_label(use_geo())
    geo$tip <- paste0("<div style='font-size:13px;line-height:1.6;'>",
      "<b>",geo$GEOID,"</b><br><b>",var_label,":</b> ",ifelse(!is.na(geo[[var]]),round(geo[[var]],2),"N/A"),"<br>",
      "<b>Mean BLL</b> ",ml,": ",ifelse(!is.na(geo$mean_bll),paste0(round(geo$mean_bll,2)," \u00b5g/dL"),"N/A"),"<br>",
      "<b>% Elevated:</b> ",ifelse(!is.na(geo$pct_elevated),paste0(round(geo$pct_elevated,1),"%"),"N/A"),"<br>",
      "<b>Tested:</b> ",ifelse(!is.na(geo$n_children),format(geo$n_children,big.mark=","),"N/A"),
      " | <b>Elevated:</b> ",ifelse(!is.na(geo$n_elevated),format(geo$n_elevated,big.mark=","),"N/A"),"</div>")
    proxy <- leafletProxy("census_map") |> clearShapes() |> clearControls() |>
      addPolygons(data=geo, fillColor=~fill_color, fillOpacity=input$map_opacity/100,
                  color="white", weight=0.5, layerId=~GEOID, label=lapply(geo$tip, HTML))
    if (input$census_quintile) proxy |> addLegend("bottomright", title=var_label, opacity=0.8, colors=leg_c, labels=leg_l)
    else proxy |> addLegend("bottomright", title=var_label, opacity=0.8,
                            pal=colorNumeric("Blues",range(vals,na.rm=TRUE)), values=seq(min(vals),max(vals),length.out=5))
  }, ignoreNULL = FALSE)

  observeEvent(input$census_map_shape_click, { rv$census_id <- input$census_map_shape_click$id })

  output$census_click_info <- renderPrint({
    if (is.null(rv$census_id)) return(cat("Click an area for details"))
    cd <- get_census_geo_data(); row <- cd$combined |> filter(geo_id == rv$census_id)
    if (nrow(row) == 0) return(cat("No data")); ml <- bll_mean_label(use_geo()); var <- input$census_var
    dat <- active_data(); state <- state_summary(dat, use_geo=use_geo())
    fn <- function(x) format(x, big.mark=","); fp <- function(x) if(!is.na(x)) sprintf("%.2f%%",x) else "N/A"
    fb <- function(x) if(!is.na(x)) sprintf("%.2f",x) else "N/A"
    cat("Area:", rv$census_id, "\n", get_label(var), ":", if(var %in% names(row) && !is.na(row[[var]])) round(row[[var]],2) else "N/A", "\n")
    cat(sprintf("\n%-22s %-18s\n","--- Selected ---","--- Nebraska ---"))
    cat(sprintf("Elevated:   %-12s Elevated:   %s\n", fn(row$n_elevated), fn(state$n_elevated)))
    cat(sprintf("Tested:     %-12s Tested:     %s\n", fn(row$n_children), fn(state$n_children)))
    cat(sprintf("%% Elevated: %-12s %% Elevated: %s\n", fp(row$pct_elevated), fp(state$pct_elevated)))
    cat(sprintf("Avg BLL:    %-12s Avg BLL:    %s\nBLL Average: %s\n", fb(row$mean_bll), fb(state$mean_bll), ml))
  })

  output$census_state_overall <- renderUI({
    dat <- active_data(); s <- state_summary(dat, use_geo=use_geo()); ml <- bll_mean_label(use_geo())
    sel_year <- input$census_year %||% 2023
    census_df <- if (!is.null(tract_census_ts)) {
      if (input$census_geo=="tract") tract_census_ts |> filter(acs_year==sel_year) |> select(-acs_year)
      else bg_census_ts |> filter(acs_year==sel_year) |> select(-acs_year)
    } else { if (input$census_geo=="tract") tract_census else bg_census }
    sv <- if (!is.null(census_df) && input$census_var %in% names(census_df)) round(mean(census_df[[input$census_var]], na.rm=TRUE),1) else "N/A"
    div(class="state-avg-box",
        p(strong("Tested:"), format(s$n_children,big.mark=","), " | ", strong("Elevated:"), format(s$n_elevated,big.mark=",")),
        p(strong("% Elevated:"), sprintf("%.2f%%",s$pct_elevated), " | ", strong("Avg BLL"),ml,":", sprintf("%.2f",s$mean_bll)),
        p(strong(get_label(input$census_var), paste0("(state avg, ACS ",sel_year,"):")), sv))
  })

  get_quintile_data <- reactive({
    cd <- get_census_geo_data(); var <- input$census_var; d <- cd$combined
    if (!var %in% names(d)) return(NULL)
    d <- d |> filter(!is.na(.data[[var]]), !is.na(pct_elevated)); if (nrow(d) < 10) return(NULL)
    d$quintile <- ntile(d[[var]], 5)
    sel_year <- input$census_year %||% 2023
    geo_level <- if (input$census_geo=="tract") "tract" else "bg"
    u6_src <- if (!is.null(acs_u6_ts)) acs_u6_ts[[geo_level]] else NULL
    if (!is.null(u6_src) && nrow(u6_src) > 0) {
      matched <- intersect(acs_year_map(sel_year), unique(u6_src$acs_year))
      if (length(matched)==0) matched <- unique(u6_src$acs_year)[which.min(abs(unique(u6_src$acs_year)-sel_year))]
      u6_agg <- u6_src |> filter(acs_year %in% matched) |> summarize(n_children_u6=sum(n_children_u6,na.rm=TRUE), .by=geo_id)
      if (nrow(u6_agg) > 0) d <- d |> left_join(u6_agg, by="geo_id")
    }
    if ((!"n_children_u6" %in% names(d) || all(is.na(d$n_children_u6))) && geo_level=="bg") {
      u6_tr <- if (!is.null(acs_u6_ts)) acs_u6_ts[["tract"]] else NULL
      if (!is.null(u6_tr)) {
        u6_yr <- u6_tr |> filter(acs_year==sel_year) |> select(geo_id, n_children_u6)
        if (nrow(u6_yr) > 0) d <- d |> select(-any_of("n_children_u6")) |>
          mutate(.tid=substr(geo_id,1,11)) |> left_join(u6_yr, by=c(".tid"="geo_id")) |> select(-.tid)
      }
    }
    if (!"n_children_u6" %in% names(d)) d$n_children_u6 <- NA_real_
    d |> group_by(quintile) |> summarize(
      pct_elev=mean(pct_elevated,na.rm=TRUE), avg_bll=calc_mean_bll(mean_bll,use_geo()),
      n_areas=n(), tested=sum(n_children,na.rm=TRUE), elevated=sum(n_elevated,na.rm=TRUE),
      acs_u6=sum(n_children_u6,na.rm=TRUE),
      elev_per_1k=if_else(sum(n_children_u6,na.rm=TRUE)>0, round(sum(n_elevated,na.rm=TRUE)/sum(n_children_u6,na.rm=TRUE)*1000,1), NA_real_),
      range_lo=round(min(.data[[var]],na.rm=TRUE),1), range_hi=round(max(.data[[var]],na.rm=TRUE),1), .groups="drop")
  })

  output$quintile_pct_plot <- renderPlotly({
    qd <- get_quintile_data(); if (is.null(qd)) return(NULL); vl <- get_label(input$census_var)
    plot_ly(qd, x=~paste0("Q",quintile), y=~pct_elev, type="bar", marker=list(color="#2c7fb8"),
            text=~paste0(round(pct_elev,1),"%"), textposition="outside", textfont=list(size=11),
            hovertext=~paste0("Q",quintile," (",range_lo,"-",range_hi,")\n% Elevated: ",round(pct_elev,2),"\nAreas: ",n_areas), hoverinfo="text") |>
      layout(title=list(text=paste("% Elevated by",vl),font=list(size=12)),
             xaxis=list(title=paste(vl,"Quintile")), yaxis=list(visible=FALSE,range=c(0,max(qd$pct_elev,na.rm=TRUE)*1.25)), margin=list(l=20,b=70,t=40))
  })

  output$quintile_bll_plot <- renderPlotly({
    qd <- get_quintile_data(); if (is.null(qd)) return(NULL); vl <- get_label(input$census_var)
    has_acs <- "elev_per_1k" %in% names(qd) && !all(is.na(qd$elev_per_1k))
    if (has_acs) {
      plot_ly(qd, x=~paste0("Q",quintile), y=~elev_per_1k, type="bar", marker=list(color="#d73027"),
              text=~paste0(round(elev_per_1k,1)," per 1k"), textposition="outside", textfont=list(size=10),
              hovertext=~paste0("Q",quintile,"\nElev/1k: ",round(elev_per_1k,1),"\nACS U6: ",format(acs_u6,big.mark=",")), hoverinfo="text") |>
        layout(title=list(text=paste("Elevated per 1k U6 by",vl),font=list(size=12)),
               xaxis=list(title=paste(vl,"Quintile")), yaxis=list(visible=FALSE,range=c(0,max(qd$elev_per_1k,na.rm=TRUE)*1.25)), margin=list(l=20,b=70,t=40))
    } else {
      ml <- bll_mean_label(use_geo())
      plot_ly(qd, x=~paste0("Q",quintile), y=~avg_bll, type="bar", marker=list(color="#d73027"),
              text=~paste0(round(avg_bll,2)," \u00b5g/dL"), textposition="outside", textfont=list(size=10),
              hovertext=~paste0("Q",quintile,"\nAvg BLL ",ml,": ",round(avg_bll,2)), hoverinfo="text") |>
        layout(title=list(text=paste("Avg BLL by",vl),font=list(size=12)),
               xaxis=list(title=paste(vl,"Quintile")), yaxis=list(visible=FALSE,range=c(0,max(qd$avg_bll,na.rm=TRUE)*1.25)), margin=list(l=20,b=70,t=40))
    }
  })

  output$census_quintile_table <- renderDT({
    qd <- get_quintile_data(); if (is.null(qd)) return(NULL)
    datatable(qd |> transmute(Quintile=paste0("Q",quintile), Range=paste0(range_lo,"-",range_hi),
                               Areas=n_areas, Tested=tested, Elevated=elevated, `% Elevated`=round(pct_elev,2),
                               `ACS Under-6`=acs_u6, `Elev/1k U6`=elev_per_1k, `Avg BLL`=round(avg_bll,2)),
              options=list(dom="t", pageLength=5), rownames=FALSE)
  })
  set_dl("dl_quintile", function() { qd <- get_quintile_data(); if(is.null(qd)) data.frame() else qd }, "quintile_summary")

  output$census_area_table <- renderDT({
    cd <- get_census_geo_data(); d <- cd$combined |> filter(!is.na(n_children)); ml <- bll_mean_label(use_geo())
    if (!is.null(geo_county_lookup)) d <- d |> left_join(geo_county_lookup |> select(GEOID,County), by=c("geo_id"="GEOID")) else d$County <- NA
    d |> transmute(County, Area=geo_id, Tested=n_children, Elevated=n_elevated, `% Elevated`=round(pct_elevated,2),
                   !!paste("Avg BLL",ml):=round(mean_bll,2)) |> make_dt(10)
  })
  set_dl("dl_area_detail", function() {
    cd <- get_census_geo_data(); d <- cd$combined |> filter(!is.na(n_children))
    if (!is.null(geo_county_lookup)) d <- d |> left_join(geo_county_lookup |> select(GEOID,County), by=c("geo_id"="GEOID")) else d$County <- NA
    d |> select(County, geo_id, n_children, n_elevated, pct_elevated, mean_bll)
  }, "area_details")

  # ===== DATA TABLES =====

  make_geo_dt <- function(geo_col, id_name, show_county=FALSE, show_lhd=FALSE) {
    yr <- input$table_years; unit <- input$table_unit; dat <- active_data()
    ft <- dat |> filter(sample_year >= yr[1], sample_year <= yr[2])
    d <- summarize_lead(ft, geo_col, input$suppress_n, use_geo()); ml <- bll_mean_label(use_geo())
    geo_level_name <- switch(geo_col, lhd="LHD", county="County", tract_geoid="Tract", bg_geoid="Block Group", "Tract")
    d <- attach_rates(d, geo_level_name, yr)
    if (show_county && !is.null(geo_county_lookup)) {
      d <- d |> left_join(geo_county_lookup |> select(GEOID,County), by=c("geo_id"="GEOID"))
      d$county <- coalesce(d$County, d$county); d$County <- NULL
    }
    if (show_lhd) {
      if (!is.null(geo_lhd_lookup) && geo_col %in% c("tract_geoid","bg_geoid"))
        d <- d |> left_join(geo_lhd_lookup |> select(GEOID,lhd=LHD) |> distinct(GEOID,.keep_all=TRUE), by=c("geo_id"="GEOID"))
      else if (!is.null(lhd_map_df)) {
        d <- d |> left_join(lhd_map_df, by=c("geo_id"="county"))
        if (!"lhd" %in% names(d) && "lhd.y" %in% names(d)) d <- d |> rename(lhd=lhd.y) |> select(-any_of("lhd.x"))
      }
    }
    list(table=format_table(d, id_name, unit, show_county=show_county, show_lhd=show_lhd, mean_label=ml), data=d)
  }

  output$lhd_table    <- renderDT({ make_dt(make_geo_dt("lhd","LHD")$table, as.integer(input$table_rows)) })
  output$county_table <- renderDT({ make_dt(make_geo_dt("county","County",show_lhd=TRUE)$table, as.integer(input$table_rows)) })
  output$geo_table <- renderDT({
    g <- if(input$table_geo=="tract") "tract_geoid" else "bg_geoid"
    n <- if(input$table_geo=="tract") "Tract" else "Block Group"
    make_dt(make_geo_dt(g, n, show_county=TRUE, show_lhd=TRUE)$table, as.integer(input$table_rows))
  })

  output$addr_table <- renderDT({
    yr <- input$table_years; ml <- bll_mean_label(use_geo())
    astats <- get_address_stats() |> filter(last_year >= yr[1], first_year <= yr[2])
    if (use_geo()) {
      geo_means <- active_data() |> filter(sample_year >= yr[1], sample_year <= yr[2]) |>
        group_by(AddressIdentifier) |> summarize(mean_bll_geo=calc_mean_bll(result,TRUE), .groups="drop")
      astats <- astats |> left_join(geo_means, by="AddressIdentifier") |>
        mutate(mean_bll=coalesce(mean_bll_geo, mean_bll)) |> select(-mean_bll_geo)
    }
    if (!is.null(rv$full_preds)) {
      astats <- astats |> left_join(rv$full_preds |> select(AddressIdentifier,pred_prob=pred,single_child_flag) |>
                                      distinct(AddressIdentifier,.keep_all=TRUE), by="AddressIdentifier")
    } else { astats$pred_prob <- NA_real_; astats$single_child_flag <- FALSE }
    astats$single_child_flag[is.na(astats$single_child_flag)] <- FALSE
    tbl <- astats |> transmute(Address=AddressIdentifier, LHD=lhd, County=county,
      Tested=n_children, Elevated=n_elevated, `% Elevated`=round(100*n_elevated/pmax(n_children,1),2),
      !!paste("Avg BLL",ml):=round(mean_bll,2), `Max BLL`=round(max_bll,2),
      `First Year`=first_year, `Last Year`=last_year,
      `Pred Risk`=if_else(!is.na(pred_prob), paste0(sprintf("%.1f%%",pred_prob*100),if_else(single_child_flag,"*","")),"\u2014"),
      `View`=paste0('<button class="btn btn-xs btn-info act-view-tests-addr" data-addr="',htmltools::htmlEscape(AddressIdentifier,TRUE),'">View</button>'),
      `Risk`=paste0('<button class="btn btn-xs btn-warning act-risk-breakdown" data-addr="',htmltools::htmlEscape(AddressIdentifier,TRUE),'"',
                    if_else(!is.na(pred_prob),'>',' disabled>'),if_else(!is.na(pred_prob),'Explain','\u2014'),'</button>'))
    datatable(tbl, escape=FALSE, options=list(pageLength=as.integer(input$table_rows), scrollX=TRUE), rownames=FALSE)
  })
  observeEvent(input$view_tests_addr_click, show_address_tests(input$view_tests_addr_click))

  # Download handlers
  set_dl("dl_lhd_disp", function() make_geo_dt("lhd","LHD")$table, "lhd_data")
  set_dl("dl_lhd_all", function() make_geo_dt("lhd","LHD")$table, "lhd_all")
  set_dl("dl_county_disp", function() make_geo_dt("county","County",show_lhd=TRUE)$table, "county_data")
  set_dl("dl_county_all", function() make_geo_dt("county","County",show_lhd=TRUE)$table, "county_all")
  set_dl("dl_geo_disp", function() {
    g <- if(input$table_geo=="tract") "tract_geoid" else "bg_geoid"
    n <- if(input$table_geo=="tract") "Tract" else "Block Group"
    make_geo_dt(g, n, show_county=TRUE, show_lhd=TRUE)$table
  }, "geo_data")
  set_dl("dl_geo_all", function() {
    g <- if(input$table_geo=="tract") "tract_geoid" else "bg_geoid"
    n <- if(input$table_geo=="tract") "Tract" else "Block Group"
    make_geo_dt(g, n, show_county=TRUE, show_lhd=TRUE)$table
  }, "geo_all")
  set_dl("dl_addr_disp", function() {
    get_address_stats() |> filter(last_year>=input$table_years[1], first_year<=input$table_years[2]) |>
      select(AddressIdentifier, lhd, county, n_children, n_elevated, mean_bll, max_bll, first_year, last_year)
  }, "addresses")
  set_dl("dl_addr_all", function() {
    get_address_stats() |> select(AddressIdentifier, lhd, county, n_children, n_elevated, mean_bll, max_bll, first_year, last_year)
  }, "addresses_all")

  # ===== RISK MODEL =====

  # Risk breakdown modal
  show_risk_breakdown <- function(addr_id) {
    fit <- rv$fit; preds <- rv$full_preds
    if (is.null(fit) || is.null(preds)) { showNotification("Run the Risk Model first.", type="warning"); return() }
    row <- preds |> filter(AddressIdentifier == addr_id)
    if (nrow(row) == 0) { showNotification("Address not found.", type="warning"); return() }
    row <- row[1, ]; is_sc <- isTRUE(row$single_child_flag)
    use_re <- inherits(fit, "glmerMod")
    coefs <- if (use_re) fixef(fit) else coef(fit)
    intercept <- coefs["(Intercept)"]; beta <- coefs[names(coefs) != "(Intercept)"]
    fac_vars <- c("income_cat","first_sex","first_sample_type","acs_year","first_year_cat","first_age_yr_cat")
    train_means <- if (!is.null(rv$train_sample)) {
      sapply(names(beta), function(nm) {
        bv <- nm; for (fv in fac_vars) if (startsWith(nm,fv) && nchar(nm)>nchar(fv)) { bv <- fv; break }
        if (bv %in% names(rv$train_sample)) mean(as.numeric(rv$train_sample[[bv]]), na.rm=TRUE) else NA_real_
      })
    } else setNames(rep(NA_real_, length(beta)), names(beta))
    contrib_rows <- list(); running <- intercept
    for (nm in names(beta)) {
      b <- beta[nm]; bv <- nm; val <- NA_real_; vd <- ""; fmatch <- FALSE
      for (fv in fac_vars) {
        if (startsWith(nm,fv) && nchar(nm)>nchar(fv)) {
          bv <- fv; lv <- substring(nm, nchar(fv)+1); av <- as.character(row[[fv]]); active <- (av==lv)
          val <- if(active) 1 else 0; vd <- paste0(av, if(active) paste0(" (= ",lv,")") else paste0(" (\u2260 ",lv,")"))
          fmatch <- TRUE; break
        }
      }
      if (!fmatch) { if (nm %in% names(row) && !is.na(row[[nm]])) { val <- as.numeric(row[[nm]]); vd <- round(val,3) } else { val <- 0; vd <- "N/A" } }
      contrib <- b*val; running <- running + contrib
      avg_v <- train_means[nm]; dfa <- if(!is.na(avg_v) && !fmatch) round(val-avg_v,3) else NA_real_
      ab <- if(!is.na(dfa)) { if(dfa>0.001) "\u25b2 Above" else if(dfa< -0.001) "\u25bc Below" else "\u2248 At Avg" } else "\u2014"
      contrib_rows[[length(contrib_rows)+1]] <- list(name=get_label(bv), value=vd,
        diff=if(!is.na(dfa)) sprintf("%+.3f",dfa) else "\u2014", contribution=round(contrib,4),
        direction=if(contrib>0.005) "\u2191 increases" else if(contrib< -0.005) "\u2193 decreases" else "\u2194 minimal", ab=ab)
    }
    contrib_rows <- contrib_rows[order(-sapply(contrib_rows, \(x) abs(x$contribution)))]
    re_text <- NULL
    if (use_re) {
      rev <- names(ranef(fit))[1]; rval <- row[[rev]]
      if (!is.null(rval) && !is.na(rval)) {
        re_adj <- tryCatch(ranef(fit)[[rev]][as.character(rval),1], error=function(e) 0)
        if (!is.na(re_adj) && re_adj != 0) {
          running <- running + re_adj
          re_text <- div(style="background:#f0f7ff;padding:12px;border-radius:6px;margin:10px 0;",
            p(strong("\u2699 Neighborhood:"), sprintf("%s %s adjusts by %+.3f", rev, rval, re_adj), style="margin:0;font-size:13px;"))
        }
      }
    }
    fp <- 1/(1+exp(-running)); fpct <- round(fp*100,1)
    rc <- if(fpct>=30) "#d73027" else if(fpct>=15) "#fc8d59" else if(fpct>=5) "#fee08b" else "#91cf60"
    rl <- if(fpct>=30) "High" else if(fpct>=15) "Moderate" else if(fpct>=5) "Low-Moderate" else "Low"
    th <- '<table style="width:100%;border-collapse:collapse;font-size:13px;margin:14px 0;"><tr style="background:#f5f5f5;font-weight:bold;border-bottom:2px solid #ddd;"><td style="padding:10px 12px;">Factor</td><td style="padding:10px 12px;">This Address</td><td style="padding:10px 12px;">Diff</td><td style="padding:10px 12px;text-align:right;">Impact</td><td style="padding:10px 12px;">Direction</td><td style="padding:10px 12px;">vs Avg</td></tr>'
    tr <- ""
    for (cr in contrib_rows) {
      bg <- if(cr$contribution>0.005) "#fff5f5" else if(cr$contribution< -0.005) "#f0fff0" else "#fff"
      ac <- if(grepl("\u25b2",cr$ab)) "#d73027" else if(grepl("\u25bc",cr$ab)) "#2c7fb8" else "#888"
      tr <- paste0(tr,'<tr style="border-bottom:1px solid #eee;background:',bg,';"><td style="padding:9px 12px;">',cr$name,
        '</td><td style="padding:9px 12px;">',cr$value,'</td><td style="padding:9px 12px;font-family:monospace;">',cr$diff,
        '</td><td style="padding:9px 12px;text-align:right;font-weight:bold;">',sprintf("%+.4f",cr$contribution),
        '</td><td style="padding:9px 12px;font-size:12px;">',cr$direction,'</td><td style="padding:9px 12px;color:',ac,';font-weight:600;font-size:12px;">',cr$ab,'</td></tr>')
    }
    showModal(modalDialog(title=paste0("Risk Breakdown: ",addr_id), size="l", easyClose=TRUE,
      div(style="background:#eef6ff;padding:14px;border-radius:8px;margin-bottom:14px;",
        p(strong("What is this?"), style="margin:0 0 5px 0;"),
        p("Shows how the model estimated this address's risk.", style="margin:0;font-size:13px;")),
      if(is_sc) div(style="background:#fff8e1;border:1px solid #ffc107;padding:12px;border-radius:6px;margin-bottom:14px;",
        p(HTML("<b>\u26a0 Single-child address:</b> predicted risk shown but not in validation."), style="margin:0;font-size:12px;")) else NULL,
      div(style=paste0("text-align:center;padding:18px;border-radius:8px;border:2px solid ",rc,";margin-bottom:16px;"),
        h3(paste0(fpct,"% estimated risk"), style=paste0("color:",rc,";margin:0 0 4px 0;")),
        p(HTML(paste0("Risk Level: <b>",rl,"</b>")), style="margin:0;font-size:14px;")),
      h5("Factor Contributions:", style="color:#2c7fb8;margin-top:16px;"),
      p("Baseline:", strong(sprintf("%.3f",intercept)), "| Red rows increase risk; green decrease.", style="font-size:13px;"),
      HTML(paste0(th,tr,'</table>')), if(!is.null(re_text)) re_text else NULL,
      div(style="background:#f9f9f9;padding:12px;border-radius:6px;margin:12px 0;",
        p(HTML(paste0("Score: <b>",sprintf("%.3f",running),"</b> \u2192 Risk = <b>",fpct,"%</b>")), style="font-size:13px;margin:0;text-align:center;")),
      footer=modalButton("Close")))
  }
  observeEvent(input$risk_breakdown_click, show_risk_breakdown(input$risk_breakdown_click))

  # --- Model fitting (broken into phases) ---

  prepare_model_inputs <- function() {
    re <- input$random_effect; use_re <- re != "none"
    acs_mode <- input$acs_mode %||% "nof"; use_max <- isTRUE(input$test_selection == "max")
    use_all_addr <- isTRUE(input$model_addr_pop == "all")
    geo_key <- if (re == "bg_geoid") "bg" else "tract"
    geo_col <- if (geo_key=="bg") "bg_geoid" else "tract_geoid"
    census_ts <- if (geo_key=="bg") bg_census_ts else tract_census_ts
    endpoint_acs_yr <- if (acs_mode=="endpoint") acs_year_map(max(years)) else NULL

    md <- resolve_model_data(acs_mode, geo_key, use_max, use_all_addr)
    if (acs_mode == "endpoint") {
      base_md <- resolve_model_data("base", geo_key, use_max, use_all_addr)
      if (is.null(base_md)) { showNotification("Model data not loaded", type="error"); return(NULL) }
      train_end_acs <- acs_year_map(input$train_years[2])
      test_end_acs <- acs_year_map(max(years))
      if (is.null(census_ts)) { showNotification("Census TS not loaded", type="error"); return(NULL) }
      if (input$train_all) {
        md <- base_md |> left_join(census_ts |> filter(acs_year==train_end_acs) |> select(-acs_year), by=geo_col)
      } else {
        tp <- base_md |> filter(first_year>=input$train_years[1], first_year<=input$train_years[2]) |>
          left_join(census_ts |> filter(acs_year==train_end_acs) |> select(-acs_year), by=geo_col)
        tt <- base_md |> filter(first_year>input$train_years[2]) |>
          left_join(census_ts |> filter(acs_year==test_end_acs) |> select(-acs_year), by=geo_col)
        md <- bind_rows(tp, tt)
      }
      showNotification(paste0("Endpoint: train ACS ",train_end_acs,", test ACS ",test_end_acs), type="warning", duration=6)
    }
    if (is.null(md)) { showNotification("Model data not loaded", type="error"); return(NULL) }

    # Max-BLL leakage fix
    if (use_max && !input$train_all) {
      train_end <- input$train_years[2]
      train_feat <- first_test |> filter(sample_year<=train_end) |> rebuild_addr_features("max") |>
        filter(!is.na(.data[[geo_col]]), nchar(.data[[geo_col]])>=11)
      test_feat <- first_test |> filter(sample_year>train_end) |> rebuild_addr_features("max") |>
        filter(!is.na(.data[[geo_col]]), nchar(.data[[geo_col]])>=11)
      if (!is.null(census_ts)) {
        train_feat <- join_census_for_mode(train_feat, acs_mode, geo_col, census_ts, endpoint_acs_year=endpoint_acs_yr)
        test_feat <- join_census_for_mode(test_feat, acs_mode, geo_col, census_ts, endpoint_acs_year=endpoint_acs_yr)
      }
      md <- bind_rows(train_feat, test_feat)
      showNotification(paste0("Max-BLL rebuilt: Train ",nrow(train_feat)," | Test ",nrow(test_feat)), type="warning", duration=6)
    }

    list(md=md, re=re, use_re=use_re, acs_mode=acs_mode, use_max=use_max,
         use_all_addr=use_all_addr, geo_key=geo_key, geo_col=geo_col,
         census_ts=census_ts, endpoint_acs_yr=endpoint_acs_yr)
  }

  do_fit_model <- function() {
    showNotification("Fitting model...", id="model_note", type="message", duration=NULL)
    on.exit(removeNotification("model_note"))

    inp <- prepare_model_inputs(); if (is.null(inp)) return(NULL)
    md <- inp$md; re <- inp$re; use_re <- inp$use_re; acs_mode <- inp$acs_mode
    use_all_addr <- inp$use_all_addr

    all_covars <- get_all_covars()
    # Handle categorical year/age
    if ("first_year_cat" %in% all_covars) { md$first_year_cat <- factor(md$first_year); all_covars <- unique(all_covars) }
    if ("first_age_yr_cat" %in% all_covars) { md$first_age_yr_cat <- factor(pmin(floor(md$first_age_mo/12),6)); all_covars <- unique(all_covars) }
    if (length(all_covars) == 0) { showNotification("Select at least one variable", type="error"); return(NULL) }

    fac_covars <- intersect(all_covars, c("income_cat","first_sex","first_sample_type","first_year_cat","first_age_yr_cat"))
    num_covars <- setdiff(all_covars, fac_covars)

    d <- md
    if (length(num_covars) > 0) d <- d |> filter(if_all(all_of(num_covars), ~ !is.na(.) & is.finite(.)))
    if (length(fac_covars) > 0) d <- d |> filter(if_all(all_of(fac_covars), ~ !is.na(.)))

    # Split-aware outcome computation
    train_yrs <- input$train_years
    if (input$train_all) {
      outcomes <- first_test |>
        inner_join(d |> select(AddressIdentifier, first_id), by="AddressIdentifier") |>
        filter(PATIENT_LOCAL_ID != first_id) |>
        summarize(outcome=as.integer(any(elevated==1, na.rm=TRUE)), .by=AddressIdentifier)
      d <- d |> select(-any_of("outcome")) |> left_join(outcomes, by="AddressIdentifier") |>
        filter(!is.na(outcome))
      train <- d; test <- NULL
    } else {
      train_end <- train_yrs[2]; d <- d |> select(-any_of("outcome"))
      train_out <- first_test |>
        inner_join(d |> select(AddressIdentifier, first_id, first_year), by="AddressIdentifier") |>
        filter(PATIENT_LOCAL_ID!=first_id, sample_year<=train_end) |>
        summarize(outcome=as.integer(any(elevated==1, na.rm=TRUE)), .by=AddressIdentifier)
      test_out <- first_test |>
        inner_join(d |> select(AddressIdentifier, first_id), by="AddressIdentifier") |>
        filter(PATIENT_LOCAL_ID!=first_id, sample_year>train_end) |>
        summarize(outcome=as.integer(any(elevated==1, na.rm=TRUE)), .by=AddressIdentifier)
      d_train <- d |> filter(first_year>=train_yrs[1], first_year<=train_end) |>
        left_join(train_out, by="AddressIdentifier") |> filter(!is.na(outcome))
      d_test <- d |> filter(first_year>train_end) |>
        left_join(test_out, by="AddressIdentifier") |> filter(!is.na(outcome))
      d <- bind_rows(d_train, d_test); train <- d_train; test <- d_test
      showNotification(paste0("Split at ",train_end,": Train ",nrow(train)," (",round(100*mean(train$outcome),1),
        "%) | Test ",nrow(test)," (",round(100*mean(test$outcome),1),"%)"), type="message", duration=6)
    }
    if (nrow(train) < 50) { showNotification("Too few training records", type="error"); return(NULL) }
    rv$train_sample <- train; rv$test_sample <- test

    if (isTRUE(input$use_log_bll) && "first_bll" %in% all_covars) {
      logtr <- function(df) { df$first_bll <- log(pmax(df$first_bll,0.1)); df }
      train <- logtr(train); if(!is.null(test)) test <- logtr(test); d <- logtr(d)
    }

    # Scaling
    scaled_vars <- character(0); scale_params <- list()
    never_scale <- fac_covars
    if (input$scale_vars && length(num_covars) > 0) {
      for (v in setdiff(num_covars, never_scale)) {
        m <- mean(train[[v]], na.rm=TRUE); s <- sd(train[[v]], na.rm=TRUE)
        if (!is.na(s) && s > 0) {
          scale_params[[v]] <- list(mean=m, sd=s)
          train[[v]] <- (train[[v]]-m)/s; if(!is.null(test)) test[[v]] <- (test[[v]]-m)/s; d[[v]] <- (d[[v]]-m)/s
          scaled_vars <- c(scaled_vars, v)
        }
      }
    }
    rv$is_scaled <- length(scaled_vars) > 0
    covar_str <- paste(all_covars, collapse=" + ")
    f <- if (use_re) as.formula(paste0("outcome ~ ",covar_str," + (1|",re,")"))
         else as.formula(paste0("outcome ~ ",covar_str))

    fit <- suppressWarnings(tryCatch({
      if (use_re) glmer(f, data=train, family=binomial, nAGQ=1, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=50000)))
      else glm(f, data=train, family=binomial)
    }, error=function(e) { showNotification(paste("Model error:",e$message), type="error"); NULL }))
    if (is.null(fit)) return(NULL)
    if (use_re && inherits(fit,"glmerMod") && isSingular(fit))
      showNotification("Warning: singular model (RE variance near zero).", type="warning", duration=8)

    rv$fit <- fit; rv$model_run <- TRUE

    # GLM comparison
    rv$glm_fit <- NULL; rv$glm_test_preds <- NULL; rv$glm_roc <- NULL
    tryCatch({
      gf <- glm(as.formula(paste0("outcome ~ ",covar_str)), data=train, family=binomial)
      rv$glm_fit <- gf
      if (!is.null(test) && nrow(test) > 0) { rv$glm_test_preds <- predict(gf, newdata=test, type="response"); rv$glm_roc <- compute_roc(test$outcome, rv$glm_test_preds) }
    }, error=function(e) message("[GLM] ",e$message))

    if (!is.null(test) && nrow(test) > 0) {
      test$pred <- predict(fit, newdata=test, type="response", allow.new.levels=if(use_re) TRUE else NULL)
      rv$test_preds <- test; rv$roc_obj <- compute_roc(test$outcome, test$pred)
    }
    d$pred <- predict(fit, newdata=d, type="response", allow.new.levels=if(use_re) TRUE else NULL)
    d$single_child_flag <- FALSE; rv$full_preds <- d

    # Single-child predictions (2+ mode only)
    if (!use_all_addr) {
    tryCatch({
      existing <- unique(d$AddressIdentifier)
      single_addr <- address_stats_all |> filter(!AddressIdentifier %in% existing, n_children==1)
      if (nrow(single_addr) > 0) {
        sdf <- join_census_for_mode(single_addr, acs_mode, inp$geo_col, inp$census_ts, endpoint_acs_year=inp$endpoint_acs_yr)
        if (is.null(sdf)) sdf <- single_addr
        for (mc in setdiff(all_covars, names(sdf))) sdf[[mc]] <- NA_real_
        if (length(num_covars)>0) sdf <- sdf |> filter(if_all(all_of(num_covars), ~!is.na(.)&is.finite(.)))
        if (length(fac_covars)>0) sdf <- sdf |> filter(if_all(all_of(fac_covars), ~!is.na(.)))
        if (nrow(sdf) > 0) {
          if (isTRUE(input$use_log_bll) && "first_bll" %in% all_covars) sdf$first_bll <- log(pmax(sdf$first_bll,0.1))
          for (v in scaled_vars) if (v %in% names(scale_params) && v %in% names(sdf)) {
            m <- scale_params[[v]]$mean; s <- scale_params[[v]]$sd
            if (!is.na(s) && s>0) sdf[[v]] <- (sdf[[v]]-m)/s
          }
          if (use_re && !re %in% names(sdf)) sdf[[re]] <- NA_character_
          sdf$outcome <- NA_integer_
          sdf$pred <- predict(fit, newdata=sdf, type="response", allow.new.levels=TRUE)
          sdf$single_child_flag <- TRUE
          shared <- intersect(names(rv$full_preds), names(sdf))
          rv$full_preds <- bind_rows(rv$full_preds[shared], sdf[shared])
        }
      }
    }, error=function(e) message("[Single-child] ",e$message))
    }

    # Correlation matrices
    if (length(num_covars) >= 2) {
      rv$cor_matrix_pearson <- cor(train[num_covars], use="pairwise.complete.obs", method="pearson")
      rv$cor_matrix_spearman <- cor(train[num_covars], use="pairwise.complete.obs", method="spearman")
      rv$cor_matrix <- rv$cor_matrix_pearson
    } else { rv$cor_matrix <- NULL; rv$cor_matrix_pearson <- NULL; rv$cor_matrix_spearman <- NULL }

    # Comparison stats
    if (!is.null(test) && nrow(test) > 0) {
      n_test <- nrow(test); n_elev <- sum(test$outcome==1); screen_n <- round(n_test*0.5)
      test_ranked <- test |> arrange(desc(pred))
      model_det <- sum(test_ranked$outcome[1:screen_n]==1); random_det <- round(n_elev*0.5)
      rv$comparison <- list(n_test=n_test, n_elev=n_elev, screen_n=screen_n,
        model_det=model_det, random_det=random_det,
        model_pct=round(100*model_det/max(n_elev,1),1), random_pct=50,
        improve=round(100*(model_det-random_det)/max(random_det,1),1),
        auc=round(rv$roc_obj$auc,3),
        tests_per_pos_model=round(screen_n/max(model_det,1),1),
        tests_per_pos_random=round(screen_n/max(random_det,1),1))
    }

    # XGBoost comparison
    rv$xgb_fit <- NULL; rv$xgb_roc <- NULL; rv$xgb_comparison <- NULL; rv$xgb_test_preds <- NULL
    if (!isTRUE(rv$autorun_active)) {
      showNotification("Fitting XGBoost...", id="xgb_note", type="message", duration=NULL)
      tryCatch({
        xf <- as.formula(paste0("~ ",covar_str))
        xt <- model.matrix(xf, data=train)[,-1,drop=FALSE]; dtrain <- xgb.DMatrix(data=xt, label=train$outcome)
        if (!is.null(test) && nrow(test)>0) {
          xte <- model.matrix(xf, data=test)[,-1,drop=FALSE]; dtest <- xgb.DMatrix(data=xte, label=test$outcome)
          wl <- list(train=dtrain, test=dtest)
        } else { wl <- list(train=dtrain); dtest <- NULL }
        xgb_fit <- xgb.train(params=list(objective="binary:logistic",eval_metric="auc",max_depth=4,eta=0.1,subsample=0.8),
                              data=dtrain, nrounds=200, watchlist=wl, early_stopping_rounds=20, verbose=0)
        rv$xgb_fit <- xgb_fit
        if (!is.null(test) && nrow(test)>0) {
          xp <- predict(xgb_fit, dtest); rv$xgb_test_preds <- xp; rv$xgb_roc <- compute_roc(test$outcome, xp)
          xord <- order(xp, decreasing=TRUE); n_half <- ceiling(nrow(test)/2)
          rv$xgb_comparison <- list(glmer_auc=rv$roc_obj$auc, xgb_auc=rv$xgb_roc$auc,
            glmer_capture=rv$comparison$model_pct/100, xgb_capture=sum(test$outcome[xord[1:n_half]])/sum(test$outcome),
            n_test=nrow(test), importance=xgb.importance(model=xgb_fit) |> head(10))
        }
      }, error=function(e) showNotification(paste("XGBoost:",e$message), type="warning"))
      removeNotification("xgb_note")
    }
    showNotification("Model fitting complete!", type="message", duration=3)
    return(fit)
  }

  observeEvent(input$fit_model, {
    tryCatch(do_fit_model(), error=function(e) showNotification(paste("Error:",e$message), type="error", duration=10))
  }, ignoreInit=TRUE)

  rv$autorun_active <- FALSE
  session$onFlushed(function() {
    rv$autorun_active <- TRUE
    tryCatch(do_fit_model(), error=function(e) showNotification(paste("Auto-run:",e$message), type="warning"))
    rv$autorun_active <- FALSE
  }, once=TRUE)

  # --- Model outputs ---

  output$interpretation <- renderUI({
    fit <- rv$fit
    if (is.null(fit)) return(div(class="metric-box", p("Click 'Run Model' to fit.")))
    is_glmer <- inherits(fit, "glmerMod")
    coefs <- tryCatch({
      raw <- if(is_glmer) tidy(fit, effects="fixed", conf.int=TRUE) else tidy(fit, conf.int=TRUE)
      raw |> filter(term!="(Intercept)") |> mutate(or=exp(estimate), var_name=gsub("scale\\(|\\)","",term))
    }, error=function(e) data.frame())
    is_scaled <- rv$is_scaled
    never_sv <- c("income_cat","first_sex","first_sample_type","acs_year",
                   paste0("re_",c("hispanic","white_nh","black_nh","aian_nh","asian_nh","nhopi_nh","multi_nh","other_nh")))
    findings <- if(nrow(coefs)>0) {
      lapply(1:nrow(coefs), function(i) {
        r <- coefs[i,]; pct <- round((r$or-1)*100,0)
        is_fac <- any(sapply(c("income_cat","first_sex","first_sample_type","acs_year"), \(f) startsWith(r$var_name,f)))
        is_ind <- r$var_name %in% paste0("re_",c("hispanic","white_nh","black_nh","aian_nh","asian_nh","nhopi_nh","multi_nh","other_nh"))
        is_yr <- grepl("year",r$var_name,ignore.case=TRUE)
        vu <- if(is_fac) "being in this category" else if(is_ind) "being in this group"
              else if(is_yr) "a 1-year" else if(is_scaled && !r$var_name %in% never_sv) "a 1 SD" else "a 1-unit"
        tags$li(strong(get_label(r$var_name)),": ",
          paste0(vu,if(!is_fac&&!is_ind) " increase" else ""," \u2192 ",abs(pct),"% ",
                 if(r$or>1)"higher" else "lower"," odds (OR=",round(r$or,2),")"))
      })
    } else list(tags$li("Geographic grouping only."))
    auc_val <- rv$comparison$auc
    auc_text <- if(!is.null(auc_val)) tagList(
      p(strong("Multilevel AUC = ",auc_val)),
      p(if(auc_val>=0.8)"Good discrimination." else if(auc_val>=0.7)"Acceptable." else "Modest \u2014 consider more predictors."))
    xgb_text <- if(!is.null(rv$xgb_roc)) {
      xa <- round(rv$xgb_roc$auc,3); d <- round(xa-(auc_val%||%0.5),3)
      w <- if(d>0.01)"outperforms" else if(d< -0.01)"underperforms" else "matches"
      tagList(p(strong("XGBoost AUC = ",xa)), p(paste0("XGBoost ",w," Multilevel by ",abs(d),".")))
    }
    glm_text <- if(!is.null(rv$glm_roc)) {
      ga <- round(rv$glm_roc$auc,3); d <- round(ga-(auc_val%||%0.5),3)
      w <- if(d>0.01)"outperforms" else if(d< -0.01)"underperforms" else "matches"
      tagList(p(strong("Logistic AUC = ",ga)), p(paste0("Logistic ",w," Multilevel by ",abs(d),".")))
    }
    div(class="metric-box", h4("Findings:"), tags$ul(findings),
      if(is_scaled) p(em("ORs are standardized (1 SD change)."), style="font-size:11px;color:#666;"),
      hr(), h4("Performance:"), auc_text, glm_text, xgb_text)
  })

  output$or_note <- renderUI({
    if(is.null(rv$fit)) return(NULL)
    div(class="warning-box", if(rv$is_scaled) "Standardized ORs (1 SD for continuous). Categorical/time vars show raw ORs."
        else "Unstandardized ORs (1-unit). Magnitudes not directly comparable across scales.")
  })

  output$coef_plot <- renderPlotly({
    or_model <- input$or_model_select %||% "multilevel"
    fit <- if(or_model=="glm"&&!is.null(rv$glm_fit)) rv$glm_fit else rv$fit; if(is.null(fit)) return(NULL)
    ig <- inherits(fit,"glmerMod")
    coefs <- tryCatch({
      (if(ig) tidy(fit,effects="fixed",conf.int=TRUE) else tidy(fit,conf.int=TRUE)) |>
        filter(term!="(Intercept)") |>
        mutate(or=exp(estimate), or_lo=exp(conf.low), or_hi=exp(conf.high),
               var_name=gsub("scale\\(|\\)","",term), label=sapply(var_name,get_label),
               Effect=ifelse(or>1.5,"Strong +",ifelse(or>1.01,"Moderate +",ifelse(or<0.67,"Strong -",ifelse(or<0.99,"Moderate -","Weak")))))
    }, error=function(e) data.frame()); if(nrow(coefs)==0) return(NULL)
    p <- ggplot(coefs, aes(or, reorder(label,or), text=paste0(label,"\nOR: ",round(or,2),"\nCI: ",round(or_lo,2),"-",round(or_hi,2)))) +
      geom_vline(xintercept=1, linetype="dashed", color="gray50") +
      geom_errorbar(aes(xmin=or_lo, xmax=or_hi), width=0.2, color="gray40", orientation="y") +
      geom_point(aes(color=Effect), size=4) +
      scale_color_manual(values=c("Strong +"="#d73027","Moderate +"="#fc8d59","Weak"="#999",
                                  "Moderate -"="#91bfdb","Strong -"="#4575b4")) +
      labs(title="Odds Ratios (95% CI)", x="Odds Ratio", y=NULL) + theme_minimal() + theme(legend.position="bottom")
    ggplotly(p, tooltip="text") |> layout(legend=list(orientation="h",y=-0.2))
  })

  output$cor_table <- renderDT({
    cm <- if(isTRUE(input$cor_method=="spearman")) rv$cor_matrix_spearman else rv$cor_matrix_pearson
    if(is.null(cm)) cm <- rv$cor_matrix; if(is.null(cm)) return(NULL)
    cd <- as.data.frame(round(cm,3)); cd$Variable <- sapply(rownames(cd), get_label)
    cd <- cd |> select(Variable, everything()); names(cd)[-1] <- sapply(names(cd)[-1], get_label)
    datatable(cd, options=list(dom="t",pageLength=15,scrollX=TRUE), rownames=FALSE) |>
      formatStyle(columns=2:ncol(cd), backgroundColor=styleInterval(c(-0.7,-0.5,0.5,0.7,0.99),
                                                                     c('#f8d7da','#fff3cd','white','#fff3cd','#f8d7da','white')))
  })

  output$vif_table <- renderDT({
    cm <- if(isTRUE(input$cor_method=="spearman")) rv$cor_matrix_spearman else rv$cor_matrix_pearson
    if(is.null(cm)) cm <- rv$cor_matrix; if(is.null(cm)||ncol(cm)<2) return(NULL)
    tryCatch({
      inv <- solve(cm); vifs <- diag(inv)
      vdf <- data.frame(Variable=sapply(names(vifs),get_label), VIF=round(vifs,2)) |> arrange(desc(VIF))
      datatable(vdf, options=list(dom="t",pageLength=15), rownames=FALSE) |>
        formatStyle("VIF", backgroundColor=styleInterval(c(5,10), c("white","#fff3cd","#f8d7da")))
    }, error=function(e) NULL)
  })

  output$model_summary <- renderPrint({ if(is.null(rv$fit)) cat("No model fitted") else summary(rv$fit) })

  # Model comparison UI
  output$model_comparison_ui <- renderUI({
    rc <- rv$comparison; if(is.null(rc)) return(p("Run model first."))
    has_xgb <- !is.null(rv$xgb_roc); has_glm <- !is.null(rv$glm_roc)
    ma <- rc$auc; mc <- rc$model_pct; md <- rc$model_det; mt <- rc$tests_per_pos_model
    ga <- gc <- gd <- gt <- NA
    if(has_glm && !is.null(rv$glm_test_preds) && !is.null(rv$test_preds)) {
      ga <- round(rv$glm_roc$auc,3); go <- order(rv$glm_test_preds,decreasing=TRUE)
      gd <- sum(rv$test_preds$outcome[go[1:rc$screen_n]]==1); gc <- 100*gd/max(rc$n_elev,1)
      gt <- round(rc$screen_n/max(gd,1),1)
    }
    xa <- xc <- xd <- xt <- NA
    if(has_xgb && !is.null(rv$xgb_comparison)) {
      xa <- round(rv$xgb_roc$auc,3); xc <- rv$xgb_comparison$xgb_capture*100
      if(!is.null(rv$xgb_test_preds)&&!is.null(rv$test_preds)) {
        xo <- order(rv$xgb_test_preds,decreasing=TRUE)
        xd <- sum(rv$test_preds$outcome[xo[1:rc$screen_n]]==1); xt <- round(rc$screen_n/max(xd,1),1)
      }
    }
    mp <- function(nm,col,bg,det,cap,tpp)
      div(style=paste0("text-align:center;padding:15px;background:",bg,";border-radius:6px;"),
        h5(nm,style=paste0("color:",col,";margin:0;")),
        h2(if(!is.na(det)) format(det,big.mark=",") else "N/A", style=paste0("margin:10px 0;color:",col,";")),
        p("elevated found"), p(strong(if(!is.na(cap)) sprintf("%.1f%%",cap) else "N/A")," of all elevated"),
        hr(), p(strong(if(!is.na(tpp)) tpp else "N/A")," tests per positive"))
    ac <- function(nm,col,bg,auc,cap)
      div(style=paste0("background:",bg,";padding:14px;border-radius:6px;border-left:4px solid ",col,";margin-bottom:10px;"),
        h4(nm,style=paste0("margin-top:0;color:",col,";")),
        p(strong("AUC:"),if(!is.na(auc))auc else"N/A"), p(strong("Capture@50%:"),if(!is.na(cap))sprintf("%.1f%%",cap)else"N/A"))
    tagList(
      div(class="compare-box", h4("Screening Efficiency: 50% Screen Rate", style="margin-top:0;"),
          p("Testing ",format(rc$screen_n,big.mark=",")," of ",format(rc$n_test,big.mark=",")," addresses:"),
          hr(),
          fluidRow(column(6, mp("Multilevel","#2c7fb8","#e8f4f8",md,mc,mt)),
                   column(6, mp("Logistic","#1a9850","#f0f7f0",gd,gc,gt))),
          fluidRow(column(6, mp("XGBoost","#d73027","#fff5f5",if(!is.na(xd))xd else NA,xc,xt)),
                   column(6, mp("Random","#666","#f5f5f5",rc$random_det,50.0,rc$tests_per_pos_random)))),
      hr(),
      fluidRow(column(6, ac("Multilevel","#2c7fb8","#f0f7fb",ma,mc)),
               column(6, ac("Logistic","#1a9850","#f0f7f0",ga,gc))),
      fluidRow(column(6, ac("XGBoost","#d73027","#fff5f5",xa,xc)),
               column(6, ac("Random","#999","#f5f5f5",0.5,50.0))),
      if(has_xgb && !is.null(rv$xgb_comparison$importance) && nrow(rv$xgb_comparison$importance)>0)
        tagList(hr(), h5("XGBoost Top Features"), DTOutput("xgb_importance_table")))
  })
  output$xgb_importance_table <- renderDT({
    comp <- rv$xgb_comparison; if(is.null(comp)||is.null(comp$importance)) return(NULL)
    datatable(comp$importance |> mutate(across(where(is.numeric), ~round(.,4))), options=list(dom="t",pageLength=10), rownames=FALSE)
  })

  # Confusion Matrix
  output$confusion_matrix_ui <- renderUI({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(p("Run model first."))
    ccm <- function(labels, preds, roc_obj) {
      if(is.null(roc_obj)) return(NULL)
      yi <- which.max(roc_obj$sensitivities+roc_obj$specificities-1)
      th <- roc_obj$thresholds[yi]; if(!is.finite(th)) th <- 0.5
      pc <- as.integer(preds>=th)
      tp <- sum(pc==1&labels==1); fp <- sum(pc==1&labels==0); fn <- sum(pc==0&labels==1); tn <- sum(pc==0&labels==0)
      list(tp=tp,fp=fp,fn=fn,tn=tn,sens=tp/(tp+fn),spec=tn/(tn+fp),acc=(tp+tn)/(tp+fp+fn+tn),
           ppv=if((tp+fp)>0)tp/(tp+fp) else NA,npv=if((tn+fn)>0)tn/(tn+fn) else NA,threshold=round(th,4))
    }
    cmp <- function(nm,col,cm) {
      if(is.null(cm)) return(NULL)
      div(style=paste0("border:2px solid ",col,";border-radius:8px;padding:14px;margin-bottom:12px;"),
        h5(nm,style=paste0("color:",col,";margin-top:0;")),
        p(em(paste0("Threshold: ",cm$threshold)), style="font-size:11px;color:#666;margin-bottom:8px;"),
        tags$table(style="width:100%;border-collapse:collapse;text-align:center;font-size:13px;margin-bottom:10px;",
          tags$tr(tags$th(""), tags$th("Pred +",style="padding:6px;background:#f5f5f5;"), tags$th("Pred -",style="padding:6px;background:#f5f5f5;")),
          tags$tr(tags$td(strong("Act +"),style="padding:6px;background:#f5f5f5;"),
                  tags$td(format(cm$tp,big.mark=","),style="padding:6px;background:#e8f5e9;font-weight:bold;"),
                  tags$td(format(cm$fn,big.mark=","),style="padding:6px;background:#ffebee;")),
          tags$tr(tags$td(strong("Act -"),style="padding:6px;background:#f5f5f5;"),
                  tags$td(format(cm$fp,big.mark=","),style="padding:6px;background:#ffebee;"),
                  tags$td(format(cm$tn,big.mark=","),style="padding:6px;background:#e8f5e9;font-weight:bold;"))),
        fluidRow(column(4,p(strong("Sens:"),sprintf("%.1f%%",cm$sens*100))),
                 column(4,p(strong("Spec:"),sprintf("%.1f%%",cm$spec*100))),
                 column(4,p(strong("Acc:"),sprintf("%.1f%%",cm$acc*100)))),
        fluidRow(column(4,p(strong("PPV:"),sprintf("%.1f%%",cm$ppv*100))),
                 column(4,p(strong("NPV:"),sprintf("%.1f%%",cm$npv*100))),
                 column(4,p(strong("N:"),format(cm$tp+cm$fp+cm$fn+cm$tn,big.mark=",")))))
    }
    tagList(fluidRow(
      column(4, cmp("Multilevel","#2c7fb8",ccm(test$outcome,test$pred,rv$roc_obj))),
      column(4, cmp("Logistic","#1a9850",if(!is.null(rv$glm_test_preds)) ccm(test$outcome,rv$glm_test_preds,rv$glm_roc))),
      column(4, cmp("XGBoost","#d73027",if(!is.null(rv$xgb_test_preds)) ccm(test$outcome,rv$xgb_test_preds,rv$xgb_roc)))))
  })

  # Capture overlay
  output$capture_overlay_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(NULL)
    n <- nrow(test); ne <- sum(test$outcome); if(ne==0) return(NULL)
    idx <- unique(c(1, round(seq(1,n,length.out=500)), n)); pct <- idx/n*100
    gc <- cumsum(test$outcome[order(test$pred,decreasing=TRUE)])/ne
    p <- plot_ly() |> add_trace(x=pct, y=gc[idx]*100, type="scatter", mode="lines", name="Multilevel",
                                 line=list(color="#2c7fb8",width=3), hovertemplate="Multilevel<br>Screen: %{x:.1f}%<br>Capture: %{y:.1f}%<extra></extra>")
    if(!is.null(rv$xgb_test_preds)) {
      xc <- cumsum(test$outcome[order(rv$xgb_test_preds,decreasing=TRUE)])/ne
      p <- p |> add_trace(x=pct, y=xc[idx]*100, type="scatter", mode="lines", name="XGBoost",
                           line=list(color="#d73027",width=3), hovertemplate="XGBoost<br>Screen: %{x:.1f}%<br>Capture: %{y:.1f}%<extra></extra>")
    }
    if(!is.null(rv$glm_test_preds)) {
      lc <- cumsum(test$outcome[order(rv$glm_test_preds,decreasing=TRUE)])/ne
      p <- p |> add_trace(x=pct, y=lc[idx]*100, type="scatter", mode="lines", name="Logistic",
                           line=list(color="#1a9850",width=3,dash="dot"), hovertemplate="Logistic<br>Screen: %{x:.1f}%<br>Capture: %{y:.1f}%<extra></extra>")
    }
    p |> add_trace(x=c(0,100),y=c(0,100),type="scatter",mode="lines",name="Random",line=list(color="gray60",width=2,dash="dash")) |>
      layout(title=list(text="Cumulative Capture",font=list(size=15)),
             xaxis=list(title="% Screened",ticksuffix="%",range=c(0,100)),
             yaxis=list(title="% Captured",ticksuffix="%",range=c(0,100)),
             legend=list(orientation="h",y=-0.15), margin=list(t=50,b=60))
  })

  # ROC with Youden
  output$roc_youden_plot <- renderPlotly({
    roc <- rv$roc_obj; if(is.null(roc)) return(NULL)
    fpr <- 1-roc$specificities; tpr <- roc$sensitivities; ys <- roc$youden$sens; yf <- 1-roc$youden$spec
    p <- plot_ly() |>
      add_trace(x=fpr,y=tpr,type="scatter",mode="lines",name=paste0("Multi (AUC=",round(roc$auc,3),")"),line=list(color="#2c7fb8",width=2.5))
    if(!is.null(rv$xgb_roc)) p <- p |> add_trace(x=1-rv$xgb_roc$specificities,y=rv$xgb_roc$sensitivities,type="scatter",mode="lines",
      name=paste0("XGB (",round(rv$xgb_roc$auc,3),")"),line=list(color="#d73027",width=2.5))
    if(!is.null(rv$glm_roc)) p <- p |> add_trace(x=1-rv$glm_roc$specificities,y=rv$glm_roc$sensitivities,type="scatter",mode="lines",
      name=paste0("GLM (",round(rv$glm_roc$auc,3),")"),line=list(color="#1a9850",width=2.5,dash="dot"))
    p |> add_trace(x=yf,y=ys,type="scatter",mode="markers",name=paste0("Youden (S=",round(ys,2),",Sp=",round(roc$youden$spec,2),")"),
                   marker=list(size=12,color="#ff7f00",symbol="star")) |>
      add_trace(x=c(0,1),y=c(0,1),type="scatter",mode="lines",line=list(color="gray60",width=1,dash="dash"),showlegend=FALSE,hoverinfo="skip") |>
      layout(title=list(text="ROC Curve",font=list(size=14)), xaxis=list(title="FPR",range=c(0,1)),
             yaxis=list(title="TPR",range=c(0,1)), legend=list(orientation="h",y=-0.2), margin=list(t=50))
  })

  # AUC by year
  output$auc_year_overlay_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(NULL)
    ty <- sort(unique(test$first_year))
    af <- function(preds) sapply(ty, function(yr) {
      i <- which(test$first_year==yr); if(sum(test$outcome[i])<3||sum(!test$outcome[i])<3) NA else compute_roc(test$outcome[i],preds[i])$auc
    })
    ym <- af(test$pred); if(sum(!is.na(ym))<2) return(NULL)
    p <- add_model_traces(plot_ly(), rv, ty, ym,
           if(!is.null(rv$xgb_test_preds)) af(rv$xgb_test_preds),
           if(!is.null(rv$glm_test_preds)) af(rv$glm_test_preds))
    p |> add_trace(x=range(ty),y=c(0.5,0.5),type="scatter",mode="lines",line=list(color="gray60",width=1,dash="dash"),showlegend=FALSE,hoverinfo="skip") |>
      layout(title=list(text="AUC by Year",font=list(size=14)), xaxis=list(title="Year"),
             yaxis=list(title="AUC",range=c(max(0.4,min(ym[!is.na(ym)])-0.05),1)),
             legend=list(orientation="h",y=-0.2), margin=list(t=40))
  })

  # Calibration over time
  output$cal_time_overlay_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(NULL)
    cs <- input$cal_models %||% c("multi","glm","xgb"); ty <- sort(unique(test$first_year))
    obs <- sapply(ty, \(yr) mean(test$outcome[test$first_year==yr]))
    pf <- function(preds) sapply(ty, \(yr) mean(preds[test$first_year==yr]))*100
    p <- plot_ly() |> add_trace(x=ty, y=obs*100, type="scatter", mode="lines+markers", name="Observed",
                                 line=list(color="#333",width=2.5), marker=list(color="#333",size=7))
    p <- add_model_traces(p, rv, ty,
      if("multi" %in% cs) pf(test$pred), if("xgb" %in% cs && !is.null(rv$xgb_test_preds)) pf(rv$xgb_test_preds),
      if("glm" %in% cs && !is.null(rv$glm_test_preds)) pf(rv$glm_test_preds))
    p |> layout(title=list(text="Calibration Over Time",font=list(size=14)), xaxis=list(title="Year"),
                yaxis=list(title="Elevated Rate (%)",ticksuffix="%"), legend=list(orientation="h",y=-0.2), margin=list(t=40))
  })

  # Calibration: Risk Decile
  output$cal_decile_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(NULL)
    cs <- input$cal_models %||% c("multi","glm","xgb"); test$decile <- ntile(test$pred,10)
    cal <- test |> group_by(decile) |> summarize(predicted=mean(pred)*100, observed=mean(outcome)*100, n=n(), .groups="drop")
    if(nrow(cal)==0) return(NULL); mv <- max(c(cal$predicted,cal$observed),na.rm=TRUE)*1.1
    p <- plot_ly(); ceq <- ""
    if("multi" %in% cs) {
      cl <- tryCatch(lm(observed~predicted, data=cal), error=function(e)NULL)
      if(!is.null(cl)) ceq <- paste0("y=",round(coef(cl)[1],2),"+",round(coef(cl)[2],2),"x (R\u00b2=",round(summary(cl)$r.squared,3),")")
      p <- p |> add_trace(data=cal, x=~predicted, y=~observed, type="scatter", mode="markers+text",
                marker=list(size=12,color="#2c7fb8"), text=~paste0("D",decile), textposition="top center",
                hovertext=~paste0("D",decile," Pred:",round(predicted,1),"% Obs:",round(observed,1),"%"), hoverinfo="text", name="Multilevel")
    }
    if("glm" %in% cs && !is.null(rv$glm_test_preds)) {
      test$gd <- ntile(rv$glm_test_preds,10)
      cg <- test |> group_by(decile=gd) |> summarize(predicted=mean(rv$glm_test_preds[cur_group_rows()])*100, observed=mean(outcome)*100, n=n(), .groups="drop")
      mv <- max(mv, max(c(cg$predicted,cg$observed),na.rm=TRUE)*1.1)
      p <- p |> add_trace(data=cg, x=~predicted, y=~observed, type="scatter", mode="markers",
                marker=list(size=10,color="#1a9850",symbol="triangle-up"), name="Logistic")
    }
    if("xgb" %in% cs && !is.null(rv$xgb_test_preds)) {
      test$xd <- ntile(rv$xgb_test_preds,10)
      cx <- test |> group_by(decile=xd) |> summarize(predicted=mean(rv$xgb_test_preds[cur_group_rows()])*100, observed=mean(outcome)*100, n=n(), .groups="drop")
      mv <- max(mv, max(c(cx$predicted,cx$observed),na.rm=TRUE)*1.1)
      p <- p |> add_trace(data=cx, x=~predicted, y=~observed, type="scatter", mode="markers",
                marker=list(size=10,color="#d73027",symbol="square"), name="XGBoost")
    }
    p |> add_trace(x=c(0,mv),y=c(0,mv),type="scatter",mode="lines",line=list(color="gray60",width=1,dash="dash"),showlegend=FALSE,hoverinfo="skip") |>
      layout(title=list(text="Calibration by Decile",font=list(size=14)),
             xaxis=list(title="Predicted (%)",ticksuffix="%"), yaxis=list(title="Observed (%)",ticksuffix="%"),
             legend=list(orientation="h",y=-0.15), margin=list(t=50),
             annotations=if(nchar(ceq)>0)list(list(x=0.02,y=0.98,xref="paper",yref="paper",text=ceq,showarrow=FALSE,font=list(size=10,color="#666")))else list())
  })

  # Calibration: 10% bins
  output$cal_pctbin_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(NULL)
    d10 <- quantile(test$pred, 0.9, na.rm=TRUE)
    test$pct_bin <- cut(test$pred, breaks=seq(0,1,0.1), include.lowest=TRUE,
                        labels=paste0(seq(0,90,10),"-",seq(10,100,10),"%"))
    cal <- test |> group_by(pct_bin) |>
      summarize(predicted=mean(pred)*100, observed=mean(outcome)*100, n=n(),
                min_pred=min(pred), .groups="drop") |>
      filter(n>=1) |> mutate(in_d10=min_pred>=d10)
    if(nrow(cal)==0) return(NULL); mv <- max(c(cal$predicted,cal$observed),na.rm=TRUE)*1.1
    p <- plot_ly() |>
      add_trace(data=cal|>filter(!in_d10), x=~predicted, y=~observed, type="scatter", mode="markers",
                marker=list(size=12,color="#2c7fb8"), hovertext=~paste0(pct_bin,"\nPred:",round(predicted,1),"%\nObs:",round(observed,1),"%"), hoverinfo="text", name="Standard") |>
      add_trace(data=cal|>filter(in_d10), x=~predicted, y=~observed, type="scatter", mode="markers",
                marker=list(size=14,color="#d73027",symbol="diamond"), hovertext=~paste0(pct_bin," [D10]\nPred:",round(predicted,1),"%\nObs:",round(observed,1),"%"), hoverinfo="text", name="In D10") |>
      add_trace(x=c(0,mv),y=c(0,mv),type="scatter",mode="lines",line=list(color="gray60",width=1,dash="dash"),showlegend=FALSE,hoverinfo="skip")
    cl <- tryCatch(lm(observed~predicted,data=cal), error=function(e)NULL)
    ceq <- if(!is.null(cl)) paste0("y=",round(coef(cl)[1],2),"+",round(coef(cl)[2],2),"x (R\u00b2=",round(summary(cl)$r.squared,3),")") else ""
    p |> layout(title=list(text="Calibration by 10% Bins",font=list(size=14)),
                xaxis=list(title="Predicted (%)",ticksuffix="%"), yaxis=list(title="Observed (%)",ticksuffix="%"),
                legend=list(orientation="h",y=-0.15), margin=list(t=50),
                annotations=if(nchar(ceq)>0)list(list(x=0.02,y=0.98,xref="paper",yref="paper",text=ceq,showarrow=FALSE,font=list(size=10,color="#666")))else list())
  })

  # Mean BLL by screening quintile
  output$mean_bll_screened_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0||!"first_bll" %in% names(test)) return(NULL)
    n <- nrow(test); om <- mean(test$first_bll, na.rm=TRUE)
    bb <- function(ord, mn) {
      bo <- test$first_bll[ord]; bi <- pmin(ceiling(seq_len(n)/n*5), 5)
      data.frame(bin=1:5, label=paste0("Q",1:5," (",(0:4)*20,"-",(1:5)*20,"%)"),
                 mean_bll=tapply(bo,bi,mean,na.rm=TRUE), model=mn, stringsAsFactors=FALSE)
    }
    pd <- bb(order(test$pred,decreasing=TRUE), "Multilevel")
    if(!is.null(rv$xgb_test_preds)) pd <- rbind(pd, bb(order(rv$xgb_test_preds,decreasing=TRUE), "XGBoost"))
    if(!is.null(rv$glm_test_preds)) pd <- rbind(pd, bb(order(rv$glm_test_preds,decreasing=TRUE), "Logistic"))
    mc <- c(Multilevel="#2c7fb8",XGBoost="#d73027",Logistic="#1a9850"); md <- c(Multilevel="solid",XGBoost="solid",Logistic="dot")
    p <- plot_ly()
    for(m in unique(pd$model)) {
      d <- pd[pd$model==m,]
      p <- p |> add_trace(x=d$label,y=d$mean_bll,name=m,type="scatter",mode="lines+markers",
                           line=list(color=mc[m],width=2.5,dash=md[m]),marker=list(color=mc[m],size=8),
                           hovertext=paste0(m," | ",d$label,"\nBLL: ",round(d$mean_bll,2)), hoverinfo="text")
    }
    p |> add_trace(x=c(pd$label[1],pd$label[5]),y=c(om,om),type="scatter",mode="lines",name="Overall",
                   line=list(color="gray60",width=2,dash="dash")) |>
      layout(title=list(text="Mean BLL by Screening Quintile",font=list(size=13)),
             xaxis=list(title="Bin (highest risk first)",categoryorder="array",
                        categoryarray=paste0("Q",1:5," (",(0:4)*20,"-",(1:5)*20,"%)")),
             yaxis=list(title="Mean BLL (\u00b5g/dL)"), legend=list(orientation="h",y=-0.2), margin=list(t=50))
  })

  # Capture by year
  output$capture_by_year_plot <- renderPlotly({
    test <- rv$test_preds; if(is.null(test)||nrow(test)==0) return(NULL)
    ty <- sort(unique(test$first_year))
    cf <- function(preds) sapply(ty, function(yr) {
      i <- which(test$first_year==yr); ne <- sum(test$outcome[i]==1); if(ne<2) return(NA)
      nh <- ceiling(length(i)/2); o <- order(preds[i],decreasing=TRUE); 100*sum(test$outcome[i[o[1:nh]]])/ne
    })
    ym <- cf(test$pred); if(sum(!is.na(ym))<2) return(NULL)
    p <- add_model_traces(plot_ly(), rv, ty, ym,
           if(!is.null(rv$xgb_test_preds)) cf(rv$xgb_test_preds),
           if(!is.null(rv$glm_test_preds)) cf(rv$glm_test_preds))
    p |> add_trace(x=range(ty),y=c(50,50),type="scatter",mode="lines",line=list(color="gray60",width=1,dash="dash"),showlegend=FALSE,hoverinfo="skip") |>
      layout(title=list(text="50% Screen Capture by Year",font=list(size=13)),
             xaxis=list(title="Year"), yaxis=list(title="% Captured",ticksuffix="%"),
             legend=list(orientation="h",y=-0.2), margin=list(t=50))
  })

  # Sample tab plots
  get_model_data_for_plot <- reactive({
    acs_mode <- input$acs_mode %||% "nof"; use_max <- isTRUE(input$test_selection=="max")
    use_all_addr <- isTRUE(input$model_addr_pop=="all")
    geo_key <- if(input$random_effect=="bg_geoid") "bg" else "tract"
    resolve_model_data(acs_mode, geo_key, use_max, use_all_addr)
  })

  output$sample_crosstab_table <- renderDT({
    md <- bind_rows(rv$train_sample, rv$test_sample)
    if(is.null(md)||nrow(md)==0||!"outcome" %in% names(md)) return(NULL)
    te <- if(!input$train_all) input$train_years[2] else max(md$first_year,na.rm=TRUE)
    md$Split <- if_else(md$first_year<=te,"Train","Test")
    ct <- md |> group_by(Split) |>
      summarize(Addresses=n(), `Outcome=1`=sum(outcome==1,na.rm=TRUE), `Outcome=0`=sum(outcome==0,na.rm=TRUE),
                `% Elevated`=round(100*mean(outcome,na.rm=TRUE),2), `Mean BLL`=round(mean(first_bll,na.rm=TRUE),2),
                `Median BLL`=round(median(first_bll,na.rm=TRUE),2), .groups="drop") |> arrange(desc(Split))
    datatable(ct, options=list(dom="t",pageLength=5), rownames=FALSE)
  })

  output$sample_plot_addr <- renderPlotly({
    if(is.null(addr_by_year)||nrow(addr_by_year)==0) return(NULL)
    p <- ggplot(addr_by_year, aes(sample_year, n_addr, fill=type,
                text=paste0("Year: ",sample_year,"\n",type,": ",format(n_addr,big.mark=",")))) +
      geom_col(position="stack") + scale_fill_manual(values=c("2+ Children"="#2c7fb8","1 Child"="#b0d5e8")) +
      scale_y_continuous(labels=comma) + labs(title="Addresses by Year",x="Year",y="Addresses",fill="") +
      theme_minimal(base_size=11) + theme(legend.position="bottom")
    ggplotly(p, tooltip="text") |> layout(autosize=TRUE)
  })

  output$sample_plot_tested <- renderPlotly({
    md <- get_model_data_for_plot(); if(is.null(md)) return(NULL)
    if(!"outcome" %in% names(md)) { md <- bind_rows(rv$train_sample,rv$test_sample); if(is.null(md)||!"outcome" %in% names(md)) return(NULL) }
    d <- md |> group_by(sample_year=first_year) |>
      summarize(Elevated=sum(outcome==1,na.rm=TRUE), `Not Elevated`=sum(outcome==0,na.rm=TRUE), .groups="drop")
    te <- if(!input$train_all) input$train_years[2] else max(d$sample_year); d$split <- if_else(d$sample_year<=te,"Train","Test")
    plot_ly(d, x=~sample_year) |>
      add_bars(y=~Elevated, name="Elevated", marker=list(color="#d73027"), hoverinfo="text",
               hovertext=~paste0(sample_year,"\nElevated: ",format(Elevated,big.mark=","),"\n",split), textposition="none") |>
      add_bars(y=~`Not Elevated`, name="Not Elevated", marker=list(color="#b0d5e8"), hoverinfo="text",
               hovertext=~paste0(sample_year,"\nNot Elev: ",format(`Not Elevated`,big.mark=","),"\n",split), textposition="none") |>
      layout(barmode="stack", title=list(text="Addresses: Elevated vs Not",font=list(size=13)),
             xaxis=list(title="Year"), yaxis=list(title="Count"), legend=list(orientation="h",y=-0.2), autosize=TRUE)
  })

  output$sample_plot_bll <- renderPlotly({
    md <- bind_rows(rv$train_sample, rv$test_sample)
    if(is.null(md)||nrow(md)==0||!"first_bll" %in% names(md)) { md <- get_model_data_for_plot(); if(is.null(md)) return(NULL) }
    d <- md |> group_by(sample_year=first_year) |> summarize(bll=mean(first_bll,na.rm=TRUE), .groups="drop")
    te <- if(!input$train_all) input$train_years[2] else max(d$sample_year)
    d$split <- if_else(d$sample_year<=te,"Train","Test"); use_log <- isTRUE(input$sample_bll_log)
    p <- ggplot(d, aes(sample_year,bll,group=1, text=paste0("Year: ",sample_year,"\nBLL: ",round(bll,2),"\n",split))) +
      geom_line(linewidth=0.5,color="#555") + geom_point(aes(color=split),size=2) +
      scale_color_manual(values=c(Train="#2c7fb8",Test="#d73027")) +
      labs(title="Mean First BLL by Year",x="Year",y="\u00b5g/dL",color="") + theme_minimal(base_size=11) + theme(legend.position="bottom")
    if(use_log) p <- p + scale_y_log10()
    ggplotly(p, tooltip="text") |> layout(autosize=TRUE)
  })

  output$sample_plot_elev <- renderPlotly({
    md <- get_model_data_for_plot(); if(is.null(md)) return(NULL)
    if(!"outcome" %in% names(md)) { md <- bind_rows(rv$train_sample,rv$test_sample); if(is.null(md)||!"outcome" %in% names(md)) return(NULL) }
    d <- md |> group_by(sample_year=first_year) |> summarize(pct=100*mean(outcome,na.rm=TRUE), .groups="drop")
    te <- if(!input$train_all) input$train_years[2] else max(d$sample_year); d$split <- if_else(d$sample_year<=te,"Train","Test")
    p <- ggplot(d, aes(sample_year,pct,group=1, text=paste0("Year: ",sample_year,"\n% Elev: ",round(pct,1),"%\n",split))) +
      geom_line(linewidth=0.5,color="#555") + geom_point(aes(color=split),size=2) +
      scale_color_manual(values=c(Train="#2c7fb8",Test="#d73027")) +
      labs(title="% Elevated (Outcome)",x="Year",y="%",color="") + theme_minimal(base_size=11) + theme(legend.position="bottom")
    ggplotly(p, tooltip="text") |> layout(autosize=TRUE)
  })

  output$dl_sample_table <- downloadHandler(
    filename=function() paste0("sample_summary_",Sys.Date(),".csv"),
    content=function(file) {
      md <- get_model_data_for_plot(); if(is.null(md)) { write.csv(data.frame(),file); return() }
      d <- md |> group_by(Year=first_year) |>
        summarize(Addresses=n(), `% Elevated`=round(100*mean(outcome,na.rm=TRUE),2),
                  `Mean BLL`=round(mean(first_bll,na.rm=TRUE),2), .groups="drop")
      write.csv(d, file, row.names=FALSE)
    })

  # ===== PREDICTION MAP =====

  observeEvent(c(input$pred_type, input$pred_geo, input$pred_class, input$pred_palette,
                 input$pred_model, rv$model_run, input$map_opacity), {
    pred_model <- input$pred_model %||% "glmer"
    preds <- if(input$pred_type=="test") rv$test_preds else rv$full_preds
    if(is.null(preds)||nrow(preds)==0) return()
    alt <- switch(pred_model, glm=rv$glm_test_preds, xgb=rv$xgb_test_preds, NULL)
    if(!is.null(alt) && input$pred_type=="test" && length(alt)==nrow(preds)) preds$pred <- alt
    else if(pred_model!="glmer" && !is.null(alt)) return()
    PC <- c(county="county",lhd="lhd",tract="tract_geoid",bg="bg_geoid")
    PS <- list(county=ne_counties,lhd=lhd_bounds,tract=ne_tracts,bg=ne_bgs)
    PI <- c(county="NAME",lhd="lhd",tract="GEOID",bg="GEOID")
    gc <- PC[input$pred_geo]; gs <- PS[[input$pred_geo]]; gi <- PI[input$pred_geo]
    if(!gc %in% names(preds)||is.null(gs)) return()
    pc <- if("pred" %in% names(preds)) "pred" else if("pred_prob" %in% names(preds)) "pred_prob" else return()
    agg <- preds |> filter(!is.na(.data[[gc]])) |>
      group_by(geo_id=.data[[gc]]) |>
      summarize(mean_pred=mean(.data[[pc]],na.rm=TRUE), n_addr=n(), pct_elevated=100*mean(outcome,na.rm=TRUE), .groups="drop")
    geo <- gs |> left_join(agg, by=setNames("geo_id",gi))
    vals <- geo$mean_pred[!is.na(geo$mean_pred)]; if(length(vals)==0) return()
    cm <- input$pred_class %||% "quintile"
    if(cm %in% c("quintile","decile")) {
      nb <- if(cm=="quintile") 5 else 10; qb <- quantile(vals, probs=seq(0,1,1/nb), na.rm=TRUE)
      rc <- get_risk_palette(input$pred_palette, nb)
      geo$bin <- pmin(pmax(findInterval(geo$mean_pred, qb, rightmost.closed=TRUE),1), nb)
      geo$fill <- ifelse(is.na(geo$mean_pred), "#e0e0e0", rc[geo$bin])
      pfx <- if(cm=="quintile") "Q" else "D"
      ll <- paste0(pfx,1:nb,": ",sprintf("%.1f%%",qb[1:nb]*100)," \u2013 ",sprintf("%.1f%%",qb[2:(nb+1)]*100))
    } else {
      mx <- max(0.20, max(vals,na.rm=TRUE)); rp <- colorNumeric(get_risk_palette(input$pred_palette,100), domain=c(0,mx), na.color="#e0e0e0")
      geo$fill <- ifelse(is.na(geo$mean_pred), "#e0e0e0", sapply(geo$mean_pred, rp))
    }
    gl <- geo[[gi]]
    geo$risk_cat <- case_when(is.na(geo$mean_pred)~"No Data", geo$mean_pred>=0.15~"\u26a0 High",
                              geo$mean_pred>=0.08~"Moderate", geo$mean_pred>=0.04~"Low-Moderate", TRUE~"Lower")
    geo$tip <- paste0("<div style='font-size:13px;line-height:1.6;'><b>",gl,"</b><br>",
      "<span style='font-size:18px;font-weight:bold;color:",
      ifelse(is.na(geo$mean_pred),"#999",ifelse(geo$mean_pred>=0.15,"#bd0026",ifelse(geo$mean_pred>=0.08,"#fc4e2a","#333"))),";'>",
      ifelse(!is.na(geo$mean_pred), sprintf("%.1f%%",geo$mean_pred*100), "N/A"), "</span> risk<br>",
      "<b>",geo$risk_cat,"</b><br>Observed: ",ifelse(!is.na(geo$pct_elevated),sprintf("%.1f%%",geo$pct_elevated),"N/A"),"<br>",
      "Addresses: ",ifelse(!is.na(geo$n_addr),format(geo$n_addr,big.mark=","),"0"),"</div>")
    proxy <- leafletProxy("pred_map") |> clearShapes() |> clearControls() |>
      addPolygons(data=geo, fillColor=~fill, fillOpacity=input$map_opacity/100, color="#444", weight=0.3, opacity=0.5,
                  highlightOptions=highlightOptions(weight=2.5,color="#000",fillOpacity=min(1,input$map_opacity/100+0.1)),
                  label=lapply(geo$tip, HTML))
    if(cm %in% c("quintile","decile"))
      proxy |> addLegend("bottomright",colors=c(rc,"#e0e0e0"),labels=c(ll,"No Data"),title=HTML("<b>Predicted Risk</b>"),opacity=0.9)
    else
      proxy |> addLegend("bottomright",pal=rp,values=seq(0,mx,length.out=6),title=HTML("<b>Predicted Risk</b>"),
                         labFormat=labelFormat(suffix="%",transform=\(x)round(x*100,1)),opacity=0.9)
  })

  output$pred_map <- renderLeaflet({
    leaflet() |>
      addProviderTiles(providers$CartoDB.Positron, group="Light") |>
      addProviderTiles(providers$CartoDB.DarkMatter, group="Dark") |>
      addProviderTiles(providers$Esri.WorldImagery, group="Satellite") |>
      addLayersControl(baseGroups=c("Light","Dark","Satellite"), position="topright") |>
      setView(-99.5, 41.5, zoom=7)
  })
}

shinyApp(ui, server)
