# Steps
# 1. CTRL + A
# 2. CTRL + enter (or click run)
# 3. Use app. Wait ~60 seconds upon start for default risk model to run
# If maps aren't loading or appear stuck, click on a different setting to refresh

# Run if these packages not already installed for you (very likely)
packages <- c(
  "shiny", "shinydashboard", "leaflet", "sf", "dplyr",
  "ggplot2", "DT", "plotly", "lme4", "broom.mixed", "pROC"
)

install.packages(setdiff(packages, rownames(installed.packages())))
invisible(lapply(packages, library, character.only = TRUE))


library(shiny)
library(shinydashboard)
library(leaflet)
library(sf)
library(dplyr)
library(ggplot2)
library(DT)
library(plotly)
library(lme4)
library(broom.mixed)
library(pROC)

DATA_PATH  <- "K:/CLPPP/TrentonS/shiny_data"
BLL_CUTOFF <- 3.5 # 2021 CDC change from 5 (started in 2011 from 10)

# --- Helpers ---

# Aggregate lead testing data by a geographic column.
# Returns child counts, address counts, mean BLL, pct elevated, suppression flag.
summarize_lead <- function(data, geo_col, sup = 0) {
  data <- data %>% filter(!is.na(.data[[geo_col]]))
  
  child <- data %>%
    group_by(geo_id = .data[[geo_col]]) %>%
    summarize(n_children = n(), n_elevated = sum(elevated, na.rm = TRUE),
              mean_bll = mean(result, na.rm = TRUE), county = first(county), .groups = "drop")
  
  addr <- data %>%
    group_by(geo_id = .data[[geo_col]], address_id) %>%
    summarize(e = max(elevated, na.rm = TRUE), .groups = "drop") %>%
    group_by(geo_id) %>%
    summarize(n_addresses = n(), addr_elevated = sum(e > 0, na.rm = TRUE), .groups = "drop")
  
  left_join(child, addr, by = "geo_id") %>%
    mutate(pct_elevated = 100 * n_elevated / n_children, suppressed = n_children < sup)
}

# State-level summary (children or addresses, by year or overall)
state_summary <- function(data, by_year = FALSE, by_address = FALSE) {
  if (by_address) {
    grp <- if (by_year) data %>% group_by(sample_year, address_id) else data %>% group_by(address_id)
    addr <- grp %>%
      summarize(e = max(elevated, na.rm = TRUE), bll = mean(result, na.rm = TRUE), .groups = "drop")
    if (by_year) addr <- addr %>% group_by(sample_year)
    addr %>%
      summarize(mean_bll = mean(bll, na.rm = TRUE), n_children = n(),
                n_elevated = sum(e > 0, na.rm = TRUE),
                pct_elevated = 100 * mean(e > 0, na.rm = TRUE), .groups = "drop")
  } else {
    if (by_year) data <- data %>% group_by(sample_year)
    data %>%
      summarize(mean_bll = mean(result, na.rm = TRUE), n_children = n(),
                n_elevated = sum(elevated, na.rm = TRUE),
                pct_elevated = 100 * mean(elevated, na.rm = TRUE), .groups = "drop")
  }
}

# Format a summarize_lead result for display in a DataTable
format_table <- function(df, id_name, unit = "Children", show_county = FALSE) {
  tested   <- if (unit == "Addresses") df$n_addresses  else df$n_children
  elevated <- if (unit == "Addresses") df$addr_elevated else df$n_elevated
  out <- tibble(!!id_name := df$geo_id, Tested = tested, Elevated = elevated,
                `% Elevated` = round(100 * elevated / pmax(tested, 1), 2),
                `Average BLL` = round(df$mean_bll, 2), Suppressed = df$suppressed)
  if (show_county) out <- bind_cols(tibble(County = df$county), out)
  out
}

# --- Load Data ---

load_rds <- function(name) {
  path <- file.path(DATA_PATH, paste0(name, ".rds"))
  if (file.exists(path)) readRDS(path) else NULL
}

first_test    <- load_rds("first_test")
model_data    <- load_rds("model_data")
address_stats <- load_rds("address_stats")
county_stats  <- load_rds("county_stats")
lhd_stats     <- load_rds("lhd_stats")
bg_stats      <- load_rds("bg_stats")
tract_stats   <- load_rds("tract_stats")
ne_counties   <- load_rds("ne_counties")
lhd_bounds    <- load_rds("lhd_bounds")
ne_bgs        <- load_rds("ne_bgs")
ne_tracts     <- load_rds("ne_tracts")
bg_census     <- load_rds("bg_census")
tract_census  <- load_rds("tract_census")
lhd_map       <- load_rds("lhd_map")

# Defensive fix: str_to_title("MCPHERSON") -> "Mcpherson" in data_prep
fix_mc <- function(df) {
  if (!is.null(df) && "county" %in% names(df))
    df %>% mutate(county = if_else(county == "Mcpherson", "McPherson", county))
  else df
}
first_test   <- fix_mc(first_test)
county_stats <- fix_mc(county_stats)
lhd_map      <- fix_mc(lhd_map)

years     <- if (!is.null(first_test)) sort(unique(first_test$sample_year)) else 2010:2023
lhd_names <- if (!is.null(lhd_stats)) sort(unique(lhd_stats$lhd)) else character(0)

# Pre-compute state-level summaries (4 combinations: year/all × children/addresses)
state_avgs         <- state_summary(first_test, by_year = TRUE)
state_avg_all      <- state_summary(first_test)
state_avgs_addr    <- state_summary(first_test, by_year = TRUE, by_address = TRUE)
state_avg_all_addr <- state_summary(first_test, by_address = TRUE)

VAR_LABELS <- c(
  pct_pre1980 = "Pre-1980 Housing (%)", pct_poverty = "Poverty Rate (%)",
  pct_renter = "Renter-Occupied (%)", median_income = "Median Income ($)",
  pct_nonwhite = "Non-White (%)", log_housing_density = "Housing Density (log)",
  housing_units = "Housing Units", pct_children_u6 = "Children Under 6 (%)",
  first_bll = "First Child's BLL", first_year = "Year of First Test"
)

METRIC_LABELS <- c(
  n_elevated = "Elevated Count", n_children = "Total Tested",
  pct_elevated = "Percent Elevated (%)", mean_bll = "Average BLL (\u00b5g/dL)",
  n_addresses = "Addresses Tested", addr_elevated = "Addresses Elevated"
)

get_label <- function(v) {
  if (v %in% names(VAR_LABELS)) unname(VAR_LABELS[v])
  else if (v %in% names(METRIC_LABELS)) unname(METRIC_LABELS[v])
  else gsub("_", " ", tools::toTitleCase(v))
}

get_bll_color <- function(bll) {
  ifelse(bll >= 20, "#67000d",
         ifelse(bll >= 10, "#a50f15",
                ifelse(bll >= 5,  "#d73027",
                       ifelse(bll >= 3.5, "#fc8d59", "#c9a84c"))))
}

# ==============================================================================
# UI
# ==============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "NE Lead Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("home")),
      menuItem("Surveillance Map", tabName = "map", icon = icon("map")),
      menuItem("Priority Addresses", tabName = "ranked", icon = icon("list-ol")),
      menuItem("Census Data", tabName = "census", icon = icon("chart-bar")),
      menuItem("Data Tables", tabName = "tables", icon = icon("table")),
      menuItem("Risk Model", tabName = "model", icon = icon("chart-line"))
    ),
    hr(),
    numericInput("suppress_n", "Suppress areas with fewer than:", 10, min = 0, max = 50, step = 5),
    p(em("Applies to maps and tables"), style = "padding: 0 15px; font-size: 11px;")
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML("
      .section-title { font-size: 18px; font-weight: 600; margin-bottom: 5px; }
      .section-desc { color: #444; margin-bottom: 15px; font-size: 13px; line-height: 1.5; }
      .metric-box { background: #f8f9fa; border-left: 3px solid #2c7fb8;
                    padding: 12px; margin: 8px 0; border-radius: 4px; }
      .state-avg-box { background: #e8f4f8; border: 1px solid #2c7fb8; padding: 8px;
                       border-radius: 4px; margin-top: 10px; font-size: 12px; }
      .compare-box { background: #f0f7f0; border: 1px solid #4a9; padding: 15px;
                     border-radius: 6px; margin: 10px 0; }
      .warning-box { background: #fff3cd; border: 1px solid #ffc107; padding: 10px;
                     border-radius: 4px; margin: 5px 0; font-size: 12px; }
      .cor-high { background-color: #f8d7da !important; }
      .cor-borderline { background-color: #fff3cd !important; }
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
                strong("Surveillance Map:"), " View testing metrics by Local Health Department, county, census tract,
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
                         "Children Elevated", icon = icon("exclamation-triangle"), color = "blue", width = 6)
              ),
              fluidRow(
                valueBox(format(n_distinct(first_test$address_id), big.mark = ","), "Addresses Tested",
                         icon = icon("home"), color = "green", width = 6),
                valueBox({
                  addr_any <- first_test %>% group_by(address_id) %>%
                    summarize(e = max(elevated, na.rm = TRUE), .groups = "drop")
                  paste0(format(sum(addr_any$e > 0, na.rm = TRUE), big.mark = ","),
                         " (", round(100 * mean(addr_any$e > 0, na.rm = TRUE), 2), "%)")
                }, "Addresses with Elevated Child", icon = icon("home"), color = "olive", width = 6)
              ),
              fluidRow(
                valueBox(format(nrow(address_stats), big.mark = ","), "Addresses with 2+ Children Tested",
                         icon = icon("users"), color = "yellow", width = 6),
                valueBox({
                  paste0(format(sum(address_stats$n_elevated > 0, na.rm = TRUE), big.mark = ","), " (",
                         round(100 * mean(address_stats$n_elevated > 0, na.rm = TRUE), 2), "%)")
                }, "2+ Child Addresses with Elevated", icon = icon("users"), color = "orange", width = 6)
              )
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
                    leafletOutput("surv_map", height = 400),
                    hr(),
                    plotlyOutput("state_avg_plot", height = 200),
                    hr(),
                    h5("Annual Summary Statistics"),
                    DTOutput("trend_summary_table")
                ),
                box(width = 3, title = "Map Options",
                    selectInput("geo_level", "Geographic Level:",
                                c("Health Department" = "LHD", "County" = "County",
                                  "Census Tract" = "Tract", "Block Group" = "Block Group")),
                    selectInput("fill_var", "Measure:",
                                c("Elevated Count" = "n_elevated",
                                  "Total Tested" = "n_children",
                                  "Percent Elevated (%)" = "pct_elevated",
                                  "Average BLL (\u00b5g/dL)" = "mean_bll")),
                    checkboxInput("use_quintile", "Display by Quintile", FALSE),
                    radioButtons("unit_type", "Count By:", c("Children", "Addresses"), inline = TRUE),
                    selectInput("year_filter", "Year:", c("All Years" = "All", years)),
                    hr(),
                    div(class = "state-avg-box",
                        h5("Nebraska", style = "margin: 0 0 5px 0;"),
                        uiOutput("state_avg_text")
                    ),
                    hr(),
                    h5("Selected Area:"),
                    verbatimTextOutput("click_info")
                )
              )
      ),
      
      # ===== PRIORITY ADDRESSES =====
      tabItem(tabName = "ranked",
              h3(class = "section-title", "Priority Addresses for Outreach"),
              p(class = "section-desc",
                "Addresses ranked by elevated count or average BLL. Use this to prioritize outreach, inspection referrals,
           or remediation assistance. Colors indicate BLL severity. Numbers show rank within the selected health department."),
              fluidRow(
                box(width = 9,
                    leafletOutput("ranked_map", height = 380),
                    hr(),
                    downloadButton("dl_ranked", "Download Displayed Addresses"),
                    br(), br(),
                    div(style = "height: 250px; overflow-y: auto;", DTOutput("ranked_table"))
                ),
                box(width = 3, title = "Options",
                    selectInput("ranked_lhd", "Health Department:", c("Statewide" = "All", lhd_names)),
                    sliderInput("top_n", "Show Top N per LHD:", 5, 50, 20, 5),
                    radioButtons("rank_by", "Rank By:",
                                 c("Elevated Count" = "n_elevated", "Average BLL" = "mean_bll")),
                    hr(),
                    radioButtons("map_tiles", "Base Map:", c("Street" = "street", "Satellite" = "sat"), inline = TRUE),
                    hr(),
                    h5("BLL Color Key:"),
                    HTML('<span style="color:#67000d">\u25cf</span> \u226520 \u00b5g/dL<br>
                  <span style="color:#a50f15">\u25cf</span> 10-20<br>
                  <span style="color:#d73027">\u25cf</span> 5-10<br>
                  <span style="color:#fc8d59">\u25cf</span> 3.5-5<br>
                  <span style="color:#c9a84c">\u25cf</span> <3.5')
                )
              )
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
                    hr(),
                    h5("Lead Outcomes by Quintile"),
                    downloadButton("dl_quintile", "Download Quintile Summary"),
                    DTOutput("census_quintile_table"),
                    hr(),
                    h5("Area Details"),
                    downloadButton("dl_area_detail", "Download Area Details"),
                    DTOutput("census_area_table")
                ),
                box(width = 3, title = "Options",
                    selectInput("census_var", "Map Variable:",
                                c("Pre-1980 Housing (%)" = "pct_pre1980",
                                  "Poverty Rate (%)" = "pct_poverty",
                                  "Renter-Occupied (%)" = "pct_renter",
                                  "Median Income ($)" = "median_income",
                                  "Non-White (%)" = "pct_nonwhite",
                                  "Housing Density (log)" = "log_housing_density",
                                  "Housing Units" = "housing_units",
                                  "Children Under 6 (%)" = "pct_children_u6")),
                    radioButtons("census_geo", "Level:", c("Tract" = "tract", "Block Group" = "bg"), inline = TRUE),
                    checkboxInput("census_quintile", "Display by Quintile", FALSE),
                    hr(),
                    h5("Selected Area:"),
                    verbatimTextOutput("census_info")
                )
              )
      ),
      
      # ===== DATA TABLES =====
      tabItem(tabName = "tables",
              h3(class = "section-title", "Data Tables"),
              p(class = "section-desc", "Browse and download datasets. Choose to download displayed rows only or all rows."),
              fluidRow(
                column(3, sliderInput("table_years", "Years:", min(years), max(years),
                                      c(min(years), max(years)), sep = "")),
                column(3, selectInput("table_rows", "Display Rows:", c(10, 25, 50, 100), selected = 25)),
                column(3, radioButtons("table_unit", "Unit:", c("Children", "Addresses"), inline = TRUE)),
                column(3, radioButtons("table_geo", "Geography:", c("Tract" = "tract", "Block Group" = "bg"), inline = TRUE))
              ),
              tabsetPanel(
                tabPanel("Health Departments", br(),
                         fluidRow(
                           column(6, downloadButton("dl_lhd_disp", "Download Displayed")),
                           column(6, downloadButton("dl_lhd_all", "Download All"))
                         ), br(), DTOutput("lhd_table")),
                tabPanel("Counties", br(),
                         fluidRow(
                           column(6, downloadButton("dl_county_disp", "Download Displayed")),
                           column(6, downloadButton("dl_county_all", "Download All"))
                         ), br(), DTOutput("county_table")),
                tabPanel("Small Area", br(),
                         fluidRow(
                           column(6, downloadButton("dl_geo_disp", "Download Displayed")),
                           column(6, downloadButton("dl_geo_all", "Download All"))
                         ), br(), DTOutput("geo_table")),
                tabPanel("Addresses", br(),
                         fluidRow(
                           column(6, downloadButton("dl_addr_disp", "Download Displayed")),
                           column(6, downloadButton("dl_addr_all", "Download All"))
                         ), br(), DTOutput("addr_table"))
              )
      ),
      
      # ===== RISK MODEL =====
      tabItem(tabName = "model",
              h3(class = "section-title", "Predictive Risk Model"),
              p(class = "section-desc",
                "This model predicts which addresses are most likely to have future children with elevated BLL.
           It uses mixed-effects logistic regression with geographic random effects to account for neighborhood clustering."),
              fluidRow(
                box(width = 4, title = "Model Settings",
                    sliderInput("train_years", "Training Years:", min(years), max(years), c(2010, 2016), sep = ""),
                    checkboxInput("train_all", "Train on all years", FALSE),
                    div(class = "warning-box",
                        em("Train on all years: Use when you want predictions for future addresses but cannot evaluate
                    model performance (no held-out test set). Useful for deployment, not validation.")),
                    hr(),
                    radioButtons("random_effect", "Random Intercept Level:",
                                 c("Census Tract" = "tract_geoid", "Block Group" = "bg_geoid"),
                                 selected = "tract_geoid", inline = TRUE),
                    p(em("Adds a neighborhood-specific adjustment. Tracts are larger (more stable estimates);
                  block groups are smaller (more local but noisier)."), style = "font-size: 11px; color: #666;"),
                    hr(),
                    checkboxInput("scale_vars", "Standardize variables", TRUE),
                    div(class = "warning-box",
                        em("Standardization: Makes odds ratios comparable across variables with different scales.
                    However, ORs then represent a 1 standard deviation change, not a 1-unit change.
                    Also helps model convergence.")),
                    hr(),
                    h5("Neighborhood Variables:"),
                    checkboxGroupInput("census_covars", NULL,
                                       choices = setNames(names(VAR_LABELS)[1:8], VAR_LABELS[1:8]),
                                       selected = "pct_pre1980"),
                    h5("Address Variables:"),
                    checkboxGroupInput("addr_covars", NULL,
                                       choices = setNames(names(VAR_LABELS)[9:10], VAR_LABELS[9:10]),
                                       selected = c("first_bll", "first_year")),
                    hr(),
                    actionButton("fit_model", "Run Model", class = "btn-primary btn-block")
                ),
                box(width = 8, title = "Results",
                    tabsetPanel(
                      tabPanel("Summary", br(), uiOutput("interpretation")),
                      tabPanel("Model vs Random", br(), uiOutput("comparison")),
                      tabPanel("Odds Ratios", br(),
                               uiOutput("or_note"),
                               plotlyOutput("coef_plot", height = 300)),
                      tabPanel("Diagnostics", br(),
                               p("The ", strong("ROC curve"), " shows discrimination ability. AUC is the probability
                          the model ranks a true positive higher than a true negative. Look for AUC > 0.7.",
                                 style = "font-size: 12px;"),
                               p("The ", strong("calibration plots"), " show whether predicted probabilities match
                          observed rates. Points should fall near the diagonal line.",
                                 style = "font-size: 12px;"),
                               fluidRow(column(4, plotOutput("roc_plot", height = 260)),
                                        column(4, plotOutput("cal_plot_decile", height = 260)),
                                        column(4, plotOutput("cal_plot_full", height = 260)))),
                      tabPanel("Correlations", br(),
                               p("High correlations between predictors can cause unstable estimates. Consider
                          removing one variable if correlation > 0.70.", style = "font-size: 12px;"),
                               p(span("Yellow = borderline (0.50-0.70)", style = "background: #fff3cd; padding: 2px 5px;"),
                                 span("Red = high (> 0.70)", style = "background: #f8d7da; padding: 2px 5px; margin-left: 10px;")),
                               DTOutput("cor_table"),
                               br(),
                               h5("Variance Inflation Factors (VIF)"),
                               p("VIF > 5 suggests problematic multicollinearity. VIF > 10 is severe.", style = "font-size: 12px;"),
                               DTOutput("vif_table")),
                      tabPanel("Technical", br(), verbatimTextOutput("model_summary"))
                    )
                )
              ),
              fluidRow(
                box(width = 12, title = "Predicted Risk Map",
                    p("Average predicted risk for addresses in each area. This map uses the geographic level selected below,
               independent of the random effect setting above.", style = "font-size: 13px; color: #555;"),
                    fluidRow(
                      column(4, radioButtons("pred_type", "Data:", c("Test Set" = "test", "All Addresses" = "full"), inline = TRUE)),
                      column(4, radioButtons("pred_geo", "Map Level:", c("Tract" = "tract", "Block Group" = "bg"), inline = TRUE))
                    ),
                    leafletOutput("pred_map", height = 450))
              )
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
    roc_obj = NULL, youden = NULL, comparison = NULL, cor_matrix = NULL,
    clicked_id = NULL, census_id = NULL, model_run = FALSE, is_scaled = FALSE
  )
  
  # Download handler factory: creates a downloadHandler on output[[id]]
  set_dl <- function(id, data_func, prefix, n_head = NULL) {
    output[[id]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".csv"),
      content = function(file) {
        d <- data_func()
        if (!is.null(n_head)) d <- head(d, n_head())
        write.csv(d, file, row.names = FALSE)
      }
    )
  }
  
  # Geo column lookup for surveillance/tables
  geo_col_for <- function(level) {
    switch(level, "LHD" = "lhd", "County" = "county",
           "Tract" = "tract_geoid", "Block Group" = "bg_geoid")
  }
  
  # ----- SHARED REACTIVE: Surveillance Stats -----
  
  get_stats <- reactive({
    yr  <- input$year_filter
    sup <- input$suppress_n
    ft  <- if (yr == "All") first_test else first_test %>% filter(sample_year == as.numeric(yr))
    
    list(
      lhd    = summarize_lead(ft, "lhd", sup),
      county = summarize_lead(ft, "county", sup),
      tract  = summarize_lead(ft %>% filter(nchar(tract_geoid) == 11), "tract_geoid", sup),
      bg     = summarize_lead(ft %>% filter(nchar(bg_geoid) == 12), "bg_geoid", sup)
    )
  })
  
  # ===== SURVEILLANCE MAP =====
  
  output$surv_map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>% setView(-99.5, 41.5, zoom = 7)
  })
  
  output$state_avg_text <- renderUI({
    unit <- input$unit_type
    var  <- input$fill_var
    yr   <- input$year_filter
    
    avg_all <- if (unit == "Addresses") state_avg_all_addr else state_avg_all
    
    fmt <- function(a, v) {
      switch(v,
             "mean_bll"     = sprintf("%.2f \u00b5g/dL", a$mean_bll),
             "n_elevated"   = format(a$n_elevated, big.mark = ","),
             "n_children"   = format(a$n_children, big.mark = ","),
             "pct_elevated" = sprintf("%.2f%%", a$pct_elevated))
    }
    
    all_txt <- paste0("<b>All Years:</b> ", fmt(avg_all, var))
    
    yr_txt <- ""
    if (yr != "All") {
      trend <- if (unit == "Addresses") state_avgs_addr else state_avgs
      avg_yr <- trend %>% filter(sample_year == as.numeric(yr))
      if (nrow(avg_yr) > 0) yr_txt <- paste0("<br><b>", yr, ":</b> ", fmt(avg_yr, var))
    }
    
    HTML(paste0(all_txt, yr_txt))
  })
  
  output$state_avg_plot <- renderPlotly({
    var    <- input$fill_var
    yr_sel <- input$year_filter
    unit   <- input$unit_type
    
    df <- if (unit == "Addresses") state_avgs_addr else state_avgs
    df <- df %>% mutate(y_val = .data[[var]])
    y_lab <- get_label(var)
    sel_year <- if (yr_sel != "All") as.numeric(yr_sel) else NA
    
    p <- ggplot(df, aes(x = sample_year, y = y_val,
                        text = paste0("Year: ", sample_year, "\n", y_lab, ": ", round(y_val, 2)))) +
      geom_line(color = "#2c7fb8", linewidth = 0.5, alpha = 0.5) +
      geom_point(color = "#2c7fb8", size = 3) +
      labs(x = "Year", y = y_lab, title = paste("Nebraska Trend:", y_lab)) +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold", size = 12))
    
    if (!is.na(sel_year)) {
      p <- p + geom_point(data = df %>% filter(sample_year == sel_year),
                          color = "#d73027", size = 5)
    }
    
    ggplotly(p, tooltip = "text") %>% layout(margin = list(t = 40, b = 40), showlegend = FALSE)
  })
  
  output$trend_summary_table <- renderDT({
    var   <- input$fill_var
    unit  <- input$unit_type
    sup   <- input$suppress_n
    g_col <- geo_col_for(input$geo_level)
    
    trend_df <- if (unit == "Addresses") state_avgs_addr else state_avgs
    
    # Determine which variable to compute StdDev on
    display_var <- var
    if (unit == "Addresses") {
      if (var == "n_children") display_var <- "n_addresses"
      if (var == "n_elevated") display_var <- "addr_elevated"
    }
    
    yr_rows <- lapply(years, function(yr) {
      state_row <- trend_df %>% filter(sample_year == yr)
      if (nrow(state_row) == 0) return(NULL)
      
      ft_yr <- first_test %>% filter(sample_year == yr)
      geo_stats <- summarize_lead(ft_yr, g_col, sup)
      
      vals <- geo_stats[[display_var]][!is.na(geo_stats[[display_var]]) & !geo_stats$suppressed]
      n_miss <- sum(is.na(geo_stats[[display_var]]) | geo_stats$suppressed)
      
      data.frame(Year = yr, Measure = round(state_row[[var]], 2),
                 `Std Dev` = if (length(vals) > 1) round(sd(vals, na.rm = TRUE), 2) else NA_real_,
                 `Missing Areas` = n_miss, check.names = FALSE)
    })
    
    datatable(do.call(rbind, yr_rows),
              options = list(dom = 't', pageLength = 20, scrollY = "150px"), rownames = FALSE)
  })
  
  observeEvent(c(input$geo_level, input$fill_var, input$unit_type, input$year_filter,
                 input$suppress_n, input$use_quintile), {
                   stats <- get_stats()
                   var   <- input$fill_var
                   unit  <- input$unit_type
                   
                   display_var <- var
                   if (unit == "Addresses") {
                     if (var == "n_children") display_var <- "n_addresses"
                     if (var == "n_elevated") display_var <- "addr_elevated"
                   }
                   
                   stat_data <- switch(input$geo_level, "LHD" = stats$lhd, "County" = stats$county,
                                       "Tract" = stats$tract, "Block Group" = stats$bg)
                   if (is.null(stat_data) || nrow(stat_data) == 0) return()
                   
                   geo <- switch(input$geo_level,
                                 "LHD"         = if (!is.null(lhd_bounds))  lhd_bounds  %>% left_join(stat_data, by = c("lhd" = "geo_id"))  else return(),
                                 "County"      = if (!is.null(ne_counties)) ne_counties %>% left_join(stat_data, by = c("NAME" = "geo_id")) else return(),
                                 "Tract"       = if (!is.null(ne_tracts))   ne_tracts   %>% left_join(stat_data, by = c("GEOID" = "geo_id")) else return(),
                                 "Block Group" = if (!is.null(ne_bgs))      ne_bgs      %>% left_join(stat_data, by = c("GEOID" = "geo_id")) else return()
                   )
                   
                   id_col    <- switch(input$geo_level, "LHD" = "lhd", "County" = "NAME", "Tract" = "GEOID", "Block Group" = "GEOID")
                   geo_label <- switch(input$geo_level, "LHD" = "Health Department", "County" = "County",
                                       "Tract" = "Census Tract", "Block Group" = "Block Group")
                   if (!display_var %in% names(geo)) return()
                   
                   vals <- geo[[display_var]][!is.na(geo[[display_var]]) & is.finite(geo[[display_var]]) &
                                                (is.null(geo$suppressed) | !geo$suppressed)]
                   if (length(vals) == 0) vals <- c(0, 1)
                   
                   var_label <- get_label(var)
                   
                   if (input$use_quintile) {
                     q <- quantile(vals, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
                     colors <- c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026")
                     
                     geo$fill_color <- sapply(seq_len(nrow(geo)), function(i) {
                       if (!is.null(geo$suppressed) && !is.na(geo$suppressed[i]) && geo$suppressed[i]) return("#666666")
                       v <- geo[[display_var]][i]
                       if (is.na(v) || !is.finite(v)) return("#cccccc")
                       bin <- findInterval(v, q, rightmost.closed = TRUE)
                       colors[max(1, min(bin, 5))]
                     })
                     
                     legend_labels <- c(paste0("Q1: ", round(q[1], 1), "-", round(q[2], 1)),
                                        paste0("Q2: ", round(q[2], 1), "-", round(q[3], 1)),
                                        paste0("Q3: ", round(q[3], 1), "-", round(q[4], 1)),
                                        paste0("Q4: ", round(q[4], 1), "-", round(q[5], 1)),
                                        paste0("Q5: ", round(q[5], 1), "-", round(q[6], 1)),
                                        "Suppressed", "No Data")
                     legend_colors <- c(colors, "#666666", "#cccccc")
                   } else {
                     if (var == "mean_bll") {
                       domain <- c(0, 5)
                     } else if (var == "pct_elevated") {
                       domain <- c(0, 30)
                     } else {
                       domain <- c(0, max(vals, na.rm = TRUE))
                     }
                     pal <- colorNumeric("YlOrRd", domain = domain, na.color = "#cccccc")
                     
                     geo$fill_color <- sapply(seq_len(nrow(geo)), function(i) {
                       if (!is.null(geo$suppressed) && !is.na(geo$suppressed[i]) && geo$suppressed[i]) return("#666666")
                       v <- geo[[display_var]][i]
                       if (is.na(v) || !is.finite(v)) return("#cccccc")
                       pal(min(v, domain[2]))
                     })
                     
                     breaks <- seq(domain[1], domain[2], length.out = 6)
                     legend_labels <- c(sapply(1:5, function(i) paste0(round(breaks[i], 1), "-", round(breaks[i+1], 1))),
                                        "Suppressed", "No Data")
                     legend_colors <- c(pal(breaks[1:5] + diff(breaks[1:2])/2), "#666666", "#cccccc")
                   }
                   
                   geo$tip <- sapply(seq_len(nrow(geo)), function(i) {
                     sup_flag <- if (!is.null(geo$suppressed) && !is.na(geo$suppressed[i]) && geo$suppressed[i]) " [Suppressed]" else ""
                     
                     if (unit == "Children") {
                       paste0("<b>", geo_label, ": ", geo[[id_col]][i], "</b>", sup_flag, "<br>",
                              "Elevated Count: ", geo$n_elevated[i], "<br>",
                              "Total Tested: ", format(geo$n_children[i], big.mark = ","), "<br>",
                              "Percent Elevated: ", if (!is.na(geo$pct_elevated[i])) sprintf("%.1f%%", geo$pct_elevated[i]) else "N/A", "<br>",
                              "Average BLL: ", if (!is.na(geo$mean_bll[i])) sprintf("%.2f \u00b5g/dL", geo$mean_bll[i]) else "N/A")
                     } else {
                       paste0("<b>", geo_label, ": ", geo[[id_col]][i], "</b>", sup_flag, "<br>",
                              "Addresses Elevated: ", geo$addr_elevated[i], "<br>",
                              "Addresses Tested: ", format(geo$n_addresses[i], big.mark = ","), "<br>",
                              "Percent Elevated: ", if (!is.na(geo$n_addresses[i]) && geo$n_addresses[i] > 0)
                                sprintf("%.1f%%", 100 * geo$addr_elevated[i] / geo$n_addresses[i]) else "N/A", "<br>",
                              "Average BLL: ", if (!is.na(geo$mean_bll[i])) sprintf("%.2f \u00b5g/dL", geo$mean_bll[i]) else "N/A")
                     }
                   })
                   
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
    
    geo_label <- switch(input$geo_level, "LHD" = "Health Dept", "County" = "County",
                        "Tract" = "Tract", "Block Group" = "Block Group")
    cat(geo_label, ":", rv$clicked_id, "\n\n")
    
    stats <- get_stats()
    d <- switch(input$geo_level, "LHD" = stats$lhd, "County" = stats$county,
                "Tract" = stats$tract, "Block Group" = stats$bg) %>%
      filter(geo_id == rv$clicked_id)
    if (nrow(d) == 0) return(cat("No data"))
    
    if (d$suppressed) cat("[Suppressed]\n\n")
    
    if (input$unit_type == "Children") {
      cat("Elevated Count:", d$n_elevated, "\n")
      cat("Total Tested:", format(d$n_children, big.mark = ","), "\n")
      cat("Percent Elevated:", if (!is.na(d$pct_elevated)) sprintf("%.2f%%", d$pct_elevated) else "N/A", "\n")
      cat("Average BLL:", if (!is.na(d$mean_bll)) sprintf("%.2f \u00b5g/dL", d$mean_bll) else "N/A")
    } else {
      cat("Addresses Elevated:", d$addr_elevated, "\n")
      cat("Addresses Tested:", format(d$n_addresses, big.mark = ","), "\n")
      pct <- if (d$n_addresses > 0) sprintf("%.2f%%", 100 * d$addr_elevated / d$n_addresses) else "N/A"
      cat("Percent Elevated:", pct, "\n")
      cat("Average BLL:", if (!is.na(d$mean_bll)) sprintf("%.2f \u00b5g/dL", d$mean_bll) else "N/A")
    }
  })
  
  # ===== PRIORITY ADDRESSES =====
  
  output$ranked_map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>%
      setView(-99.5, 41.5, zoom = 7) %>%
      addPolygons(data = lhd_bounds, fillColor = "transparent", color = "#2c7fb8", weight = 2, label = ~lhd)
  })
  
  observe({
    tiles <- if (input$map_tiles == "sat") providers$Esri.WorldImagery else providers$CartoDB.Positron
    leafletProxy("ranked_map") %>% clearTiles() %>% addProviderTiles(tiles)
  })
  
  get_ranked_addresses <- reactive({
    addr <- address_stats %>% filter(!is.na(lhd))
    if (input$ranked_lhd != "All") addr <- addr %>% filter(lhd == input$ranked_lhd)
    
    if (!is.null(rv$full_preds)) {
      addr <- addr %>% left_join(rv$full_preds %>% select(address_id, pred_prob), by = "address_id")
    } else {
      addr$pred_prob <- NA_real_
    }
    
    addr %>%
      group_by(lhd) %>%
      arrange(desc(.data[[input$rank_by]]), desc(pred_prob), .by_group = TRUE) %>%
      mutate(rank = row_number()) %>%
      ungroup() %>%
      filter(rank <= input$top_n) %>%
      mutate(
        color = get_bll_color(mean_bll),
        risk_decile = ifelse(!is.na(pred_prob), ntile(pred_prob, 10), NA_integer_)
      )
  })
  
  observeEvent(c(input$ranked_lhd, input$top_n, input$rank_by, rv$full_preds), {
    addr <- get_ranked_addresses()
    if (nrow(addr) == 0) return()
    
    addr$tip <- paste0(
      "<b>#", addr$rank, " - ", addr$street_address, "</b><br>",
      "Health Dept: ", addr$lhd, "<br>",
      "County: ", addr$county, "<br>",
      "Children Tested: ", addr$n_children, "<br>",
      "Children Elevated: ", addr$n_elevated, "<br>",
      "Average BLL: ", round(addr$mean_bll, 2), " \u00b5g/dL<br>",
      "Max BLL: ", round(addr$max_bll, 2), " \u00b5g/dL<br>",
      "First Year: ", addr$first_year, " | Last Year: ", addr$last_year,
      ifelse(!is.na(addr$pred_prob), paste0("<br>Risk Decile: ", addr$risk_decile,
                                            " | Risk: ", round(addr$pred_prob * 100, 1), "%"), "")
    )
    
    leafletProxy("ranked_map") %>% clearMarkers() %>%
      addCircleMarkers(data = addr, lng = ~lng, lat = ~lat, radius = 8,
                       fillColor = ~color, color = "#333", weight = 1, fillOpacity = 0.85,
                       label = lapply(addr$tip, HTML)) %>%
      addLabelOnlyMarkers(data = addr, lng = ~lng, lat = ~lat,
                          label = as.character(addr$rank),
                          labelOptions = labelOptions(noHide = TRUE, direction = "center",
                                                      textOnly = TRUE, textsize = "10px",
                                                      style = list("font-weight" = "bold", "color" = "white")))
  }, ignoreNULL = FALSE)
  
  get_ranked_table_data <- reactive({
    addr <- get_ranked_addresses()
    if (nrow(addr) == 0) return(NULL)
    
    if (!is.null(rv$fit) && !is.null(rv$full_preds)) {
      coefs <- tidy(rv$fit, effects = "fixed") %>% filter(term != "(Intercept)")
      
      addr <- addr %>%
        left_join(rv$full_preds %>% select(-pred_prob), by = "address_id", suffix = c("", ".pred"))
      
      addr$`Risk Explanation` <- sapply(seq_len(nrow(addr)), function(i) {
        if (is.na(addr$pred_prob[i])) return("")
        
        parts <- sapply(seq_len(nrow(coefs)), function(j) {
          var <- gsub("scale\\(|\\)", "", coefs$term[j])
          if (var %in% names(addr) && !is.na(addr[[var]][i])) {
            val <- addr[[var]][i]
            lab <- get_label(var)
            if (grepl("pct_|pct$", var)) paste0(lab, ": ", round(val, 1), "%")
            else if (var == "median_income") paste0(lab, ": $", format(round(val), big.mark = ","))
            else if (var == "first_bll") paste0(lab, ": ", round(val, 1))
            else if (var == "first_year") paste0(lab, ": ", round(val))
            else paste0(lab, ": ", round(val, 2))
          } else NULL
        })
        parts <- parts[!sapply(parts, is.null)]
        if (length(parts) > 0) paste(parts, collapse = ", ") else "Geographic effect only"
      })
    } else {
      addr$`Risk Explanation` <- ""
    }
    
    addr %>%
      arrange(lhd, rank) %>%
      transmute(Rank = rank, `Health Dept` = lhd, County = county, Address = street_address,
                Tested = n_children, Elevated = n_elevated,
                `Average BLL` = round(mean_bll, 2), `Max BLL` = round(max_bll, 2),
                `First Year` = first_year, `Last Year` = last_year,
                `Risk Decile` = risk_decile,
                `Risk %` = ifelse(!is.na(pred_prob), round(pred_prob * 100, 1), NA_real_),
                `Risk Explanation`)
  })
  
  output$ranked_table <- renderDT({
    d <- get_ranked_table_data()
    if (is.null(d)) return(NULL)
    datatable(d, options = list(dom = 't', pageLength = 500, scrollY = "230px", scrollCollapse = TRUE,
                                columnDefs = list(list(width = '200px', targets = 12))),
              rownames = FALSE)
  })
  
  # ===== CENSUS MAP =====
  
  output$census_map <- renderLeaflet({
    leaflet() %>% addProviderTiles(providers$CartoDB.Positron) %>% setView(-99.5, 41.5, zoom = 7)
  })
  
  get_census_geo_data <- reactive({
    var <- input$census_var
    
    if (input$census_geo == "tract") {
      req(tract_census, ne_tracts)
      geo <- ne_tracts %>% left_join(tract_census, by = c("GEOID" = "tract_geoid"))
      test_data <- first_test %>% filter(!is.na(tract_geoid)) %>%
        group_by(GEOID = tract_geoid) %>%
        summarize(n_tested = n(), n_elev = sum(elevated, na.rm = TRUE),
                  pct_elev = 100 * mean(elevated, na.rm = TRUE),
                  avg_bll = mean(result, na.rm = TRUE), .groups = "drop")
    } else {
      req(bg_census, ne_bgs)
      geo <- ne_bgs %>% left_join(bg_census, by = c("GEOID" = "bg_geoid"))
      test_data <- first_test %>% filter(!is.na(bg_geoid), nchar(bg_geoid) == 12) %>%
        group_by(GEOID = bg_geoid) %>%
        summarize(n_tested = n(), n_elev = sum(elevated, na.rm = TRUE),
                  pct_elev = 100 * mean(elevated, na.rm = TRUE),
                  avg_bll = mean(result, na.rm = TRUE), .groups = "drop")
    }
    
    geo %>% left_join(test_data, by = "GEOID")
  })
  
  observeEvent(c(input$census_var, input$census_geo, input$census_quintile), {
    var <- input$census_var
    geo <- get_census_geo_data()
    if (is.null(geo) || !var %in% names(geo)) return()
    
    vals <- as.numeric(geo[[var]])
    vals_clean <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals_clean) == 0) return()
    
    geo_label <- if (input$census_geo == "tract") "Census Tract" else "Block Group"
    var_label <- get_label(var)
    
    fmt_v <- function(v) {
      if (is.na(v) || !is.finite(v)) return("N/A")
      if (var == "median_income") return(paste0("$", format(round(v), big.mark = ",")))
      if (var == "housing_units") return(format(round(v), big.mark = ","))
      sprintf("%.1f", v)
    }
    
    geo$tip <- paste0(
      "<b>", geo_label, ": ", geo$GEOID, "</b><br>",
      var_label, ": ", sapply(vals, fmt_v), "<br>",
      "Total Tested: ", ifelse(is.na(geo$n_tested), "N/A", format(geo$n_tested, big.mark = ",")), "<br>",
      "Elevated Count: ", ifelse(is.na(geo$n_elev), "N/A", geo$n_elev), "<br>",
      "Percent Elevated: ", ifelse(is.na(geo$pct_elev), "N/A", sprintf("%.1f%%", geo$pct_elev)), "<br>",
      "Average BLL: ", ifelse(is.na(geo$avg_bll), "N/A", sprintf("%.2f \u00b5g/dL", geo$avg_bll))
    )
    
    if (input$census_quintile) {
      q <- quantile(vals_clean, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
      colors <- c("#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
      
      geo$quintile <- sapply(vals, function(v) {
        if (is.na(v) || !is.finite(v)) return(NA_integer_)
        findInterval(v, q, rightmost.closed = TRUE)
      })
      geo$fill_color <- sapply(geo$quintile, function(q_num) {
        if (is.na(q_num)) "#cccccc" else colors[max(1, min(q_num, 5))]
      })
      
      legend_labels <- c(paste0("Q1: ", round(q[1], 1), "-", round(q[2], 1)),
                         paste0("Q2: ", round(q[2], 1), "-", round(q[3], 1)),
                         paste0("Q3: ", round(q[3], 1), "-", round(q[4], 1)),
                         paste0("Q4: ", round(q[4], 1), "-", round(q[5], 1)),
                         paste0("Q5: ", round(q[5], 1), "-", round(q[6], 1)),
                         "No Data")
      legend_colors <- c(colors, "#cccccc")
    } else {
      min_v <- min(vals_clean, na.rm = TRUE)
      max_v <- max(vals_clean, na.rm = TRUE)
      
      pal <- colorNumeric("Blues", domain = c(min_v, max_v), na.color = "#cccccc")
      geo$fill_color <- sapply(vals, function(v) {
        if (is.na(v) || !is.finite(v)) "#cccccc" else pal(v)
      })
      
      breaks <- seq(min_v, max_v, length.out = 6)
      legend_labels <- c(sapply(1:5, function(i) paste0(round(breaks[i], 1), "-", round(breaks[i+1], 1))), "No Data")
      legend_colors <- c(pal(breaks[1:5] + diff(breaks[1:2])/2), "#cccccc")
    }
    
    leafletProxy("census_map") %>% clearShapes() %>% clearControls() %>%
      addPolygons(data = geo, fillColor = ~fill_color, fillOpacity = 0.7,
                  color = "white", weight = 0.5, layerId = ~GEOID,
                  label = lapply(geo$tip, HTML)) %>%
      addLegend("bottomright", colors = legend_colors, labels = legend_labels,
                title = var_label, opacity = 0.8)
  }, ignoreNULL = FALSE)
  
  observeEvent(input$census_map_shape_click, { rv$census_id <- input$census_map_shape_click$id })
  
  output$census_info <- renderPrint({
    req(rv$census_id)
    geo_id <- rv$census_id
    
    d <- if (input$census_geo == "tract") tract_census %>% filter(tract_geoid == geo_id)
    else bg_census %>% filter(bg_geoid == geo_id)
    if (nrow(d) == 0) return(cat("No census data"))
    
    test_d <- if (input$census_geo == "tract") {
      first_test %>% filter(tract_geoid == geo_id)
    } else {
      first_test %>% filter(bg_geoid == geo_id)
    }
    
    cat(if (input$census_geo == "tract") "Tract" else "Block Group", ":", geo_id, "\n\n")
    
    cat("--- Lead Testing ---\n")
    cat("Total Tested:", nrow(test_d), "\n")
    cat("Elevated Count:", sum(test_d$elevated, na.rm = TRUE), "\n")
    cat("Percent Elevated:", if (nrow(test_d) > 0) sprintf("%.1f%%", 100 * mean(test_d$elevated, na.rm = TRUE)) else "N/A", "\n")
    cat("Average BLL:", if (nrow(test_d) > 0) sprintf("%.2f \u00b5g/dL", mean(test_d$result, na.rm = TRUE)) else "N/A", "\n\n")
    
    cat("--- Census Variables ---\n")
    cat("Pre-1980 Housing:", if (!is.na(d$pct_pre1980)) sprintf("%.1f%%", d$pct_pre1980) else "N/A", "\n")
    cat("Poverty Rate:", if (!is.na(d$pct_poverty) && is.finite(d$pct_poverty)) sprintf("%.1f%%", d$pct_poverty) else "N/A", "\n")
    cat("Renter-Occupied:", if (!is.na(d$pct_renter)) sprintf("%.1f%%", d$pct_renter) else "N/A", "\n")
    cat("Median Income:", if (!is.na(d$median_income)) paste0("$", format(round(d$median_income), big.mark = ",")) else "N/A", "\n")
    cat("Non-White:", if (!is.na(d$pct_nonwhite)) sprintf("%.1f%%", d$pct_nonwhite) else "N/A", "\n")
    cat("Housing Units:", if (!is.na(d$housing_units)) format(round(d$housing_units), big.mark = ",") else "N/A", "\n")
    cat("Housing Density:", if (!is.na(d$log_housing_density)) sprintf("%.2f", d$log_housing_density) else "N/A", "\n")
    cat("Children <6:", if (!is.na(d$pct_children_u6)) sprintf("%.1f%%", d$pct_children_u6) else "N/A")
  })
  
  get_quintile_data <- reactive({
    var <- input$census_var
    
    census_data <- if (input$census_geo == "tract" && !is.null(tract_census)) {
      tract_census %>% rename(geoid = tract_geoid)
    } else if (!is.null(bg_census)) {
      bg_census %>% rename(geoid = bg_geoid)
    } else return(NULL)
    
    if (!var %in% names(census_data)) return(NULL)
    
    test_geo <- if (input$census_geo == "tract") {
      first_test %>% filter(!is.na(tract_geoid)) %>% rename(geoid = tract_geoid)
    } else {
      first_test %>% filter(!is.na(bg_geoid), nchar(bg_geoid) == 12) %>% rename(geoid = bg_geoid)
    }
    
    vals <- census_data[[var]][!is.na(census_data[[var]]) & is.finite(census_data[[var]])]
    if (length(vals) == 0) return(NULL)
    q <- quantile(vals, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
    
    census_data <- census_data %>%
      filter(!is.na(.data[[var]]), is.finite(.data[[var]])) %>%
      mutate(q_num = findInterval(.data[[var]], q, rightmost.closed = TRUE),
             q_num = pmax(1, pmin(q_num, 5)))
    
    test_geo %>%
      inner_join(census_data %>% select(geoid, q_num), by = "geoid") %>%
      group_by(q_num) %>%
      summarize(Areas = n_distinct(geoid), Tested = n(), Elevated = sum(elevated, na.rm = TRUE),
                `% Elevated` = round(100 * mean(elevated, na.rm = TRUE), 2),
                `Average BLL` = round(mean(result, na.rm = TRUE), 2), .groups = "drop") %>%
      mutate(Quintile = paste0("Q", q_num, " (", round(q[q_num], 1), "-", round(q[q_num + 1], 1), ")")) %>%
      select(Quintile, Areas, Tested, Elevated, `% Elevated`, `Average BLL`)
  })
  
  output$census_quintile_table <- renderDT({
    d <- get_quintile_data()
    if (is.null(d)) return(NULL)
    datatable(d, options = list(dom = 't', pageLength = 5), rownames = FALSE)
  })
  
  get_area_detail_data <- reactive({
    var <- input$census_var
    geo <- get_census_geo_data()
    if (is.null(geo)) return(NULL)
    
    vals <- as.numeric(geo[[var]])
    vals_clean <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals_clean) == 0) return(NULL)
    q <- quantile(vals_clean, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
    
    geo_df <- as.data.frame(geo) %>%
      select(GEOID, all_of(var), pct_elev, avg_bll) %>%
      filter(!is.na(.data[[var]])) %>%
      mutate(
        q_num = findInterval(.data[[var]], q, rightmost.closed = TRUE),
        q_num = pmax(1, pmin(q_num, 5)),
        Quintile = paste0("Q", q_num)
      ) %>%
      transmute(
        Area = GEOID,
        `Census Value` = round(.data[[var]], 2),
        Quintile,
        `% Elevated` = round(pct_elev, 2),
        `Average BLL` = round(avg_bll, 2)
      )
    
    geo_df
  })
  
  output$census_area_table <- renderDT({
    d <- get_area_detail_data()
    if (is.null(d)) return(NULL)
    datatable(d, options = list(pageLength = 10, scrollY = "200px"), rownames = FALSE)
  })
  
  # ===== DATA TABLES =====
  
  get_table_data <- reactive({
    yr_range <- input$table_years
    sup <- input$suppress_n
    ft <- first_test %>% filter(sample_year >= yr_range[1], sample_year <= yr_range[2])
    
    geo_data <- if (input$table_geo == "tract") {
      summarize_lead(ft %>% filter(nchar(tract_geoid) == 11), "tract_geoid", sup)
    } else {
      summarize_lead(ft %>% filter(nchar(bg_geoid) == 12), "bg_geoid", sup)
    }
    
    list(
      lhd    = summarize_lead(ft, "lhd", sup),
      county = summarize_lead(ft, "county", sup),
      geo    = geo_data
    )
  })
  
  output$lhd_table <- renderDT({
    d <- get_table_data()$lhd
    if (is.null(d) || nrow(d) == 0) return(NULL)
    datatable(format_table(d, "Health Department", unit = input$table_unit),
              options = list(pageLength = as.integer(input$table_rows)), rownames = FALSE)
  })
  
  output$county_table <- renderDT({
    d <- get_table_data()$county
    if (is.null(d) || nrow(d) == 0) return(NULL)
    datatable(format_table(d, "County", unit = input$table_unit),
              options = list(pageLength = as.integer(input$table_rows)), rownames = FALSE)
  })
  
  output$geo_table <- renderDT({
    d <- get_table_data()$geo
    if (is.null(d) || nrow(d) == 0) return(NULL)
    datatable(format_table(d, "Small Area", unit = input$table_unit, show_county = TRUE),
              options = list(pageLength = as.integer(input$table_rows)), rownames = FALSE)
  })
  
  get_addr_table <- reactive({
    yr_range <- input$table_years
    addr <- address_stats %>% filter(first_year <= yr_range[2], last_year >= yr_range[1])
    
    if (!is.null(rv$full_preds)) {
      addr <- addr %>% left_join(rv$full_preds %>% select(address_id, pred_prob), by = "address_id")
    } else {
      addr$pred_prob <- NA_real_
    }
    
    addr %>%
      mutate(risk_decile = ifelse(!is.na(pred_prob), ntile(pred_prob, 10), NA_integer_)) %>%
      transmute(Address = street_address, County = county, `Health Dept` = lhd,
                Tested = n_children, Elevated = n_elevated,
                `Average BLL` = round(mean_bll, 2), `Max BLL` = round(max_bll, 2),
                `First Year` = first_year, `Last Year` = last_year,
                `Risk Decile` = risk_decile,
                `Risk %` = ifelse(!is.na(pred_prob), round(pred_prob * 100, 2), NA_real_))
  })
  
  output$addr_table <- renderDT({
    datatable(get_addr_table(), options = list(pageLength = as.integer(input$table_rows), scrollX = TRUE), rownames = FALSE)
  })
  
  # --- All download handlers via factory ---
  n_rows <- reactive(as.integer(input$table_rows))
  
  set_dl("dl_ranked",      get_ranked_table_data, "Priority_Addresses")
  set_dl("dl_quintile",    get_quintile_data, "Quintile_Summary")
  set_dl("dl_area_detail", get_area_detail_data, "Area_Details")
  
  set_dl("dl_lhd_disp", function() format_table(get_table_data()$lhd, "Health Department", unit = input$table_unit),
         "LHD_displayed", n_head = n_rows)
  set_dl("dl_lhd_all", function() format_table(get_table_data()$lhd, "Health Department", unit = input$table_unit),
         "LHD_all")
  set_dl("dl_county_disp", function() format_table(get_table_data()$county, "County", unit = input$table_unit),
         "County_displayed", n_head = n_rows)
  set_dl("dl_county_all", function() format_table(get_table_data()$county, "County", unit = input$table_unit),
         "County_all")
  set_dl("dl_geo_disp", function() format_table(get_table_data()$geo, "Small Area", unit = input$table_unit, show_county = TRUE),
         "SmallArea_displayed", n_head = n_rows)
  set_dl("dl_geo_all", function() format_table(get_table_data()$geo, "Small Area", unit = input$table_unit, show_county = TRUE),
         "SmallArea_all")
  set_dl("dl_addr_disp", get_addr_table, "Addresses_displayed", n_head = n_rows)
  set_dl("dl_addr_all",  get_addr_table, "Addresses_all")
  
  # ===== RISK MODEL =====
  
  run_model <- function() {
    req(model_data)
    
    showNotification("Fitting model...", id = "fit_msg", duration = NULL)
    
    tryCatch({
      re_var <- input$random_effect
      census_data <- if (re_var == "tract_geoid" && !is.null(tract_census)) tract_census else bg_census
      join_var <- if (re_var == "tract_geoid" && !is.null(tract_census)) "tract_geoid" else "bg_geoid"
      if (join_var == "bg_geoid") re_var <- "bg_geoid"
      
      if (is.null(census_data)) {
        removeNotification("fit_msg")
        showNotification("Census data not available", type = "error")
        return()
      }
      
      md <- model_data %>% left_join(census_data, by = join_var)
      
      covars <- c(input$census_covars, input$addr_covars)
      census_covars_selected <- intersect(input$census_covars, names(census_data))
      if (length(census_covars_selected) > 0) {
        n_matched <- sum(!is.na(md[[census_covars_selected[1]]]))
        if (n_matched == 0) {
          removeNotification("fit_msg")
          showNotification(
            paste0("Census join failed: 0 of ", nrow(md), " rows matched. ",
                   "Re-run data_prep.R to rebuild all files in one pass."),
            type = "error", duration = 10)
          return()
        }
      }
      
      if (input$train_all) {
        train <- md; test <- md
      } else {
        train <- md %>% filter(first_year >= input$train_years[1], first_year <= input$train_years[2])
        test  <- md %>% filter(first_year > input$train_years[2])
      }
      
      if (length(covars) > 0) {
        for (cv in covars) {
          if (cv %in% names(train)) {
            train <- train %>% filter(!is.na(.data[[cv]]) & is.finite(.data[[cv]]))
            test  <- test  %>% filter(!is.na(.data[[cv]]) & is.finite(.data[[cv]]))
          }
        }
      }
      
      train <- train %>% filter(!is.na(outcome), !is.na(.data[[re_var]]))
      test  <- test  %>% filter(!is.na(outcome), !is.na(.data[[re_var]]))
      
      if (nrow(train) < 100) {
        removeNotification("fit_msg")
        showNotification("Not enough training data.", type = "error")
        return()
      }
      if (nrow(test) < 50 && !input$train_all) {
        removeNotification("fit_msg")
        showNotification("Not enough test data.", type = "error")
        return()
      }
      
      covars_in_data <- intersect(covars, names(train))
      is_scaled <- input$scale_vars
      rv$is_scaled <- is_scaled
      
      if (length(covars_in_data) > 0) {
        terms <- if (is_scaled) paste0("scale(", covars_in_data, ")") else covars_in_data
        formula <- as.formula(paste("outcome ~", paste(terms, collapse = " + "), "+ (1 |", re_var, ")"))
      } else {
        formula <- as.formula(paste("outcome ~ 1 + (1 |", re_var, ")"))
      }
      
      fit <- suppressWarnings(glmer(formula, data = train, family = binomial,
                                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))))
      
      fit_full <- suppressWarnings(glmer(formula, data = md %>% filter(!is.na(.data[[re_var]])),
                                         family = binomial,
                                         control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))))
      
      test$pred_prob <- predict(fit, newdata = test, type = "response", allow.new.levels = TRUE)
      
      addr_pred <- address_stats %>% left_join(census_data, by = join_var)
      can_pred <- !is.na(addr_pred[[re_var]])
      if (length(covars_in_data) > 0) {
        for (cv in covars_in_data) {
          can_pred <- can_pred & !is.na(addr_pred[[cv]]) & is.finite(addr_pred[[cv]])
        }
      }
      addr_pred$pred_prob <- NA_real_
      if (sum(can_pred) > 0) {
        addr_pred$pred_prob[can_pred] <- predict(fit_full, newdata = addr_pred[can_pred, ],
                                                 type = "response", allow.new.levels = TRUE)
      }
      
      roc_obj <- roc(test$outcome, test$pred_prob, quiet = TRUE)
      coords_all <- coords(roc_obj, "all", ret = c("threshold", "sensitivity", "specificity"))
      youden_idx <- which.max(coords_all$sensitivity + coords_all$specificity - 1)
      youden <- list(threshold = coords_all$threshold[youden_idx],
                     sens = coords_all$sensitivity[youden_idx],
                     spec = coords_all$specificity[youden_idx])
      
      n_test <- nrow(test)
      n_elev <- sum(test$outcome == 1)
      screen_n <- round(n_test * 0.5)
      test_ranked <- test %>% arrange(desc(pred_prob))
      model_det <- sum(test_ranked$outcome[1:screen_n] == 1)
      random_det <- round(n_elev * 0.5)
      
      tests_per_pos_model <- screen_n / max(model_det, 1)
      tests_per_pos_random <- screen_n / max(random_det, 1)
      
      rv$comparison <- list(n_test = n_test, n_elev = n_elev, screen_n = screen_n,
                            model_det = model_det, random_det = random_det,
                            model_pct = round(100 * model_det / max(n_elev, 1), 1),
                            random_pct = 50,
                            improve = round(100 * (model_det - random_det) / max(random_det, 1), 1),
                            auc = round(as.numeric(auc(roc_obj)), 3),
                            tests_per_pos_model = round(tests_per_pos_model, 1),
                            tests_per_pos_random = round(tests_per_pos_random, 1))
      
      if (length(covars_in_data) > 1) {
        rv$cor_matrix <- cor(train[covars_in_data], use = "pairwise.complete.obs")
      } else {
        rv$cor_matrix <- NULL
      }
      
      rv$fit <- fit; rv$fit_full <- fit_full; rv$test_preds <- test
      rv$full_preds <- addr_pred; rv$roc_obj <- roc_obj; rv$youden <- youden
      
      removeNotification("fit_msg")
      showNotification(paste("Model complete!", sum(!is.na(addr_pred$pred_prob)), "addresses predicted."),
                       type = "message", duration = 4)
      
    }, error = function(e) {
      removeNotification("fit_msg")
      showNotification(paste("Error:", e$message), type = "error", duration = 8)
    })
  }
  
  observeEvent(input$fit_model, { run_model() })
  
  observe({
    if (!rv$model_run) {
      rv$model_run <- TRUE
      isolate(run_model())
    }
  })
  
  output$interpretation <- renderUI({
    if (is.null(rv$fit)) return(div(class = "metric-box", p("Click 'Run Model' to fit.")))
    
    coefs <- tidy(rv$fit, effects = "fixed", conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(or = exp(estimate), var_name = gsub("scale\\(|\\)", "", term))
    
    is_scaled <- rv$is_scaled
    unit_text <- if (is_scaled) "1 standard deviation" else "1 unit"
    
    findings <- if (nrow(coefs) > 0) {
      lapply(1:nrow(coefs), function(i) {
        r <- coefs[i, ]
        pct_change <- round((r$or - 1) * 100, 0)
        abs_pct <- abs(pct_change)
        
        tags$li(
          strong(get_label(r$var_name)), ": ",
          "A ", unit_text, " ", if (pct_change > 0) "increase" else "decrease", " is associated with ",
          abs_pct, "% ", if (r$or > 1) "higher" else "lower", " odds of elevated BLL ",
          "(OR = ", round(r$or, 2), ")"
        )
      })
    } else {
      list(tags$li("Model uses geographic grouping only."))
    }
    
    div(class = "metric-box",
        h4("What the Model Found:"),
        tags$ul(findings),
        if (is_scaled) p(em("Note: Variables were standardized. ORs represent 1 SD change, not 1-unit change."),
                         style = "font-size: 11px; color: #666;"),
        hr(),
        h4("Model Performance:"),
        p(strong("AUC = ", rv$comparison$auc)),
        p("The AUC (Area Under the ROC Curve) measures how well the model discriminates between addresses
           that will vs. won't have elevated children. Values range from 0.5 (no better than chance) to 1.0 (perfect).
           An AUC of 0.7+ is generally considered acceptable; 0.8+ is good."),
        p(if (rv$comparison$auc >= 0.8) "This model shows good discrimination."
          else if (rv$comparison$auc >= 0.7) "This model shows acceptable discrimination."
          else "This model shows modest discrimination. Consider adding more predictors."))
  })
  
  output$comparison <- renderUI({
    if (is.null(rv$comparison)) return(div(class = "metric-box", p("Run model first.")))
    c <- rv$comparison
    
    div(
      div(class = "compare-box",
          h4("Screening Efficiency: Model vs Random", style = "margin-top: 0;"),
          p("If you could only test 50% of addresses (", format(c$screen_n, big.mark = ","), " addresses):"),
          hr(),
          fluidRow(
            column(6,
                   div(style = "text-align: center; padding: 15px; background: #e8f4f8; border-radius: 6px;",
                       h5("Using This Model", style = "color: #2c7fb8; margin: 0;"),
                       h2(c$model_det, style = "margin: 10px 0; color: #2c7fb8;"),
                       p("elevated addresses found"),
                       p(strong(c$model_pct, "%"), " of all elevated"),
                       hr(),
                       p(strong(c$tests_per_pos_model), " tests per positive")
                   )),
            column(6,
                   div(style = "text-align: center; padding: 15px; background: #f5f5f5; border-radius: 6px;",
                       h5("Random Selection", style = "color: #666; margin: 0;"),
                       h2(c$random_det, style = "margin: 10px 0; color: #666;"),
                       p("elevated addresses found"),
                       p(strong("50%"), " of all elevated"),
                       hr(),
                       p(strong(c$tests_per_pos_random), " tests per positive")
                   ))
          ),
          hr(),
          div(style = "text-align: center;",
              h4(paste0("+", c$improve, "% more elevated addresses found with model"),
                 style = "color: #4a9; margin: 0;"),
              p(em("Lower tests per positive = more efficient screening")))
      )
    )
  })
  
  output$or_note <- renderUI({
    if (is.null(rv$fit)) return(NULL)
    
    div(class = "warning-box",
        if (rv$is_scaled) {
          "These odds ratios are STANDARDIZED. Each OR represents the change in odds for a 1 standard deviation
           increase in the predictor. This makes ORs comparable across variables but means the magnitude
           depends on how spread out each variable is in the data."
        } else {
          "These odds ratios are UNSTANDARDIZED. Each OR represents the change in odds for a 1-unit increase
           in the predictor (e.g., 1 percentage point for percentages, $1 for income). Magnitudes are not
           directly comparable across variables with different scales."
        }
    )
  })
  
  output$coef_plot <- renderPlotly({
    req(rv$fit)
    coefs <- tidy(rv$fit, effects = "fixed", conf.int = TRUE) %>%
      filter(term != "(Intercept)") %>%
      mutate(or = exp(estimate), or_lo = exp(conf.low), or_hi = exp(conf.high),
             var_name = gsub("scale\\(|\\)", "", term),
             label = sapply(var_name, get_label),
             Effect = ifelse(or > 1.5, "Strong positive",
                             ifelse(or > 1.1, "Moderate positive",
                                    ifelse(or < 0.67, "Strong negative",
                                           ifelse(or < 0.9, "Moderate negative", "Weak/null")))))
    
    if (nrow(coefs) == 0) return(plotly_empty() %>% layout(title = "No fixed effects"))
    
    p <- ggplot(coefs, aes(or, reorder(label, or),
                           text = paste0(label, "\nOR: ", round(or, 2),
                                         "\n95% CI: ", round(or_lo, 2), "-", round(or_hi, 2)))) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = or_lo, xmax = or_hi), height = 0.2, color = "gray40") +
      geom_point(aes(color = Effect), size = 4) +
      scale_color_manual(values = c("Strong positive" = "#d73027", "Moderate positive" = "#fc8d59",
                                    "Weak/null" = "#999999",
                                    "Moderate negative" = "#91bfdb", "Strong negative" = "#4575b4")) +
      labs(title = "Odds Ratios (with 95% CI)", x = "Odds Ratio", y = NULL) +
      theme_minimal() + theme(legend.position = "bottom")
    
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "h", y = -0.2))
  })
  
  output$roc_plot <- renderPlot({
    req(rv$roc_obj, rv$youden)
    roc_df <- data.frame(fpr = 1 - rv$roc_obj$specificities, tpr = rv$roc_obj$sensitivities)
    youden_pt <- data.frame(fpr = 1 - rv$youden$spec, tpr = rv$youden$sens)
    
    ggplot(roc_df, aes(fpr, tpr)) +
      geom_line(color = "#2c7fb8", linewidth = 1.1) +
      geom_abline(slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(data = youden_pt, color = "#d73027", size = 4) +
      annotate("text", x = 0.7, y = 0.2, label = paste("AUC =", round(auc(rv$roc_obj), 3)),
               size = 4, fontface = "bold", color = "#2c7fb8") +
      labs(title = "ROC Curve", x = "False Positive Rate", y = "True Positive Rate",
           caption = "Red dot = optimal threshold (Youden)") +
      theme_minimal(base_size = 10) + coord_equal()
  })
  
  output$cal_plot_decile <- renderPlot({
    req(rv$test_preds)
    cal <- rv$test_preds %>%
      filter(!is.na(pred_prob)) %>%
      mutate(decile = ntile(pred_prob, 10)) %>%
      group_by(decile) %>%
      summarize(pred = mean(pred_prob), obs = mean(outcome), n = n(), .groups = "drop")
    
    if (nrow(cal) == 0) return(NULL)
    
    ggplot(cal, aes(pred, obs)) +
      geom_abline(slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(size = 3, color = "#d73027") +
      geom_line(color = "#d73027", alpha = 0.5) +
      geom_text(aes(label = decile), vjust = -1.2, size = 3) +
      labs(title = "Calibration by Decile", x = "Mean Predicted", y = "Mean Observed",
           caption = "Numbers = risk decile") +
      theme_minimal(base_size = 10) +
      coord_cartesian(xlim = c(0, max(cal$pred) * 1.2), ylim = c(0, max(cal$obs) * 1.2))
  })
  
  output$cal_plot_full <- renderPlot({
    req(rv$test_preds)
    
    d10_cutoff <- quantile(rv$test_preds$pred_prob, 0.9, na.rm = TRUE)
    
    cal_full <- rv$test_preds %>%
      filter(!is.na(pred_prob)) %>%
      mutate(prob_bin = cut(pred_prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE),
             is_top_decile = pred_prob >= d10_cutoff) %>%
      group_by(prob_bin) %>%
      summarize(pred = mean(pred_prob), obs = mean(outcome), n = n(),
                all_top = all(is_top_decile), .groups = "drop") %>%
      filter(n >= 5)
    
    if (nrow(cal_full) == 0) return(NULL)
    
    ggplot(cal_full, aes(pred, obs, color = all_top)) +
      geom_abline(slope = 1, linetype = "dashed", color = "gray50") +
      geom_point(size = 4) +
      geom_line(aes(group = 1), color = "gray60", alpha = 0.5) +
      scale_color_manual(values = c("FALSE" = "#2c7fb8", "TRUE" = "#d73027"),
                         labels = c("Other deciles", "10th risk decile"), name = NULL) +
      labs(title = "Calibration (Full Range)", x = "Predicted Probability", y = "Observed Rate") +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  })
  
  output$cor_table <- renderDT({
    if (is.null(rv$cor_matrix)) return(NULL)
    
    cor_df <- as.data.frame(rv$cor_matrix)
    cor_df$Variable <- sapply(rownames(cor_df), get_label)
    cor_df <- cor_df %>% select(Variable, everything())
    names(cor_df)[-1] <- sapply(names(cor_df)[-1], get_label)
    
    datatable(cor_df, options = list(dom = 't', pageLength = 10, scrollX = TRUE), rownames = FALSE) %>%
      formatRound(columns = 2:ncol(cor_df), digits = 2) %>%
      formatStyle(columns = 2:ncol(cor_df),
                  backgroundColor = styleInterval(
                    c(-0.7, -0.5, 0.5, 0.7, 0.99),
                    c('#f8d7da', '#fff3cd', 'white', '#fff3cd', '#f8d7da', 'white')))
  })
  
  output$vif_table <- renderDT({
    if (is.null(rv$cor_matrix) || ncol(rv$cor_matrix) < 2) return(NULL)
    
    vif_vals <- tryCatch({ diag(solve(rv$cor_matrix)) }, error = function(e) NULL)
    if (is.null(vif_vals)) return(NULL)
    
    vif_df <- data.frame(Variable = sapply(names(vif_vals), get_label),
                         VIF = round(vif_vals, 2), stringsAsFactors = FALSE) %>%
      arrange(desc(VIF))
    
    datatable(vif_df, options = list(dom = 't', pageLength = 10), rownames = FALSE) %>%
      formatStyle("VIF", backgroundColor = styleInterval(c(5, 10), c("white", "#fff3cd", "#f8d7da")))
  })
  
  output$model_summary <- renderPrint({
    if (is.null(rv$fit)) return(cat("Run model first"))
    summary(rv$fit)
  })
  
  # ===== PREDICTION MAP =====
  
  output$pred_map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      setView(lng = -99.5, lat = 41.5, zoom = 7)
  })
  
  observe({
    
    # Select prediction dataset
    preds <- if (input$pred_type == "test") {
      rv$test_preds
    } else {
      rv$full_preds
    }
    
    if (is.null(preds) || nrow(preds) == 0) return()
    
    # Select geography level
    geo_col <- if (input$pred_geo == "tract") "tract_geoid" else "bg_geoid"
    geo_bounds <- if (input$pred_geo == "tract") ne_tracts else ne_bgs
    
    if (is.null(geo_bounds)) return()
    if (!geo_col %in% names(preds)) return()
    
    # Aggregate predicted risk
    geo_risk <- preds %>%
      dplyr::filter(!is.na(pred_prob), !is.na(.data[[geo_col]])) %>%
      dplyr::group_by(geo_id = .data[[geo_col]]) %>%
      dplyr::summarize(
        risk = mean(pred_prob),
        n = dplyr::n(),
        .groups = "drop"
      )
    
    if (nrow(geo_risk) == 0) return()
    
    # Join to geometry
    geo <- geo_bounds %>%
      dplyr::left_join(geo_risk, by = c("GEOID" = "geo_id"))
    
    max_r <- max(geo$risk, na.rm = TRUE)
    if (!is.finite(max_r) || max_r <= 0) max_r <- 0.3
    
    pal <- colorNumeric(
      palette = "YlOrRd",
      domain = c(0, max_r),
      na.color = "#cccccc"
    )
    
    geo$label <- ifelse(
      !is.na(geo$risk),
      paste0(
        "<strong>Predicted Risk:</strong> ",
        round(geo$risk * 100, 1), "%<br>",
        "<strong>Addresses:</strong> ", geo$n
      ),
      "No prediction"
    )
    
    leafletProxy("pred_map") %>%
      clearShapes() %>%
      clearControls() %>%
      addPolygons(
        data = geo,
        fillColor = ~pal(risk),
        fillOpacity = 0.7,
        color = "white",
        weight = 0.5,
        label = lapply(geo$label, HTML)
      ) %>%
      addLegend(
        position = "bottomright",
        pal = pal,
        values = c(0, max_r),
        title = "Predicted Risk",
        opacity = 0.7,
        labFormat = labelFormat(
          transform = function(x) round(x * 100),
          suffix = "%"
        )
      )
  })
}
  # Run app
  shinyApp(ui, server)
  