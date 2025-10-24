# install.packages(c("shiny","bslib","plotly","jsonlite"))
library(shiny)
library(bslib)
library(plotly)
library(jsonlite)

options(shiny.fullstacktrace = TRUE)

MPH_TO_KMH <- 1.609344

# ---------- Tefft (2013) ----------
bmi_from_imp <- function(height_in, weight_lb) {
    h_m <- height_in * 0.0254; w_kg <- weight_lb * 0.45359237
    w_kg / (h_m^2)
}
tefft_prob <- function(speed_mph, age_years, height_in, weight_lb, bmi = NULL, truck = 0) {
    BMI <- if (is.null(bmi)) bmi_from_imp(height_in, weight_lb) else bmi
    BMI_above25 <- pmax(BMI - 25, 0)
    w_term <- (weight_lb / 100)^(-0.5)
    z <- 196.8 +
        0.202  * speed_mph +
        0.00712 * (age_years / 10)^3 -
        1.174  * height_in -
        90.5   * w_term -
        2.247  * BMI +
        1.557  * BMI_above25 +
        0.543  * truck
    plogis(z)
}
tefft_prob_guarded <- function(speed_mph, age_years, height_in, weight_lb, bmi = NULL, truck = 0,
                               fallback = "speed_only") {
    if (age_years < 15) {
        if (fallback == "speed_only") {
            alpha <- -5.406504111506091
            beta  <-  0.1331651258991648
            return(plogis(alpha + beta * speed_mph))
        } else stop("Tefft not valid for age < 15")
    }
    tefft_prob(speed_mph, age_years, height_in, weight_lb, bmi, truck)
}

# ---------- Mueller/Monfort (2025) ----------
MU_MEANS <- list(speed_kmh = 42.79, hle_cm = 85.25, sex_male = 0.59, age_years = 42.42)
MU_COEFS <- list(
    mais2p = list(b0= 0.686, bv=0.083, bh=0.034, bs= 0.047, ba=0.023, bvh=0.001),
    mais3p = list(b0=-0.460, bv=0.091, bh=0.039, bs=-0.119, ba=0.029, bvh=0.002),
    fatal  = list(b0=-3.280, bv=0.123, bh=0.047, bs= 0.319, ba=0.030, bvh=0.000)
)
mueller_prob <- function(speed_kmh, hle_cm, sex_male, age_years, severity = "fatal") {
    c <- MU_COEFS[[severity]]; if (is.null(c)) stop("severity must be fatal|mais3p|mais2p")
    X2 <- hle_cm - MU_MEANS$hle_cm
    X3 <- sex_male - MU_MEANS$sex_male
    X4 <- age_years - MU_MEANS$age_years
    sapply(speed_kmh, function(v) {
        X1 <- v - MU_MEANS$speed_kmh
        eta <- c$b0 + c$bv*X1 + c$bh*X2 + c$bs*X3 + c$ba*X4 + c$bvh*(X1*X2)
        plogis(eta)
    })
}

# ---------- Vehicle -> HLE ----------
vehicle_hle_guess <- function(vehicle_type, suv_as = "pickup", van_as = "pickup") {
    s <- tolower(trimws(as.character(vehicle_type)))
    car      <- c("car","sedan","hatchback","coupe","wagon","passenger car","saloon")
    pickup   <- c("pickup","pickup truck","truck","pick-up","ute")
    suv      <- c("suv","sport utility vehicle","crossover","cuv")
    van      <- c("van","minivan","mpv","people carrier")
    if (s %in% car)    return(75)
    if (s %in% pickup) return(109)
    if (s %in% suv)    return(if (suv_as == "pickup") 109 else 75)
    if (s %in% van)    return(if (van_as == "pickup") 109 else 75)
    stop(sprintf("Unknown vehicle_type '%s'", vehicle_type))
}

# ---------- helpers ----------
linspace <- function(min, max, step) seq(min, max, by = step)
crossing_x <- function(xs, ys, yq) {
    if (yq < min(ys) || yq > max(ys)) return(NULL)
    idx <- which(ys >= yq)[1]
    if (is.na(idx) || idx <= 1) return(NULL)
    x0 <- xs[idx-1]; x1 <- xs[idx]; y0 <- ys[idx-1]; y1 <- ys[idx]
    if (y1 == y0) return(x0)
    x0 + (yq - y0) * (x1 - x0) / (y1 - y0)
}
slope_angle_deg <- function(xs, ys, xq) {
    idx <- which(xs >= xq)[1]; if (is.na(idx)) idx <- length(xs)
    if (idx < 2) idx <- 2
    dx <- xs[idx] - xs[idx-1]; dy <- ys[idx] - ys[idx-1]
    atan2(dy, dx) * 180 / pi
}

# ---------- UI ----------
ui <- page_fluid(
    theme = bs_theme(bootswatch = "darkly"),
    tags$head(tags$title("Pedestrian Risk — Tefft (2013) & Mueller (2025)")),
    layout_sidebar(
        sidebar = sidebar(
            h4("Plot"),
            selectInput("model", "Model", c("Tefft (2013)"="tefft", "Monfort & Mueller (2025)"="mueller", "Both"="both"), "both"),
            selectInput("units", "X-axis units", c("mph","kmh"), "mph"),
            fluidRow(
                column(4, numericInput("speedMin","Min", 5, step=1)),
                column(4, numericInput("speedMax","Max", 65, step=1)),
                column(4, numericInput("speedStep","Step", 1, step=1))
            ),
            textInput("title", "Title (optional)", ""),
            hr(),
            h4("Percent markers"),
            textInput("pctMarkers", "Levels (%)", "10,25,50,75"),
            checkboxInput("showLabels","Show labels", TRUE),
            selectInput("rotMode","Label rotation", c("Auto (follow curve)"="auto","None"="none","Fixed angle"="fixed"), "auto"),
            conditionalPanel("input.rotMode == 'fixed'",
                             fluidRow(
                                 column(6, numericInput("fixedAngle","Fixed angle (deg)", 30)),
                                 column(6, textInput("labelOffset","Label offset (px)", "10,6"))
                             )
            ),
            hr(),
            h4("Color pairing"),
            checkboxInput("matchColors","Pair colors across models", TRUE),
            selectInput("matchBy","Match by", c("Profile ID"="id","Position"="position"), "id"),
            textInput("idKey", "ID key (field name)", "profile_id"),
            checkboxInput("includeId","Include ID in legend", FALSE),
            hr(),
            h4("Mueller defaults"),
            selectInput("muDefaultSeverity","Default severity", c("fatal","mais3p","mais2p"), "fatal"),
            selectInput("suvAs","SUVs treated as", c("pickup","car"), "pickup"),
            hr(),
            h4("Profiles"),
            actionButton("addTefft","+ Add Tefft profile"),
            actionButton("addMueller","+ Add Mueller profile"),
            br(), br(),
            actionButton("clearTefft","Clear Tefft"),
            actionButton("clearMueller","Clear Mueller"),
            br(), br(),
            actionButton("render","Render plot", class = "btn-primary")
        ),
        main = div(
            plotlyOutput("plot", height = "640px"),
            hr(),
            h4("Saved inputs (live view)"),
            fluidRow(
                column(6, h5("Tefft (table)"), tableOutput("tblTefft"),
                       h6("Tefft (JSON)"), verbatimTextOutput("jsonTefft", placeholder = TRUE)),
                column(6, h5("Mueller (table)"), tableOutput("tblMueller"),
                       h6("Mueller (JSON)"), verbatimTextOutput("jsonMueller", placeholder = TRUE))
            )
        )
    )
)

# ---------- SERVER ----------
server <- function(input, output, session) {
    
    rv <- reactiveValues(
        tefft = data.frame(
            profile_id=c(1,2),
            age_years=c(35,35),
            height_in=c(72,66),
            weight_lb=c(160,120),
            vehicle_is_truck=c(0,1),
            stringsAsFactors = FALSE
        ),
        mueller = data.frame(
            profile_id=c(2,1),
            age_years=c(35,35),
            vehicle_type=c("car","pickup"),
            hle_cm=c(NA,NA),
            sex_male=c(0,1),
            severity=c("fatal","fatal"),
            stringsAsFactors = FALSE
        )
    )
    
    observeEvent(input$addTefft, {
        showModal(modalDialog(
            title = "Add Tefft profile",
            fluidRow(
                column(6, numericInput("t_pid","Profile ID (optional)", NA)),
                column(6, numericInput("t_age","Age (years)", 35, min = 1))
            ),
            fluidRow(
                column(6, numericInput("t_h","Height (in)", 72, min = 30, max = 90)),
                column(6, numericInput("t_w","Weight (lb)", 160, min = 20, max = 600))
            ),
            selectInput("t_truck","Vehicle", c("Car"=0,"Light truck/SUV"=1), 0),
            footer = tagList(modalButton("Cancel"), actionButton("t_add", "Add", class = "btn-primary")),
            easyClose = TRUE
        ))
    })
    observeEvent(input$t_add, {
        removeModal()
        rv$tefft <- rbind(rv$tefft, data.frame(
            profile_id = if (is.na(input$t_pid)) NA_integer_ else as.integer(input$t_pid),
            age_years = input$t_age, height_in = input$t_h, weight_lb = input$t_w,
            vehicle_is_truck = as.integer(input$t_truck),
            stringsAsFactors = FALSE
        ))
        showNotification("Tefft profile added", type = "message")
    })
    
    observeEvent(input$addMueller, {
        showModal(modalDialog(
            title = "Add Mueller profile",
            fluidRow(
                column(6, numericInput("m_pid","Profile ID (optional)", NA)),
                column(6, numericInput("m_age","Age (years)", 35, min = 1))
            ),
            fluidRow(
                column(6, numericInput("m_hle","HLE (cm)", NA)),
                column(6, textInput("m_vtype","Vehicle type (if HLE empty)", "car"))
            ),
            fluidRow(
                column(6, selectInput("m_sev","Severity", c("fatal","mais3p","mais2p"), selected = input$muDefaultSeverity)),
                column(6, selectInput("m_sex","Sex", c("Female"=0,"Male"=1), 1))
            ),
            footer = tagList(modalButton("Cancel"), actionButton("m_add", "Add", class = "btn-primary")),
            easyClose = TRUE
        ))
    })
    observeEvent(input$m_add, {
        removeModal()
        rv$mueller <- rbind(rv$mueller, data.frame(
            profile_id = if (is.na(input$m_pid)) NA_integer_ else as.integer(input$m_pid),
            age_years = input$m_age,
            vehicle_type = if (is.na(input$m_hle) || is.null(input$m_hle) || input$m_hle=="") input$m_vtype else NA_character_,
            hle_cm = if (!is.na(input$m_hle) && input$m_hle!="") as.numeric(input$m_hle) else NA_real_,
            sex_male = as.integer(input$m_sex),
            severity = input$m_sev,
            stringsAsFactors = FALSE
        ))
        showNotification("Mueller profile added", type = "message")
    })
    
    observeEvent(input$clearTefft,  { rv$tefft  <- rv$tefft[0,];  showNotification("Cleared Tefft profiles", type="message") })
    observeEvent(input$clearMueller,{ rv$mueller<- rv$mueller[0,]; showNotification("Cleared Mueller profiles", type="message") })
    
    output$tblTefft   <- renderTable(rv$tefft,   bordered = TRUE, striped = TRUE, hover = TRUE, spacing = "xs")
    output$tblMueller <- renderTable(rv$mueller, bordered = TRUE, striped = TRUE, hover = TRUE, spacing = "xs")
    output$jsonTefft   <- renderText(toJSON(rv$tefft,   pretty = TRUE, na = "null", dataframe = "rows"))
    output$jsonMueller <- renderText(toJSON(rv$mueller, pretty = TRUE, na = "null", dataframe = "rows"))
    
    make_plot <- eventReactive(input$render, {
        tryCatch({
            model <- input$model
            units <- input$units
            xPlot <- linspace(input$speedMin, input$speedMax, input$speedStep)
            x_mph <- if (units == "mph") xPlot else xPlot / MPH_TO_KMH
            x_kmh <- if (units == "kmh") xPlot else xPlot * MPH_TO_KMH
            
            pct <- suppressWarnings(as.numeric(strsplit(gsub("\\s+","", input$pctMarkers), ",")[[1]]))
            pct <- pct[is.finite(pct) & pct >= 0 & pct <= 100] / 100
            if (length(pct) == 0) pct <- c(0.1,0.25,0.5,0.75)
            
            rotMode <- input$rotMode
            fixedAngle <- if (rotMode == "fixed") input$fixedAngle else 0
            offs <- if (!is.null(input$labelOffset) && grepl(",", input$labelOffset)) strsplit(input$labelOffset,",")[[1]] else c("10","6")
            xshift <- suppressWarnings(as.numeric(offs[1])); if (!is.finite(xshift)) xshift <- 10
            yshift <- suppressWarnings(as.numeric(offs[2])); if (!is.finite(yshift)) yshift <- 6
            
            doTefft  <- model %in% c("tefft","both")
            doMueller<- model %in% c("mueller","both")
            
            validate(
                need(!(doTefft  && nrow(rv$tefft)==0),  "Add at least one Tefft profile or switch model."),
                need(!(doMueller&& nrow(rv$mueller)==0),"Add at least one Mueller profile or switch model.")
            )
            
            colorway <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd",
                          "#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf")
            color_map <- new.env(parent = emptyenv())
            color_for_id <- function(pid) {
                if (is.na(pid) || is.null(pid)) return(NA_character_)
                key <- paste0("id_", pid)
                if (!exists(key, envir = color_map, inherits = FALSE)) {
                    idx <- length(ls(envir = color_map, all.names = TRUE)) %% length(colorway) + 1
                    assign(key, colorway[idx], envir = color_map)
                }
                get(key, envir = color_map, inherits = FALSE)
            }
            
            plt <- plot_ly()
            ann <- list()
            tefft_colors_by_idx <- character(0)
            
            # ----- TEFFT (solid) -----
            if (doTefft && nrow(rv$tefft) > 0) {
                for (i in seq_len(nrow(rv$tefft))) {
                    p <- rv$tefft[i,]
                    y <- tefft_prob_guarded(x_mph, p$age_years, p$height_in, p$weight_lb, NULL, p$vehicle_is_truck, "speed_only")
                    col <- NA
                    if (isTRUE(input$matchColors)) {
                        if (input$matchBy == "id" && !is.na(p[[input$idKey]])) col <- color_for_id(p[[input$idKey]])
                        if (input$matchBy == "position") col <- colorway[(i-1) %% length(colorway) + 1]
                    }
                    tefft_colors_by_idx[i] <- ifelse(is.na(col), "", col)
                    veh <- if (p$vehicle_is_truck==1) "Light truck/SUV" else "Car"
                    pid <- p[[input$idKey]]
                    idtxt <- if (isTRUE(input$includeId) && !is.na(pid)) paste0("[",pid,"] ") else ""
                    lbl <- sprintf("%sTefft — %dy, %d\" %d lb, %s",
                                   idtxt, round(p$age_years), round(p$height_in), round(p$weight_lb), veh)
                    plt <- add_lines(
                        plt, x = xPlot, y = y, name = lbl,
                        line = list(color = ifelse(is.na(col), NULL, col), width = 3, dash = "solid"),
                        hovertemplate = sprintf("<b>%%{y:.2%%}</b> at %%{x:.1f} %s<extra></extra>", units)
                    )
                    for (q in pct) {
                        xq <- crossing_x(xPlot, y, q); if (is.null(xq)) next
                        plt <- add_markers(plt, x = xq, y = q, showlegend = FALSE,
                                           marker = list(size = 9, color = ifelse(is.na(col), NULL, col),
                                                         line = list(color = "#ffffff", width = 1)))
                        if (isTRUE(input$showLabels)) {
                            angle <- if (rotMode == "auto") slope_angle_deg(xPlot, y, xq) else if (rotMode=="fixed") fixedAngle else 0
                            ann[[length(ann)+1]] <- list(
                                x = xq, y = q, xref = "x", yref = "y",
                                text = sprintf("%d%% @ %.1f %s", round(q*100), xq, units),
                                showarrow = FALSE, xanchor = "left", yanchor = "bottom",
                                xshift = xshift, yshift = yshift, font = list(size = 11, color = "#e5e5e5"),
                                textangle = angle, bgcolor = "rgba(0,0,0,0)"
                            )
                        }
                    }
                }
            }
            
            # ----- MUELLER (dotted) -----
            if (doMueller && nrow(rv$mueller) > 0) {
                for (j in seq_len(nrow(rv$mueller))) {
                    p <- rv$mueller[j,]
                    hle <- if (!is.na(p$hle_cm)) p$hle_cm else vehicle_hle_guess(p$vehicle_type, suv_as = input$suvAs, van_as = "pickup")
                    sev <- if (nzchar(p$severity)) p$severity else input$muDefaultSeverity
                    y <- mueller_prob(x_kmh, hle, as.numeric(p$sex_male), p$age_years, sev)
                    col <- NA
                    if (isTRUE(input$matchColors)) {
                        if (input$matchBy == "id" && !is.na(p[[input$idKey]])) {
                            col <- color_for_id(p[[input$idKey]])
                        } else if (input$matchBy == "position" && j <= length(tefft_colors_by_idx)) {
                            col <- tefft_colors_by_idx[j]; if (identical(col,"")) col <- NA
                        }
                    }
                    sexlab <- if (p$sex_male==1) "Male" else "Female"
                    pid <- p[[input$idKey]]
                    idtxt <- if (isTRUE(input$includeId) && !is.na(pid)) paste0("[",pid,"] ") else ""
                    sev_map <- c(fatal = "Fatal", mais3p = "MAIS 3+F", mais2p = "MAIS 2+F")
                    lbl <- sprintf("%sMueller — %dy, HLE %d cm, %s, %s",
                                   idtxt, round(p$age_years), round(hle), sexlab, sev_map[[sev]])
                    plt <- add_lines(
                        plt, x = xPlot, y = y, name = lbl,
                        line = list(color = ifelse(is.na(col), NULL, col), width = 3, dash = "dot"),
                        hovertemplate = sprintf("<b>%%{y:.2%%}</b> at %%{x:.1f} %s<extra></extra>", units)
                    )
                    for (q in pct) {
                        xq <- crossing_x(xPlot, y, q); if (is.null(xq)) next
                        plt <- add_markers(plt, x = xq, y = q, showlegend = FALSE,
                                           marker = list(size = 9, color = ifelse(is.na(col), NULL, col),
                                                         line = list(color = "#ffffff", width = 1)))
                        if (isTRUE(input$showLabels)) {
                            angle <- if (rotMode == "auto") slope_angle_deg(xPlot, y, xq) else if (rotMode=="fixed") fixedAngle else 0
                            ann[[length(ann)+1]] <- list(
                                x = xq, y = q, xref = "x", yref = "y",
                                text = sprintf("%d%% @ %.1f %s", round(q*100), xq, units),
                                showarrow = FALSE, xanchor = "left", yanchor = "bottom",
                                xshift = xshift, yshift = yshift, font = list(size = 11, color = "#e5e5e5"),
                                textangle = angle, bgcolor = "rgba(0,0,0,0)"
                            )
                        }
                    }
                }
            }
            
            ttl <- if (nzchar(input$title)) input$title else
                if (input$model == "both") "Pedestrian risk vs. speed — Tefft (2013) & Monfort/Mueller (2025)"
            else if (input$model == "tefft") "Pedestrian death risk vs. speed — Tefft (2013)"
            else "Pedestrian injury risk vs. speed — Monfort & Mueller (2025)"
            
            plt <- layout(
                plt,
                template = "plotly_dark",
                paper_bgcolor = "#111111", plot_bgcolor = "#111111",
                title = list(text = ttl, x = 0.02, xanchor = "left"),
                xaxis = list(title = paste0("Impact speed (", input$units, ")"), gridcolor="#333333"),
                yaxis = list(title = "Probability (death/injury)", gridcolor="#333333", range = c(0,1)),
                legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0, bgcolor="rgba(0,0,0,0)"),
                margin = list(l=64, r=24, t=64, b=64),
                annotations = ann
            )
            plt
        }, error = function(e) {
            validate(need(FALSE, paste("Plot error:", conditionMessage(e))))
            NULL
        })
    })
    
    output$plot <- renderPlotly(make_plot())
}

shinyApp(ui, server)


# install.packages('rsconnect')
# 
# rsconnect::setAccountInfo(name='brianalmdale',
#                           token='FF800070A9B10AE8A34C9E445EA1373E',
#                           secret='K1viwK/gRt+Aq7ROap/xr2TGsqOlMa7c+vJImGKg')
# 
# library(rsconnect)
# rsconnect::deployApp('/Users/balmdale/code/ped_fatality_risk')
