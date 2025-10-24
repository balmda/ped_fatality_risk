# install.packages(c("shiny","plotly"))
library(shiny)
library(plotly)

MPH_TO_KMH <- 1.609344

# ----- Tefft (2013) -----
bmi_from_imp <- function(height_in, weight_lb) {
    h_m <- height_in * 0.0254; w_kg <- weight_lb * 0.45359237; w_kg / (h_m^2)
}
tefft_prob <- function(speed_mph, age_years, height_in, weight_lb, bmi = NULL, truck = 0) {
    BMI <- if (is.null(bmi)) bmi_from_imp(height_in, weight_lb) else bmi
    BMI_above25 <- pmax(BMI - 25, 0)
    w_term <- (weight_lb / 100)^(-0.5)
    z <- 196.8 + 0.202*speed_mph + 0.00712*(age_years/10)^3 - 1.174*height_in - 90.5*w_term - 2.247*BMI + 1.557*BMI_above25 + 0.543*truck
    plogis(z)
}
tefft_prob_guarded <- function(speed_mph, age_years, height_in, weight_lb, bmi = NULL, truck = 0) {
    if (age_years < 15) {
        alpha <- -5.406504111506091; beta <- 0.1331651258991648
        return(plogis(alpha + beta * speed_mph))
    }
    tefft_prob(speed_mph, age_years, height_in, weight_lb, bmi, truck)
}

# ----- Mueller/Monfort (2025) -----
MU_MEANS <- list(speed_kmh = 42.79, hle_cm = 85.25, sex_male = 0.59, age_years = 42.42)
MU_COEFS <- list(
    fatal  = list(b0=-3.280, bv=0.123, bh=0.047, bs= 0.319, ba=0.030, bvh=0.000),
    mais3p = list(b0=-0.460, bv=0.091, bh=0.039, bs=-0.119, ba=0.029, bvh=0.002),
    mais2p = list(b0= 0.686, bv=0.083, bh=0.034, bs= 0.047, ba=0.023, bvh=0.001)
)
mueller_prob <- function(speed_kmh, hle_cm, sex_male, age_years, severity = "fatal") {
    c <- MU_COEFS[[severity]]; stopifnot(!is.null(c))
    X2 <- hle_cm - MU_MEANS$hle_cm
    X3 <- sex_male - MU_MEANS$sex_male
    X4 <- age_years - MU_MEANS$age_years
    sapply(speed_kmh, function(v) {
        X1 <- v - MU_MEANS$speed_kmh
        eta <- c$b0 + c$bv*X1 + c$bh*X2 + c$bs*X3 + c$ba*X4 + c$bvh*(X1*X2)
        plogis(eta)
    })
}
vehicle_hle_guess <- function(vehicle_type, suv_as = "pickup", van_as = "pickup") {
    s <- tolower(trimws(as.character(vehicle_type)))
    car <- c("car","sedan","hatchback","coupe","wagon","passenger car","saloon")
    pickup <- c("pickup","pickup truck","truck","pick-up","ute")
    suv <- c("suv","sport utility vehicle","crossover","cuv")
    van <- c("van","minivan","mpv","people carrier")
    if (s %in% car) return(75)
    if (s %in% pickup) return(109)
    if (s %in% suv) return(if (suv_as=="pickup") 109 else 75)
    if (s %in% van) return(if (van_as=="pickup") 109 else 75)
    stop(sprintf("Unknown vehicle_type '%s'", vehicle_type))
}

ui <- fluidPage(
    titlePanel("Pedestrian risk — Tefft (2013) & Mueller (2025)"),
    sidebarLayout(
        sidebarPanel(
            selectInput("model", "Model", c("Tefft"="tefft","Mueller"="mueller","Both"="both"), "both"),
            selectInput("units", "X-axis units", c("mph","kmh"), "mph"),
            numericInput("speedMin","Min speed", 5, step=1),
            numericInput("speedMax","Max speed", 65, step=1),
            numericInput("speedStep","Step", 1, step=1),
            tags$hr(),
            h4("Tefft (single)"),
            numericInput("t_age","Age (y)", 35, min=1),
            numericInput("t_h","Height (in)", 72, min=30, max=90),
            numericInput("t_w","Weight (lb)", 160, min=20, max=600),
            selectInput("t_truck","Vehicle", c("Car"=0,"Light truck/SUV"=1), 0),
            tags$hr(),
            h4("Mueller (single)"),
            numericInput("m_age","Age (y)", 35, min=1),
            selectInput("m_sex","Sex", c("Female"=0,"Male"=1), "1"),
            selectInput("m_sev","Severity", c("fatal","mais3p","mais2p"), "fatal"),
            textInput("m_vtype","Vehicle type (car|pickup|suv|van)", "car"),
            numericInput("m_hle","HLE (cm) — leave 0 to use type", 0, min=0, step=1),
            actionButton("render","Render plot", class="btn-primary")
        ),
        mainPanel(
            plotlyOutput("plot", height = "600px"),
            tags$hr(),
            strong("Status:"), textOutput("status", inline = TRUE)
        )
    )
)

server <- function(input, output, session) {
    output$status <- renderText(Sys.time())
    
    make_plot <- eventReactive(input$render, {
        tryCatch({
            units <- input$units
            x <- seq(input$speedMin, input$speedMax, by = input$speedStep)
            x_mph <- if (units=="mph") x else x / MPH_TO_KMH
            x_kmh <- if (units=="kmh") x else x * MPH_TO_KMH
            
            doTefft <- input$model %in% c("tefft","both")
            doMueller <- input$model %in% c("mueller","both")
            
            p <- plot_ly()
            if (doTefft) {
                y <- tefft_prob_guarded(x_mph, input$t_age, input$t_h, input$t_w, NULL, as.integer(input$t_truck))
                p <- add_lines(p, x = x, y = y, name = sprintf("Tefft — %dy, %d\" %d lb, %s",
                                                               round(input$t_age), round(input$t_h), round(input$t_w),
                                                               if (as.integer(input$t_truck)==1) "Light truck/SUV" else "Car"),
                               line = list(width=3))
            }
            if (doMueller) {
                hle <- if (isTRUE(input$m_hle > 0)) input$m_hle else vehicle_hle_guess(input$m_vtype)
                y2 <- mueller_prob(x_kmh, hle, as.numeric(input$m_sex), input$m_age, input$m_sev)
                p <- add_lines(p, x = x, y = y2, name = sprintf("Mueller — %dy, HLE %d cm, %s, %s",
                                                                round(input$m_age), round(hle),
                                                                if (as.integer(input$m_sex)==1) "Male" else "Female", input$m_sev),
                               line = list(width=3, dash="dot"))
            }
            
            title <- if (input$model=="both") "Tefft & Mueller"
            else if (input$model=="tefft") "Tefft (2013)"
            else "Mueller (2025)"
            layout(p,
                   title = title,
                   xaxis = list(title = paste0("Impact speed (", units, ")")),
                   yaxis = list(title = "Probability", range = c(0,1))
            )
        }, error = function(e) {
            validate(need(FALSE, paste("Plot error:", conditionMessage(e))))
            NULL
        })
    })
    
    # render once on app start too
    observeEvent(TRUE, { session$sendCustomMessage("dummy", NULL); }, once = TRUE)
    output$plot <- renderPlotly(make_plot())
}

shinyApp(ui, server)


