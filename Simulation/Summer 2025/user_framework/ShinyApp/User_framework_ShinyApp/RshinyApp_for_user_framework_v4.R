library(shiny)

source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/Summer 2025/user_framework/my_functions_v2.R")
ui <- fluidPage(
  titlePanel("User Framework for Statistical Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("test_type", "Test Type:",
                  choices = c("Independent t-test" = "ttest",
                              "Regression" = "regression",
                              "ANOVA" = "anova")),
      selectInput(
        "metric",
        "Analysis Metric:",
        choices = c("Power", "Type I Error")
      ),
      
      fileInput("gen_data_file", "Upload Data Generation Function"),
      fileInput("get_params_file", "Upload Parameter Function"),
      fileInput("norm_obj_file", "Upload Normality Object Function"),
      fileInput("test1_file", "Upload Test 1 Function"),
      fileInput("test2_file", "Upload Test 2 Function"),
      
      selectInput("norm_test", "Normality Test:",
                  choices = c("Shapiro-Wilk" = "SW",
                              "Shapiro-Francia" = "SF",
                              "Kolmogorov-Smirnov" = "KS",
                              "Jarque-Bera" = "JB",
                              "D'Agostino Pearson" = "DAP",
                              "Anderson-Darling" = "AD",
                              "Cramer-Von-Mises" = "CVM")),
      numericInput("alpha", "Significance Level (α):", 0.05, min = 0.01, max = 0.2, step = 0.01),
      sliderInput("sample_sizes", "Sample Sizes Range:", min = 5, max = 200, value = c(10, 50), step = 5),
      numericInput("n_sim", "Number of Simulations:", value =  1000, min = 100, max = 10000),
      actionButton("run", "Run Simulation", class = "btn-primary"),
      
      tags$hr(),
      tags$p("Each uploaded script must define exactly one function with these names:"),
      tags$ul(
        tags$li("gen_data: Data generation function"),
        tags$li("get_parameters: Parameter setup function"),
        tags$li("fn_get_norm_obj: Normality object extractor"),
        tags$li("fn_for_ds_test_1: Test 1 function"),
        tags$li("fn_for_ds_test_2: Test 2 function")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Power Curves", plotOutput("power_plot")),
        tabPanel("AUC Results", tableOutput("auc_table")),
        tabPanel("Documentation",
                 tags$h3("User Framework for Statistical Analysis"),
                 tags$p("This app simulates power/Type I error for different statistical test approaches."),
                 tags$h4("Inputs:"),
                 tags$ul(
                   tags$li("Downstream test: t-test, ANOVA, or regression, etc"),
                   tags$li("Normality test: SW, KS, JB, etc."),
                   tags$li("Significance levels,α"),
                   tags$li("Sample size range"),
                   tags$li("Number of simulations(N)")
                 ),
                 tags$h4("Outputs:"),
                 tags$ul(
                   tags$li("Power/Type I error rate curves vs. sample size"),
                   tags$li("Area Under the Curve (AUC) comparison")
                 ),
                 tags$h4("Function Requirements:"),
                 tags$ul(
                   tags$li("gen_data(n, ...): Returns dataset for analysis"),
                   tags$li("get_parameters(n): Returns list of parameters for gen_data"),
                   tags$li("fn_get_norm_obj(data): Returns residuals or groups for normality test"),
                   tags$li("fn_for_ds_test_1(data): Returns p-value for parametric test"),
                   tags$li("fn_for_ds_test_2(data): Returns p-value for nonparametric test")
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  gen_data_override <- function(n) {
    
    # grab the user’s parameter
    user_args <- get_parameters(n)
    
    # otherwise grab the default arguments of gen_data()
    default_args <- formals(gen_data)
    
    # merge them 
    final_args <- modifyList(as.list(default_args), user_args)
    
    # call gen_data() 
    do.call(gen_data, final_args)
  }
  
  results <- eventReactive(input$run, {
    # range of sample sizes
    sizes <- seq(input$sample_sizes[1], input$sample_sizes[2], by = 5)
    
    # Create environment for user functions
    user_env <- new.env()
    
    # Source files with error handling
    try_source <- function(file, env) {
      if (is.null(file)) stop("No file uploaded")
      source(file$datapath, local = env)
    }
    
    tryCatch({
      # Source all user functions
      if (!is.null(input$gen_data_file)) try_source(input$gen_data_file, user_env)
      if (!is.null(input$get_params_file)) try_source(input$get_params_file, user_env)
      if (!is.null(input$norm_obj_file)) try_source(input$norm_obj_file, user_env)
      if (!is.null(input$test1_file)) try_source(input$test1_file, user_env)
      if (!is.null(input$test2_file)) try_source(input$test2_file, user_env)
    }, error = function(e) {
      showNotification(paste("Error loading functions:", e$message), type = "error")
      return(NULL)
    })
    
    # Verify required functions exist
    required_funcs <- c("gen_data", "get_parameters", "fn_get_norm_obj", 
                        "fn_for_ds_test_1", "fn_for_ds_test_2")
    
    missing_funcs <- setdiff(required_funcs, ls(envir = user_env))
    
    if (length(missing_funcs) > 0) {
      showNotification(paste("Missing functions:", paste(missing_funcs, collapse = ", ")), type = "error")
      return(NULL)
    }
    
    # Extract functions
    gen_data       <- user_env$gen_data
    get_parameters <- user_env$get_parameters
    fn_get_norm_obj <- user_env$fn_get_norm_obj
    fn_ds_test_1   <- user_env$fn_for_ds_test_1
    fn_ds_test_2   <- user_env$fn_for_ds_test_2
    
    # total work steps = 3 modes
    withProgress(message = "Running simulation...", value = 0, {
      # 1/3 per mode
      step <- 1/3
      
      # 1) Parametric
      incDetail <- function(pct) {
        # round to integer percent
        detail_text <- paste0(round(pct * 100), "% complete")
        setProgress(detail = detail_text)
      }
      
      # initial detail
      incDetail(0)
      
    # Run power analysis
    power_param <- calculate_power(
      sample_sizes       = sizes,
      n_sim              = input$n_sim,
      alpha              = input$alpha,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_get_norm_obj,
      fn_for_norm_test   = normality_test,
      fn_for_ds_test_1   = fn_ds_test_1,
      fn_for_ds_test_2   = fn_ds_test_2,
      test_method        = input$norm_test,
      mode               = "parametric"
    )
    
    incProgress(step); incDetail(step)
    
    power_nonpar <- calculate_power(
      sample_sizes       = sizes,
      n_sim              = input$n_sim,
      alpha              = input$alpha,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_get_norm_obj,
      fn_for_norm_test   = normality_test,
      fn_for_ds_test_1   = fn_ds_test_1,
      fn_for_ds_test_2   = fn_ds_test_2,
      test_method        = input$norm_test,
      mode               = "nonparametric"
    )
    
    incProgress(step); incDetail(2 * step)
    
    power_adapt <- calculate_power(
      sample_sizes       = sizes,
      n_sim              = input$n_sim,
      alpha              = input$alpha,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_get_norm_obj,
      fn_for_norm_test   = normality_test,
      fn_for_ds_test_1   = fn_ds_test_1,
      fn_for_ds_test_2   = fn_ds_test_2,
      test_method        = input$norm_test,
      mode               = "adaptive"
    )
    
    incProgress(step); incDetail(1)
    
    
    auc_df <- data.frame(
      Method = c("Test 1", "Test 2", "Test 3(Adaptive)"),
      AUC = c(
        compute_area(sizes, power_param),
        compute_area(sizes, power_nonpar),
        compute_area(sizes, power_adapt)
      )
    )
    
    list(
      sizes = sizes,
      power_param = power_param,
      power_nonpar = power_nonpar,
      power_adapt = power_adapt,
      auc = auc_df
    )
    })  # end withProgress
  })
  
  
  output$power_plot <- renderPlot({
    res <- results()
    if (is.null(res)) return(NULL)
    
    # isolate create a non-reactive scope
    metric    <- isolate(input$metric)
    test_type <- isolate(input$test_type)
    norm_test <- isolate(input$norm_test)
    
    # dynamically decide y-limits and labels
    is_power   <- input$metric == "Power"
    threshold  <- if (is_power) 0.8 else 0.05
    legend_pos <- if (is_power) "bottomright" else "topleft"
    y_limits <- if (is_power) c(0, 1) else c(0, 0.1)
    y_label <- if (is_power) "Power" else "Type I error"
    
    # unpack
    sizes     <- res$sizes
    p_param   <- res$power_param
    p_nonpar  <- res$power_nonpar
    p_adapt   <- res$power_adapt
    auc_vals  <- res$auc$AUC
    
    # main + subtitle
    main_title <- paste(input$metric, "vs Sample Size")
    sub_title  <- paste("Test Type:", input$test_type, 
                        "| Normality Test:", input$norm_test)
    
    # create the plot 
    plot(sizes, p_param, type = "b", pch = 19,
         col  = "red",
         ylim = y_limits,
         xlab = "Sample Size", ylab = y_label,
         main = main_title, sub = sub_title)
    
    # add the others
    lines(sizes, p_nonpar, type = "b", pch = 17, col = "blue")
    lines(sizes, p_adapt,  type = "b", pch = 15, col = "green")
    
    # reference line
    abline(h = threshold, lty = 2, col = "gray50")
    
    # legend with AUC
legend(legend_pos,
      legend = paste0(
      c("Test1(Parametric)","Test2(Nonparametric)","Test3(Adaptive)"),
      " (AUC=", formatC(res$auc$AUC, format = "f", digits = 4),")"
      ),
      title    = paste("AUC for", input$metric),
      col      = c("red","blue","green"),
      pch      = c(19,17,15),
      lty      = 1,
      bty      = "o",        
      box.lwd  = 1.5         
    )
  })
  
  output$auc_table <- renderTable({
    res <- results()
    if (is.null(res)) return(NULL)
    res$auc
  }, striped = TRUE, digits = 4)
}

shinyApp(ui, server)