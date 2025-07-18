library(shiny)

# Load helper source files
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/user_framework/ROC Curves/ROC_curves_for_ds_test_v1.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/user_framework/ShinyApp/general_functions/my_functions_v2.R")
compute_area <- compute_area
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/user_framework/ShinyApp")

# Modify plotting functions to support NULL filename (display in Shiny)
plot_power_error_tradeoff <- function(alpha_pretest, metrics, file_name = NULL) {
  if (!is.null(file_name)) pdf(file_name, width = 12, height = 10)
  par(mfrow = c(2, 2))
  plot(alpha_pretest, 
       metrics$EPL, 
       type = "l", 
       col = "blue", 
       lwd = 2,
       ylab = "Expected Power Loss (EPL)", 
       xlab = expression(alpha),
       main = "Power Loss (Normal)")
  
  plot(alpha_pretest, 
       metrics$EPG, 
       type = "l", 
       col = "red", 
       lwd = 2,
       ylab = "Expected Power Gain (EPG)", 
       xlab = expression(alpha),
       main = "Power Gain (LogNormal)")
  
  plot(alpha_pretest, 
       metrics$Expected_Inflation_normal, 
       type = "l", 
       col = "orange", 
       lwd = 2,
       ylab = "Type I Error Inflation", 
       xlab = expression(alpha),
       main = "Inflation (Normal)")
  
  plot(alpha_pretest, 
       metrics$Expected_Inflation_lognormal, 
       type = "l", 
       col = "green", 
       lwd = 2,
       ylab = "Type I Error Inflation", 
       xlab = expression(alpha),
       main = "Inflation (LogNormal)")
  if (!is.null(file_name)) dev.off()
}

plot_roc_like_curves <- function(metrics, file_name = NULL) {
  if (!is.null(file_name)) pdf(file_name, width = 12, height = 6)
  par(mfrow = c(1, 2))
  plot(metrics$EPL, 
       metrics$EPG, 
       type = "l", 
       col = "blue", 
       lwd = 2,
       xlab = "Power Loss (Normal)", 
       ylab = "Power Gain (LogNormal)",
       main = "ROC-like Curve: EPG vs. EPL")
  legend("bottomright", 
         legend = paste("EPG Gain =", round(metrics$EPG_lognormal, 3)), 
         title = "Point Estimate")
  
  plot(metrics$Expected_Inflation_normal, 
       metrics$Expected_Inflation_lognormal, 
       type = "l", 
       col = "red", 
       lwd = 2,
       xlab = "Type I Error Inflation (Normal)",
       ylab = "Type I Error Inflation (LogNormal)",
       main = "ROC-like Curve: Type I Inflation")
  if (!is.null(file_name)) dev.off()
}

# ROC Curves for Normality Tests
plot_ROC <- function(FPR, 
                     TPR, 
                     tests_to_plot = rownames(FPR), 
                     alpha = NULL, 
                     title = "ROC Curves") {
  
  colors <- rainbow(length(tests_to_plot))
  plot_chars <- 1:length(tests_to_plot)
  
  par(mar = c(5, 5, 4, 2))
  plot(0, 0, 
       type = "n", 
       xlim = c(0, 1), 
       ylim = c(0, 1),
       xlab = "False Positive Rate (FPR)", 
       ylab = "True Positive Rate (TPR)",
       main = title)
  for (i in seq_along(tests_to_plot)) {
    test <- tests_to_plot[i]
    lines(FPR[test, ], TPR[test, ], 
          col = colors[i], lwd = 2,
          lty = 1)
    points(FPR[test, ], TPR[test, ], 
           col = colors[i], 
           pch = plot_chars[i], 
           cex = 0.75)
  }
  abline(0, 1, lty = 2, col = "gray")
  legend("bottomright", 
         legend = tests_to_plot, 
         col = colors, 
         pch = plot_chars, 
         lty = 1,
         lwd = 2, 
         title = "Tests")
}

# UI
ui <- fluidPage(
  titlePanel("User Framework for Statistical Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("test_type", "Test Type:",
                  choices = c("Independent t-test" = "ttest",
                              "Regression" = "regression",
                              "ANOVA" = "anova")),
      selectInput("metric","Analysis Metric:", 
                  choices = c("Power", "Type I Error")),
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
      numericInput("alpha", "Significance Level (α):", value = 0.05, min = 0.01, max = 0.2, step = 0.01),
      sliderInput("sample_sizes", "Sample Sizes Range:", min = 5, max = 200, value = c(10, 50), step = 5),
      numericInput("n_sim", "Number of Simulations:", value = 1000, min = 10, max = 100000),
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
        tabPanel("Power/Type I error Curves", plotOutput("power_plot")),
        tabPanel("AUC Results", tableOutput("auc_table")),
        tabPanel("ROC Curve Analysis",
                 fluidRow(
                   column(
                     4,
                     checkboxGroupInput("roc_tests", "Select Normality Tests:",
                                        choices  = c("SW","KS","AD","DAP","SF","JB","CVM","SKEW"),
                                        selected = c("SW","KS","AD","DAP")),
                     
                     numericInput("roc_sample_size", "Sample Size for ROC:", 10, 5, 200),
                     numericInput("roc_Nsim",        "Number of Simulations:", 1000, 100, 1e6),
                     
                     sliderInput("alpha_range", HTML("&alpha;<sub>pretest</sub> range"),
                                 min = 0.001, max = 1, value = c(0.001, 1), step = 0.001),
                     numericInput("alpha_by", HTML("&alpha;<sub>pretest</sub> step"),
                                  value = 0.01, min = 0.001, max = 0.5, step = 0.001),
                     
                     sliderInput("sig_range", "sig_level range",
                                 min = 0.01,  max = 1, value = c(0.01, 1), step = 0.01),
                     numericInput("sig_by", "sig_level step",
                                  value = 0.01, min = 0.001, max = 0.5, step = 0.001),
                     
                     actionButton("run_roc", "Run ROC Analysis", class = "btn-success")
                   ),
                   column(8,
                          tabsetPanel(
                            tabPanel("ROC Curves", plotOutput("roc_curve_plot")),
                            tabPanel("Power/Error Tradeoff", plotOutput("power_error_plot")),
                            tabPanel("ROC-like Metrics", plotOutput("roc_like_metrics_plot"))
                          )
                   )
                 )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # --- Simulation Results ---
  results <- eventReactive(input$run, {
    sizes <- seq(input$sample_sizes[1], input$sample_sizes[2], by = 5)
    user_env <- new.env()
    
    try_source <- function(file, env) {
      if (is.null(file)) stop("No file uploaded")
      source(file$datapath, local = env)
    }
    
    tryCatch({
      if (!is.null(input$gen_data_file)) try_source(input$gen_data_file, user_env)
      if (!is.null(input$get_params_file)) try_source(input$get_params_file, user_env)
      if (!is.null(input$norm_obj_file)) try_source(input$norm_obj_file, user_env)
      if (!is.null(input$test1_file)) try_source(input$test1_file, user_env)
      if (!is.null(input$test2_file)) try_source(input$test2_file, user_env)
    }, error = function(e) {
      showNotification(paste("Error loading functions:", e$message), type = "error")
      return(NULL)
    })
    
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
    
    # unpack values
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
  
  
  # --- ROC Curve Analysis ---
  roc_results <- eventReactive(input$run_roc, {
    tests <- input$roc_tests
    Nsim <- input$roc_Nsim
    n <- input$roc_sample_size
    alpha_pretest <- seq(from = input$alpha_range[1], to   = input$alpha_range[2], by   = input$alpha_by
    )
    
    sig_level <- seq(from = input$sig_range[1], to   = input$sig_range[2], by   = input$sig_by
    )
    
    
    withProgress(message = "Running ROC Analysis...", value = 0, {
      roc_res <- ROC_curve_function(sample_size = n,
                                    alpha_pretest = alpha_pretest,
                                    tests = tests,
                                    Nsim = Nsim)
      incProgress(1/3)
      
      sim_data <- perform_analysis(sample_size = n,
                                   N = Nsim,
                                   distributions = c("Normal", "LogNormal"),
                                   effect_size = 0.75,
                                   alpha_pretest = alpha_pretest)
      incProgress(1/3)
      
      metrics <- compute_roc_metrics(sim_data$power_results, sim_data$error_results)
      incProgress(1/3)
      
      roc_data <- generate_roc_data(n, Nsim,
                                    c("Normal", "LogNormal"),
                                    effect_size = 0.75,
                                    sig_levels = sig_level)
      
      list(FPR = roc_res$FPR,
           TPR = roc_res$TPR,
           alpha = alpha_pretest,
           tests = tests,
           sim_data = sim_data,
           metrics = metrics,
           roc_data = roc_data)
    })
  })
  
  output$roc_curve_plot <- renderPlot({
    res <- roc_results()
    if (is.null(res)) return(NULL)
    plot_ROC(FPR = res$FPR, TPR = res$TPR, tests_to_plot = res$tests,
             alpha = res$alpha, title = "ROC Curves for Normality Tests")
  })
  
  output$power_error_plot <- renderPlot({
    res <- roc_results()
    if (is.null(res)) return(NULL)
    plot_power_error_tradeoff(alpha_pretest = res$alpha, metrics = res$metrics)
  })
  
  output$roc_like_metrics_plot <- renderPlot({
    res <- roc_results()
    if (is.null(res)) return(NULL)
    plot_roc_like_curves(metrics = res$metrics)
  })
}

shinyApp(ui, server)
