library(shiny)

# Load helper source files
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/user_framework/ShinyApp")

# -----------------------------------------------------------------
# UI (Designs the user interface layout)
# -----------------------------------------------------------------
ui <- fluidPage(
  titlePanel("User Framework for Statistical Analysis"),
  sidebarLayout(
    sidebarPanel(
      selectInput("test_type", "Test Type:",
                  choices = c("Independent t-test" = "ttest","Regression" = "regression","ANOVA" = "anova")),
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
      
      numericInput("alpha", "Significance Level (α):", 
                   value = 0.05, 
                   min = 0.01, 
                   max = 0.2, 
                   step = 0.01),
      sliderInput("sample_sizes", "Sample Sizes Range:", 
                  min = 5, 
                  max = 200, 
                  value = c(10, 50), 
                  step = 5),
      numericInput("n_sim", "Number of Simulations:", 
                   value = 1000, 
                   min = 10, 
                   max = 100000),
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
        tabPanel("Power/Type I error Curves", 
                 plotOutput("power_plot"),
                 uiOutput("power_caption")
        ),
        tabPanel("AUC Results", 
                 tableOutput("auc_table"),
                 uiOutput("auc_table_caption")
        ),
        tabPanel("ROC Curve Analysis",
                 fluidRow(
                   column(
                     4,
                     checkboxGroupInput("roc_tests", "Select Normality Tests:",
                                        choices  = c("SW","KS","AD","DAP","SF","JB","CVM","SKEW"),
                                        selected = c("SW","KS")),
                     
                     numericInput("roc_sample_size", "Sample Size for ROC:", 
                                  value = 10, 
                                  min = 5, 
                                  max = 200),
                     numericInput("roc_Nsim", "Number of Simulations:", 
                                  value = 1000, 
                                  min = 100, 
                                  max = 1e6),
                     numericInput("effect_size", "Effect Size:", 
                                  value = 0.5, 
                                  min =  0.0, 
                                  max = 1),
                     numericInput("test_alpha", "Fixed Significance Level (α):", 
                                  value = 0.05, 
                                  min =  0.0, 
                                  max = 0.1),
                     
                     sliderInput("alpha_range", HTML("&alpha;<sub>pretest</sub> range"),
                                 min = 0.001, max = 1, value = c(0.001, 1), step = 0.001),
                     numericInput("alpha_by", HTML("&alpha;<sub>pretest</sub> step"),
                                  value = 0.001, min = 0.001, max = 0.5, step = 0.001),
                     
                     sliderInput("sig_range", "sig_level range",
                                 min = 0.001,  max = 1, value = c(0.001, 1), step = 0.001),
                     numericInput("sig_by", "sig_level step",
                                  value = 0.001, min = 0.001, max = 0.5, step = 0.001),
                     
                     actionButton("run_roc", "Run ROC Analysis", class = "btn-success")
                   ),
                   column(8,
                          tabsetPanel(
                            tabPanel("ROC Curves", 
                                     plotOutput("roc_curve_plot"),
                                     uiOutput("roc_caption")
                            ),
                            tabPanel("Power/Error Tradeoff", 
                                     plotOutput("power_error_plot"),
                                     uiOutput("tradeoff_caption")
                            ),
                            tabPanel("ROC-like Metrics", 
                                     plotOutput("roc_like_metrics_plot"),
                                     uiOutput("roc_like_caption"))
                          )
                   )
                 )
        ),
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

# --------------------------------------------------------------------
# server : Put all R functions here
# --------------------------------------------------------------------
server <- function(input, output, session){
 
  # Ensure all required user-defined functions are loaded
  user_funcs <- reactive({
    req(input$gen_data_file, input$get_params_file, input$norm_obj_file, input$test1_file, input$test2_file)
    
    # Open a new environment for files
    user_env <- new.env()
    
    # Check for missing files and source them
    try_source <- function(file, env){
      if(is.null(file)) stop("No file uploaded")
      source(file$datapath, local = env)
    }
    
    # ensure each required function is properly loaded in user_env
    tryCatch({
      try_source(input$gen_data_file, user_env)
      try_source(input$get_params_file, user_env)
      try_source(input$norm_obj_file, user_env)
      try_source(input$test1_file, user_env)
      try_source(input$test2_file, user_env)
    },
    error = function(e){
      showNotification(paste("Error loading functions:", e$message), type = "error")
      return(NULL)
    })
    
    # Defines a list of function names that must be present.
    required_funcs <- c("gen_data", "get_parameters", "fn_get_norm_obj", "fn_for_ds_test_1", "fn_for_ds_test_2")
    missing_funcs <- setdiff(required_funcs, ls(envir = user_env))
    if(length(missing_funcs) > 0){
      showNotification(paste("Missing functions:", paste(missing_funcs, collapse = ", ")), type = "error")
      return(NULL)
    }
    # Return all functions
    list(
      gen_data = user_env$gen_data,
      get_parameters = user_env$get_parameters,
      fn_get_norm_obj = user_env$fn_get_norm_obj,
      fn_ds_test_1 = user_env$fn_for_ds_test_1,
      fn_ds_test_2 = user_env$fn_for_ds_test_2
    )
  })
  
  # -------------------------------------------------------------------
  #-------------- Define all other functions here ---------------------
  # -------------------------------------------------------------------
  #----------   function for computing auc ----------------------------
  compute_area <- function(sample_sizes, y) {
    if (is.null(y) || length(y) != length(sample_sizes)) return(NA_real_)
    sum(diff(sample_sizes) * (head(y, -1) + tail(y, -1)) / 2) /
      (max(sample_sizes) - min(sample_sizes))
  }
  
  
  # --------------------------------------------------------
  # ---------- Normality check function --------------------
  # --------------------------------------------------------
  
  normality_test <- function(data, test = "SW", alpha = 0.05) {
    pvals <- NULL
    all_normality_satisfied <- NULL
    
    # Case 1: Single numeric vector
    if (is.numeric(data) && is.atomic(data) && is.null(dim(data))) {
      pval <- generate_tests(data, test = test)$p.value
      pvals <- pval
      all_normality_satisfied <- pval > alpha
    }
    
    # Case 2: List of numeric vectors
    else if (is.list(data) && !is.data.frame(data)) {
      pvals <- sapply(data, function(sample) {
        generate_tests(sample, test = test)$p.value
      })
      names(pvals) <- names(data) %||% paste0("Sample", seq_along(pvals))
      all_normality_satisfied <- all(pvals > alpha)
    }
    
    # Case 3: Wide-format data frame
    else if (is.data.frame(data) && all(sapply(data, is.numeric))) {
      pvals <- sapply(as.list(data), function(sample) generate_tests(sample, test = test)$p.value)
      names(pvals) <- names(data)
      all_normality_satisfied <- all(pvals > alpha)
    }
    
    # Case 4: Long-format with group labels
    else if ((is.data.frame(data) || is.matrix(data)) && ncol(data) >= 2) {
      # We assume first column is group, second is value
      grouped_samples <- split(data[[2]], data[[1]])
      pvals <- sapply(grouped_samples, function(sample) {
        generate_tests(sample, test = test)$p.value
      })
      all_normality_satisfied <- all(pvals > alpha)
    }
    else {
      stop("Unsupported input type: must be numeric vector, list of vectors, or grouped data.")
    }
    
    return(list(p_values = pvals, normality_satisfied = all_normality_satisfied))
  }
  # ------------------------------------------------------------------------------ 
  # ------------------ power/Type I error function -------------------------------
  user_interface <- function(
    gen_data = generate_regression_data,
    fn_to_get_norm_obj = regression_residuals,
    fn_for_norm_test = normality_test,
    fn_for_ds_test_1 = simple_linear_regression,
    fn_for_ds_test_2 = bootstrap_regression,
    paras = NULL,
    alpha = 0.05,
    test_method = "SW",
    ...
  ) {
    # generate dataset
    data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
    
    # residuals & normality check
    resid     <- fn_to_get_norm_obj(data)
    norm_res  <- fn_for_norm_test(data = resid, test = test_method, alpha = alpha)
    
    # choose test
    if (isTRUE(norm_res$normality_satisfied)) {
      result <- fn_for_ds_test_1(data)
    } else {
      result <- fn_for_ds_test_2(data, ...)
    }
    
    return(result$p.value)
  }
  
  calculate_power <- function(
    sample_sizes    = c(10, 20, 30, 40, 50),
    n_sim           = 1e3,
    alpha           = 0.05,
    n_boot          = NULL,
    gen_data        = generate_regression_data,
    get_parameters  = function(n) list(n = n),
    fn_to_get_norm_obj   = regression_residuals,
    fn_for_norm_test     = normality_test,
    fn_for_ds_test_1     = simple_linear_regression,
    fn_for_ds_test_2     = bootstrap_regression,
    test_method    = "SW",
    mode           = c("adaptive", "parametric", "nonparametric"),
    ...
  ) {
    mode <- match.arg(mode)
    power_results <- numeric(length(sample_sizes))
    names(power_results) <- paste0("n=", sample_sizes)
    
    for (i in seq_along(sample_sizes)) {
      n          <- sample_sizes[i]
      rejections <- 0
      
      for (sim in seq_len(n_sim)) {
        paras <- get_parameters(n)
        # generate data
        data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
        
        p_value <- tryCatch({
          if (mode == "adaptive") {
            user_interface(
              gen_data           = gen_data,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras              = paras,
              alpha              = alpha,
              test_method        = test_method,
              ...
            )
          } else if (mode == "parametric") {
            fn_for_ds_test_1(data)$p.value
          } else {
            fn_for_ds_test_2(data, ...)$p.value
          }
        }, error = function(e) NA)
        
        if (!is.na(p_value) && p_value < alpha) {
          rejections <- rejections + 1
        }
      }
      
      power_results[i] <- rejections / n_sim
      cat("Completed n =", n, "| mode =", mode, "| Power =", power_results[i], "\n")
    }
    
    return(power_results)
  }
  # ------------------------------------------------------------------------------
  # ------- Main simulation results: power/Typer I error analysis ---------------
  results <- eventReactive(input$run, {
    sizes <- seq(input$sample_sizes[1], input$sample_sizes[2], by = 5)
    funcs <- user_funcs()
    if (is.null(funcs)) return(NULL)
    
    withProgress(message = "Running simulation...", value = 0, {
      # Initialize progress tracking
      n_steps <- 3
      step <- 1/n_steps
      
      # Run power analysis for the three test methods
      power_param <- calculate_power(
        sample_sizes       = sizes,
        n_sim              = input$n_sim,
        alpha              = input$alpha,
        gen_data           = funcs$gen_data,
        get_parameters     = funcs$get_parameters,
        fn_to_get_norm_obj = funcs$fn_get_norm_obj,
        fn_for_norm_test   = normality_test,
        fn_for_ds_test_1   = funcs$fn_ds_test_1,
        fn_for_ds_test_2   = funcs$fn_ds_test_2,
        test_method        = input$norm_test,
        mode               = "parametric"
      )
      
      incProgress(step, detail = paste("Step 1 of", n_steps, "complete"))
      
      power_nonpar <- calculate_power(
        sample_sizes       = sizes,
        n_sim              = input$n_sim,
        alpha              = input$alpha,
        gen_data           = funcs$gen_data,
        get_parameters     = funcs$get_parameters,
        fn_to_get_norm_obj = funcs$fn_get_norm_obj,
        fn_for_norm_test   = normality_test,
        fn_for_ds_test_1   = funcs$fn_ds_test_1,
        fn_for_ds_test_2   = funcs$fn_ds_test_2,
        test_method        = input$norm_test,
        mode               = "nonparametric"
      )
      
      incProgress(step, detail = paste("Step 2 of", n_steps, "complete"))
      
      power_adapt <- calculate_power(
        sample_sizes       = sizes,
        n_sim              = input$n_sim,
        alpha              = input$alpha,
        gen_data           = funcs$gen_data,
        get_parameters     = funcs$get_parameters,
        fn_to_get_norm_obj = funcs$fn_get_norm_obj,
        fn_for_norm_test   = normality_test,
        fn_for_ds_test_1   = funcs$fn_ds_test_1,
        fn_for_ds_test_2   = funcs$fn_ds_test_2,
        test_method        = input$norm_test,
        mode               = "adaptive"
      )
      
      incProgress(step, detail = paste("Step 3 of", n_steps, "complete"))
      
      # calculate auc
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
    })
  })
  
  
  # -------------------------------------------------------------------
  # ----------------- Plot Power/Type I error results -----------------
  plot_power_curve <- function(res, metric, test_type, norm_test) {
    is_power   <- metric == "Power"
    threshold  <- if (is_power) 0.8 else 0.05
    legend_pos <- if (is_power) "bottomright" else "topleft"
    y_limits   <- if (is_power) c(0, 1) else c(0, 0.1)
    y_label    <- if (is_power) "Power" else "Type I error"
    
    sizes    <- res$sizes
    p_param  <- res$power_param
    p_nonpar <- res$power_nonpar
    p_adapt  <- res$power_adapt
    auc_vals <- res$auc$AUC
    
    main_title <- paste(metric, "vs Sample Size")
    sub_title  <- paste("Test Type:", test_type,
                        "| Normality Test:", norm_test)
    
    plot(sizes, p_param,
         type = "b",
         pch = 19,
         col = "red",
         ylim = y_limits,
         xlab = "Sample Size",
         ylab = y_label,
         main = main_title,
         sub = sub_title)
    
    lines(sizes, p_nonpar, type = "b", pch = 17, col = "blue")
    lines(sizes, p_adapt,  type = "b", pch = 15, col = "green")
    
    abline(h = threshold, lty = 2, col = "gray50")
    
    legend(legend_pos,
           legend = paste0(c("Test1(Parametric)", "Test2(Nonparametric)", "Test3(Adaptive)"),
                           " (AUC=", formatC(auc_vals, format = "f", digits = 4), ")"),
           title   = paste("AUC for", metric),
           col     = c("red", "blue", "green"),
           pch     = c(19, 17, 15),
           lty     = 1,
           bty     = "o",
           box.lwd = 1.5)
  }
  
  # ----------------- Execute the plot function ---------------------------
  output$power_plot <- renderPlot({
    res <- results()
    if (is.null(res)) return(NULL)
    
    metric    <- isolate(input$metric)
    test_type <- isolate(input$test_type)
    norm_test <- isolate(input$norm_test)
    
    plot_power_curve(res, metric, test_type, norm_test)
  })
  
  # add caption
  output$power_caption <- renderUI({
    req(results())
    HTML(paste(
      "Showing ", input$metric, " comparison for ", input$test_type, 
      " among parametric, nonparamteric and an adaptive procedure based on ", 
      input$norm_test, " for normality for sample sizes from ", paste(input$sample_sizes[1], "to", input$sample_sizes[2]), " at significance level (α)  of", input$alpha
    ))
  })  
  
  # ------------------ Generate AUC table ----------------------------------
  output$auc_table <- renderTable({
    res <- results()
    if (is.null(res)) return(NULL)
    res$auc
  }, striped = TRUE, digits = 4)
  
  # -----------------------------------------------------------------------------
  # ---------------------- Function to compute FPR and TPR 
  # -----------------------------------------------------------------------------
  ROC_curve_function <- function(sample_size, alpha_pretest, tests, Nsim = 1e3) {
    req(user_funcs())
    funcs <- user_funcs()
    
    FPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
    TPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
    rownames(FPR) <- rownames(TPR) <- tests
    colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha_pretest)
    
    for (i in seq_along(tests)) {
      test_name <- tests[i]
      
      for (j in seq_along(alpha_pretest)) {
        alpha <- alpha_pretest[j]
        pval_norm <- pval_non_normal <- numeric(Nsim)
        
        for (k in 1:Nsim) {
          normal_data <- generate_data(sample_size, dist = "Normal", par = NULL)
          non_normal_data <- generate_data(sample_size, dist = "LogNormal", par = NULL)
          
          pval_norm[k] <- generate_tests(normal_data, test_name)$p.value
          pval_non_normal[k] <- generate_tests(non_normal_data, test_name)$p.value
        }
        
        FPR[i, j] <- mean(pval_norm < alpha)
        TPR[i, j] <- mean(pval_non_normal < alpha)
      }
    }
    
    return(list(FPR = FPR, TPR = TPR, alpha = alpha_pretest))
  }
  
 
  # ROC curve for normality test methods ---------------
  plot_ROC <- function(FPR, TPR, tests_to_plot = rownames(FPR), alpha = NULL, title = "ROC Curves") {
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
      lines(FPR[test, ], TPR[test, ], col = colors[i], lwd = 2)
      points(FPR[test, ], TPR[test, ], col = colors[i], pch = plot_chars[i], cex = 0.75)
    }
    abline(0, 1, lty = 2, col = "gray")
    legend("bottomright", 
           legend = tests_to_plot, 
           col = colors, 
           pch = plot_chars, 
           lwd = 2, 
           title = "Tests")
  }
  
  # Generate the ROC curve plot
  output$roc_curve_plot <- renderPlot({
    res <- roc_results()
    if (is.null(res)) return(NULL)
    plot_ROC(FPR = res$FPR, 
             TPR = res$TPR, 
             tests_to_plot = res$tests,
             alpha = res$alpha, 
             title = "ROC Curves for Normality Tests")
  })
  
  # Add caption
  output$roc_caption <- renderUI({
    req(roc_results())
    HTML(paste0(
      "ROC Curve Comparing ", paste(input$roc_tests, collapse = ", "), 
      " normality test methods for Sample size of ", input$roc_sample_size, 
      " for", input$roc_Nsim, " iterations"
    ))
  })
  
  # ------------------------------------------------------------
  # Function ROC like curves for downstream tests
  # ------------------------------------------------------------
  # Function to generate all required p-values 
  generate_pval <- function(sample_size, N, dist, effect_size) {
    # Initialize vectors
    pval_t.test_H0 <- pval_wilcox.test_H0 <- numeric(N)
    pval_t.test_H1 <- pval_wilcox.test_H1 <- numeric(N)
    p_sw_x <- p_sw_y <- numeric(N)
    
    # Initialize progress bar
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    
    for (i in 1:N) {
      x <- generate_data(sample_size, dist)
      y <- generate_data(sample_size, dist)
      
      # Pretest (Shapiro-Wilk)
      p_sw_x[i] <- shapiro.test(x)$p.value
      p_sw_y[i] <- shapiro.test(y)$p.value
      
      # Type I error p-values (under H0)
      pval_t.test_H0[i] <- t.test(x, y, var.equal = TRUE)$p.value
      pval_wilcox.test_H0[i] <- wilcox.test(x, y, exact = FALSE)$p.value
      
      # Power p-values (under H1)
      pval_t.test_H1[i] <- t.test(x, y + effect_size, var.equal = TRUE)$p.value
      pval_wilcox.test_H1[i] <- wilcox.test(x, y + effect_size, exact = FALSE)$p.value
      
      # Update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    return(list(
      p_sw_x = p_sw_x,
      p_sw_y = p_sw_y,
      pval_t.test_H0 = pval_t.test_H0,
      pval_wilcox.test_H0 = pval_wilcox.test_H0,
      pval_t.test_H1 = pval_t.test_H1,
      pval_wilcox.test_H1 = pval_wilcox.test_H1
    ))
  }
  
  # ------------------------------------------------------------
  # Function to perform complete analysis
  perform_analysis <- function(sample_size, N, distributions, effect_size, alpha_pretest, test_alpha) {
    results <- list()
    power_results <- list()
    error_results <- list()
    prob_sw_test <- list()
    
    for (dist in distributions) {
      cat("Processing distribution:", dist, "\n")
      results[[dist]] <- generate_pval(sample_size, N, dist, effect_size)
      
      # Calculate fixed test results
      power_results[[dist]] <- list(
        power_t.test = mean(results[[dist]]$pval_t.test_H1 < test_alpha),
        power_wilcox.test = mean(results[[dist]]$pval_wilcox.test_H1 < test_alpha)
      )
      
      error_results[[dist]] <- list(
        error_t.test = mean(results[[dist]]$pval_t.test_H0 < test_alpha),
        error_wilcox.test = mean(results[[dist]]$pval_wilcox.test_H0 < test_alpha)
      )
      
      # Initialize adaptive test results
      power_results[[dist]]$adaptive_wilcox <- numeric(length(alpha_pretest))
      error_results[[dist]]$adaptive_wilcox <- numeric(length(alpha_pretest))
      prob_sw_test[[dist]]$pr_sw_vec <- numeric(length(alpha_pretest))
      
      # Evaluate adaptive test at each alpha_pretest level
      for (j in seq_along(alpha_pretest)) {
        alpha <- alpha_pretest[j]
        
        # Decision rule: use t-test if both samples appear normal
        decision_region <- results[[dist]]$p_sw_x > alpha & results[[dist]]$p_sw_y > alpha
        
        # Adaptive test p-values
        adaptive_error_pvals <- ifelse(decision_region,
                                       results[[dist]]$pval_t.test_H0,
                                       results[[dist]]$pval_wilcox.test_H0)
        
        adaptive_power_pvals <- ifelse(decision_region,
                                       results[[dist]]$pval_t.test_H1,
                                       results[[dist]]$pval_wilcox.test_H1)
        
        # Store results
        error_results[[dist]]$adaptive_wilcox[j] <- mean(adaptive_error_pvals < test_alpha)
        power_results[[dist]]$adaptive_wilcox[j] <- mean(adaptive_power_pvals < test_alpha)
        prob_sw_test[[dist]]$pr_sw_vec[j] <- mean(!decision_region)
      }
    }
    
    return(list(
      results = results,
      power_results = power_results,
      error_results = error_results,
      prob_sw_test = prob_sw_test
    ))
  }
  
  # ------------------------------------------------------------
  # Function to generate ROC data for power vs type I error
  generate_roc_data <- function(sample_size, N, distributions, effect_size, sig_levels) {
    
    roc_data <- data.frame()
    
    for (dist in distributions) {
      cat("Generating ROC data for:", dist, "\n")
      res <- generate_pval(sample_size, N, dist, effect_size)
      
      for (alpha in sig_levels) {
        # Compute power & Type I error
        power_t <- mean(res$pval_t.test_H1 < alpha)
        power_w <- mean(res$pval_wilcox.test_H1 < alpha)
        error_t <- mean(res$pval_t.test_H0 < alpha)
        error_w <- mean(res$pval_wilcox.test_H0 < alpha)
        
        # Append results
        roc_data <- rbind(roc_data, 
                          data.frame(
                            Distribution = dist,
                            Method = "t-test",
                            Alpha = alpha,
                            Power = power_t,
                            TypeIError = error_t
                          ),
                          data.frame(
                            Distribution = dist,
                            Method = "Wilcoxon",
                            Alpha = alpha,
                            Power = power_w,
                            TypeIError = error_w
                          ))
      }
    }
    
    return(roc_data)
  }
  
  # ------------------------------------------------------------
  # Function to compute ROC-like metrics
  compute_roc_metrics <- function(power_results, error_results) {
    EPG <- power_results$LogNormal$adaptive_wilcox - power_results$LogNormal$power_t.test
    EPL <- power_results$Normal$power_t.test - power_results$Normal$adaptive_wilcox
    
    Expected_Inflation_lognormal <- error_results$LogNormal$adaptive_wilcox - error_results$LogNormal$error_t.test
    Expected_Inflation_normal <- error_results$Normal$adaptive_wilcox - error_results$Normal$error_t.test
    
    # Point estimates for benchmark comparison
    EPG_lognormal <- power_results$LogNormal$power_wilcox.test - power_results$LogNormal$power_t.test
    EPL_normal <- power_results$Normal$power_t.test - power_results$Normal$power_wilcox.test
    
    return(list(
      EPL = EPL,
      EPG = EPG,
      Expected_Inflation_lognormal = Expected_Inflation_lognormal,
      Expected_Inflation_normal = Expected_Inflation_normal,
      EPG_lognormal = EPG_lognormal,
      EPL_normal = EPL_normal
    ))
  }
  # ------------------------------------------------------------------------------
  # --------------- Plots the above results --------------------------------------
  # ------------------------------------------------------------------------------
  # Power and Type I error trade-off plots
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
  
  # ROC like Curve for downstream test --------------------------
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

  # Calculate ROC results
    roc_results <- eventReactive(input$run_roc, {
    # inputs
    tests <- input$roc_tests
    Nsim <- input$roc_Nsim
    n <- input$roc_sample_size
    # normality test alpha
    alpha_pretest <- seq(from = input$alpha_range[1], to   = input$alpha_range[2], by   = input$alpha_by)
    # downstream test alpha
    sig_level <- seq(from = input$sig_range[1], to   = input$sig_range[2], by   = input$sig_by)
    
    withProgress(message = "Running ROC Analysis...", value = 0, {
      roc_res <- ROC_curve_function(sample_size = n,
                                    alpha_pretest = alpha_pretest,
                                    tests = tests,
                                    Nsim = Nsim)
      incProgress(1/3)
      
      sim_data <- perform_analysis(sample_size = n,
                                   N = Nsim,
                                   distributions = c("Normal", "LogNormal"),
                                   effect_size = input$effect_size,
                                   alpha_pretest = alpha_pretest,
                                   test_alpha = input$test_alpha)
      incProgress(1/3)
      
      metrics <- compute_roc_metrics(sim_data$power_results, sim_data$error_results)
      incProgress(1/3)
      
      roc_data <- generate_roc_data(n, Nsim,
                                    c("Normal", "LogNormal"),
                                    effect_size = input$effect_size,
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
  
  # ----------------------Display the plots above -------------------------------
  # Power/Type I error trade off plots
  output$power_error_plot <- renderPlot({
    res <- roc_results()
    if (is.null(res)) return(NULL)
    plot_power_error_tradeoff(alpha_pretest = res$alpha, metrics = res$metrics)
  })
  
  # add caption
  output$tradeoff_caption <- renderUI({
    req(roc_results())
    HTML(
      paste0("Row 1 column 1 and row 1 column 2 respectively showing Expected Power loss 
             and expected power gains for pretesting for normality. Row 1 column 1 and row 1 
             column 2 respectively showing Expected Power loss and expected power gains for 
             pretesting for normality for effect size of ", input$effect_size, " significant 
             level for normality pretest ranging from ", input$alpha_range[1], " to ", input$alpha_range[2], 
             " at intervals of ", input$alpha_by
      ))
  })
  
  # Power vs Type I error ROC like curve 
  output$roc_like_metrics_plot <- renderPlot({
    res <- roc_results()
    if (is.null(res)) return(NULL)
    plot_roc_like_curves(metrics = res$metrics)
  })
  
  # add caption
  output$roc_like_caption <- renderUI({
    req(roc_results())
    HTML(paste0(
      "ROC like curve of EPG vs EPL and ROC like curve for 
      expected inflation of Type I error vs expected deflation 
      of Type I error for sample size ", input$roc_sample_size, 
      " effect size of ", input$effect_size, " for a downstream 
      significance level range from ", input$sig_range[1], " to ", 
      input$sig_range[2], " varying by ", input$sig_by.
    ))
  })
  
}


# Run the application
shinyApp(ui = ui, server = server)









