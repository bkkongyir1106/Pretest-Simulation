library(shiny)
# Load helper source files
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
setwd("~/Desktop/OSU/Research/Pretest-Simulation/User_interface/User_framework_ShinyApp")

ui <- fluidPage(
  titlePanel("User Framework for Statistical Analysis"),
  
  # Inject custom CSS to tighten up sidebar padding/margins
  tags$head(
    tags$style(HTML("
      .sidebar .well { padding: 6px; }
      .shiny-input-container { margin-bottom: 6px; }
    "))
  ),
  
  sidebarLayout(
    # make sidebar narrower: 3/12ths of the page instead of 4
    sidebarPanel(
      selectInput("test_type", "Test Type:",
                  choices = c("Independent t-test" = "ttest",
                              "Regression" = "regression",
                              "ANOVA" = "anova")),
      selectInput("metric","Analysis Metric:", 
                  choices = c("Power", "Type I Error")),
      # Load all user input function files here
      fileInput("gen_data_file",   "Upload Data Generation Function"),
      fileInput("get_params_file", "Upload Parameter Function"),
      fileInput("norm_obj_file",   "Upload Normality Object Function"),
      fileInput("test1_file",      "Upload Test 1 Function"),
      fileInput("test2_file",      "Upload Test 2 Function"),
      # Select Normality test inputs
      selectInput("norm_test", "Normality Test:",
                  choices = c("Shapiro-Wilk" = "SW",
                              "Shapiro-Francia" = "SF",
                              "Kolmogorov-Smirnov" = "KS",
                              "Jarque-Bera" = "JB",
                              "D'Agostino Pearson" = "DAP",
                              "Anderson-Darling" = "AD",
                              "Cramer-Von-Mises" = "CVM")
      ),
      # More input selection buttons
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
      textInput("effect_size_power",  
                label = "Effect Size for power",
                value = 0.5),
      
      numericInput("n_sim", "Number of Simulations:", 
                   value = 1000, 
                   min = 10, 
                   max = 100000),
      # action button to run simulation
      actionButton("run", "Run Simulation", class = "btn-primary"),
      
      # Resize sidebar Panel by changing width
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel(title = "ROC Curve Analysis",
                 fluidRow(
                   column(
                     4,
                     selectInput("distributions", "Select Distributions for ROC:",
                                 choices = c("Normal", "LogNormal", "Uniform", "Exponential", "Chi-Square", "Beta", "Gamma"),
                                 selected = c("Normal", "Uniform"),
                                 multiple = TRUE),
                     
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
                     
                     textInput("effect_size", "Means under H₁ (comma-separated):",
                               value = 0.5),
                     
                     numericInput("test_alpha", "Fixed Significance Level (α):", 
                                  value = 0.05, 
                                  min =  0.0, 
                                  max = 0.1),
                     
                     sliderInput("alpha_range", HTML("&alpha;<sub>pretest</sub> range"),
                                 min = 0.001, max = 1, value = c(0.001, 1), step = 0.001),
                     numericInput("alpha_by", HTML("&alpha;<sub>pretest</sub> step"),
                                  value = 0.05, min = 0.001, max = 0.5, step = 0.001),
                     
                     sliderInput("sig_range", "sig_level range",
                                 min = 0.001,  max = 1, value = c(0.001, 1), step = 0.001),
                     numericInput("sig_by", "sig_level step",
                                  value = 0.05, min = 0.001, max = 0.5, step = 0.001),
                     
                     actionButton("run_roc", "Run ROC Analysis", class = "btn-success")
                   ),
                   column(8,
                          tabsetPanel(
                            tabPanel(title = "Normality Test Methods", 
                                     plotOutput("roc_curve_plot" , height = "600px"),
                                     uiOutput("roc_caption")
                            ),
                            tabPanel(title = "Power/Error Tradeoff", 
                                     plotOutput("power_error_plot" , height = "600px"),
                                     uiOutput("tradeoff_caption")
                            ),
                            tabPanel(title = "ROC-like Metrics", 
                                     plotOutput("roc_like_metrics_plot" , height = "400px"),
                                     uiOutput("roc_like_caption"))
                          )
                   )
                 )
        ),
        # tabPanel(title = "Normality Test Methods", 
        #          plotOutput("norm_test_roc_plot"),
        #          uiOutput("norm_test_roc_caption")
        # ),
        tabPanel(title = "Power/Type I error Curves", 
                 plotOutput("power_plot"),
                 uiOutput("power_caption")
        ),
        tabPanel(title = "AUC Results", 
                 tableOutput("auc_table"),
                 uiOutput("auc_table_caption")
        ),
        
        # Create your documentation here
        tabPanel(title = "Documentation",
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
        ),
      )
    )
  )
)

# ------------------------------------------------------------------------------
# ------------------------------ Server Function -------------------------------
# ------------------------------------------------------------------------------
server <- function(input, output, session) {
 
  # Create an environment to store user-defined functions
  user_env <- new.env()
  
  # Reactive value to track if functions are loaded
  functions_loaded <- reactiveVal(FALSE)
  
  # Observe file uploads and load functions
  observe({
    req(input$gen_data_file, input$get_params_file, input$norm_obj_file, 
        input$test1_file, input$test2_file)
    
    tryCatch({
      # Source each file into the user environment
      source(input$gen_data_file$datapath, local = user_env)
      source(input$get_params_file$datapath, local = user_env)
      source(input$norm_obj_file$datapath, local = user_env)
      source(input$test1_file$datapath, local = user_env)
      source(input$test2_file$datapath, local = user_env)
      
      # Verify all required functions are present
      required_funcs <- c("gen_data", "get_parameters", "fn_get_norm_obj",
                          "fn_for_ds_test_1", "fn_for_ds_test_2")
      missing_funcs <- setdiff(required_funcs, ls(envir = user_env))
      
      if (length(missing_funcs) > 0 ){
        stop(paste("Missing functions:", paste(missing_funcs, collapse = ", ")))
      }
      
      functions_loaded(TRUE)
      showNotification("All functions loaded successfully!", type = "message")
    }, error = function(e) {
      functions_loaded(FALSE)
      showNotification(paste("Error loading functions:", e$message), type = "error")
    })
  })
  
  #----------   function for computing auc ----------------------------
  compute_area <- function(sample_sizes, y) {
    if (is.null(y) || length(y) != length(sample_sizes)) return(NA_real_)
    sum(diff(sample_sizes) * (head(y, -1) + tail(y, -1)) / 2) /
      (max(sample_sizes) - min(sample_sizes))
  }
  
  # --------------------------------------------------------
  # ---------- General Normality check function ------------
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
  # ------------------ Downstream Test function -------------------------------
  # ------------------------------------------------------------------------------
  ds_test_function <- function(
    gen_data = gen_data,
    get_parameters = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test = normality_test,
    fn_for_ds_test_1 = fn_for_ds_test_1,
    fn_for_ds_test_2 = fn_for_ds_test_2,
    paras            = NULL,
    alpha = 0.05,
    norm_test_method = "SW",
    ...
  ) {
    # generate dataset
    data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
    
    # get normality test object
    normality_test_object     <- fn_to_get_norm_obj(data)
    normality_test_pval_ds_test  <- fn_for_norm_test(data = normality_test_object, test = norm_test_method, alpha = alpha)
    
    # choose test
    if (isTRUE(normality_test_pval_ds_test$normality_satisfied)) {
      ds_norm_test_method <- fn_for_ds_test_1(data)
    } else {
      ds_norm_test_method <- fn_for_ds_test_2(data, ...)
    }
    
    return(ds_norm_test_method$p.value)
  }
  
  # ------------------------------------------------------------------------------
  perform_ds_func <- function(
    sample_sizes       = c(10, 20, 30, 40, 50),
    Nsim               = 1e3,
    alpha              = 0.05,
    n_boot             = NULL,
    gen_data           = gen_data,
    get_parameters     = function(n) list(n = n),
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = fn_for_norm_test,
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = "SW",
    ds_test_methods    = c("parametric", "nonparametric", "adaptive"),
    ...
  ) {
    ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
    results <- list()
    
    # # Provide Effect Size Here
    effect_size <- as.numeric(trimws(strsplit(input$effect_size_power, ",")[[1]]))
    
    # Create progress bar
    progress <- Progress$new(session, min = 0, max = length(sample_sizes) * length(ds_test_methods))
    on.exit(progress$close())
    
    for (method in ds_test_methods) {
      
      progress$set(message = paste("Running", method), detail = "Initializing...")
      
      ds_test_results <- numeric(length(sample_sizes))
      names(ds_test_results) <- paste0("n=", sample_sizes)
      
      for (i in seq_along(sample_sizes)) {
        n <- sample_sizes[i]
        rejections <- 0
        
        # Update progress
        progress$inc(1, detail = paste(method, "- n =", n))
        
        # Get parameters with actual sample size
        paras <- user_env$get_parameters(n = n, means = effect_size)
        
        for (sim in seq_len(Nsim)) {
          
          data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
          
          p_value <- tryCatch({
            if (method == "adaptive") {
              ds_test_function(
                gen_data           = gen_data,
                get_parameters     = get_parameters,
                fn_to_get_norm_obj = fn_to_get_norm_obj,
                fn_for_norm_test   = fn_for_norm_test,
                fn_for_ds_test_1   = fn_for_ds_test_1,
                fn_for_ds_test_2   = fn_for_ds_test_2,
                paras              = paras,
                alpha              = alpha,
                norm_test_method   = norm_test_method,
                ...
              )
            } else if (method == "parametric") {
              fn_for_ds_test_1(data)$p.value
            } else {
              fn_for_ds_test_2(data, ...)$p.value
            }
          }, error = function(e) NA)
          
          if (!is.na(p_value) && p_value < alpha) {
            rejections <- rejections + 1
          }
        }
        
        ds_test_results[i] <- rejections / Nsim
        cat("Completed n =", n, "| Method =", method, "| Power =", ds_test_results[i], "\n")
      }
      
      results[[method]] <- ds_test_results
    }
    
    return(results)
  }
  
  # ------------------------------------------------------------------------------
  # ----------------- Main simulation: power/Type I error analysis ---------------
  
  # store simulation results in a reactive container 
  results <- reactiveVal()
  observeEvent(input$run, {
    # Ensure functions are loaded
    req(functions_loaded())  
    
    # Add progress bar on shiny
    withProgress(message = "Running simulation...", value = 0, {
      
    # define the test methods
    ds_test_methods <- c("parametric", "nonparametric", "adaptive")
    # define the sequence of sample sizes
    sizes <- seq(input$sample_sizes[1], input$sample_sizes[2], by = 5)
    
    ds_test_results <- tryCatch({
      perform_ds_func(
        sample_sizes       = sizes,
        Nsim               = input$n_sim,
        alpha              = input$alpha,
        gen_data           = user_env$gen_data,
        get_parameters     = user_env$get_parameters,
        fn_to_get_norm_obj = user_env$fn_get_norm_obj,
        fn_for_norm_test   = normality_test,
        fn_for_ds_test_1   = user_env$fn_for_ds_test_1,
        fn_for_ds_test_2   = user_env$fn_for_ds_test_2,
        norm_test_method   = input$norm_test,
        ds_test_methods    = ds_test_methods
      )
    }, error = function(e) {
      showNotification(paste("Simulation error:", e$message), type = "error")
      return(NULL)
    })
    
    # Access individual test results
    ds_test1_results  <- ds_test_results$parametric
    ds_test2_results <- ds_test_results$nonparametric
    ds_adaptive_results  <- ds_test_results$adaptive
    
    # calculate auc
    auc_df <- data.frame(
      Method = c("Test 1", "Test 2", "Test 3(Adaptive)"),
      AUC = c(
        compute_area(sizes, ds_test1_results),
        compute_area(sizes, ds_test2_results),
        compute_area(sizes, ds_adaptive_results)
      )
    )
    
    results(list(
      sizes = sizes,
      power_param = ds_test1_results,
      power_nonpar = ds_test2_results,
      power_adapt = ds_adaptive_results,
      auc = auc_df
    ))
    
    
    # ------------------------------------------------------------------------------
    # ----------------- Plot Power/Type I error -----------------------
    
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
      plot(sizes, p_param, 
           type = "b", 
           pch = 19,
           col  = "red",
           ylim = y_limits,
           xlab = "Sample Size", 
           ylab = y_label,
           main = main_title, 
           sub = sub_title)
      
      # add the others
      lines(sizes, p_nonpar, type = "b", pch = 17, col = "blue")
      lines(sizes, p_adapt,  type = "b", pch = 15, col = "green")
      
      # reference line
      abline(h = threshold, lty = 2, col = "gray50")
      
      # legend with AUC
      legend(legend_pos,
             legend = paste0(c("Test1(Parametric)","Test2(Nonparametric)","Test3(Adaptive)"),
                             " (AUC=", formatC(res$auc$AUC, format = "f", digits = 4),")"
             ),
             title    = paste("AUC for", input$metric),
             col      = c("red","blue","green"),
             pch      = c(19, 17, 15),
             lty      = 1,
             bty      = "o",        
             box.lwd  = 1.5         
      )
    })
    
    # add caption
    output$power_caption <- renderUI({
      req(results())
      HTML(paste(
        "Showing ", input$metric, " comparison for ", input$test_type, 
        " among parametric, nonparamteric and an adaptive procedure based on ", 
        input$norm_test, " for normality for sample sizes from ", paste(input$sample_sizes[1],
        "to", input$sample_sizes[2]), " at significance level (α)  of", input$alpha
      ))
    })
    
    
    # AUC table
    output$auc_table <- renderTable({
      res <- results()
      if (is.null(res)) return(NULL)
      res$auc
    }, striped = TRUE, digits = 4)
    
    incProgress(1)
    })
  })
  
  # =============================================================================
  # -----------------------------------------------------------------------------
  # Function to Create ROC Curves for Normality Test Methods
  # -----------------------------------------------------------------------------
  
  # create  a distribution selection function
  select_dist <- reactive({
    validate(
      need(length(input$distributions) == 2, "Please select exactly two distributions."),
      need(input$distributions[1] == "Normal", "First selected distribution must be 'Normal'.")
    )
    list(dist1 = input$distributions[1], dist2 = input$distributions[2]
    )
  })
  
  # validate functions
  roc_results <- reactiveVal(NULL)
  observeEvent(input$run_roc, {
    req(functions_loaded())
    # add progress bar to shiny app
    withProgress(message = "Running ROC Analysis...", value = 0, {
    alpha_seq <- seq(input$alpha_range[1], input$alpha_range[2], by = input$alpha_by)
    
    all_tests <- c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM", "SKEW")
    
    tryCatch({
      FPR <- TPR <- matrix(0, nrow = length(all_tests), ncol = length(alpha_seq),
                           dimnames = list(all_tests, paste0("alpha=", round(alpha_seq, 3))))
      
      for (i in seq_along(all_tests)) {
        test <- all_tests[i]
        for (j in seq_along(alpha_seq)) {
          alpha <- alpha_seq[j]
          reject_H0 <- reject_H1 <- numeric(input$roc_Nsim)
          
          for (k in seq_len(input$roc_Nsim)) {
            # Under H0 
            paras_H0 <- user_env$get_parameters(n = input$roc_sample_size, dist = select_dist()$dist1)
            data_H0 <- do.call(user_env$gen_data, paras_H0)
            norm_obj_H0 <- user_env$fn_get_norm_obj(data_H0)
            reject_H0[k] <- any(normality_test(norm_obj_H0, test, alpha)$p_values < alpha)
            
            # Under H1
            paras_H1 <- user_env$get_parameters(n = input$roc_sample_size, dist = select_dist()$dist2)
            data_H1 <- do.call(user_env$gen_data, paras_H1)
            norm_obj_H1 <- user_env$fn_get_norm_obj(data_H1)
            reject_H1[k] <- any(normality_test(norm_obj_H1, test, alpha)$p_values < alpha)
          }
          
          FPR[i, j] <- mean(reject_H0, na.rm = TRUE)
          TPR[i, j] <- mean(reject_H1, na.rm = TRUE)
        }
        
        incProgress(1 / length(all_tests))  # update per test
      }
      
      roc_results(list(FPR = FPR, TPR = TPR, alpha = alpha_seq))
    }, error = function(e) {
      showNotification(paste("ROC analysis error:", e$message), type = "error")
    })
  })
  
  })
  
  # function to create ROC curve
  output$roc_curve_plot <- renderPlot({
    roc_res <- roc_results()
    req(roc_res)
    # plot only selected tests
    selected_tests <- input$roc_tests
    palette <- rainbow(length(selected_tests))
    
    plot(0, 0, 
         type = "n", 
         xlim = c(0, 1), 
         ylim = c(0, 1),
         xlab = "False Positive Rate (FPR)", 
         ylab = "True Positive Rate (TPR)",
         main = "ROC Curve: Normality Tests")
    for (i in seq_along(selected_tests)) {
      test <- selected_tests[i]
      if (test %in% rownames(roc_res$FPR)) {
        lines(roc_res$FPR[test, ], roc_res$TPR[test, ], type = "b", lwd = 2, col = palette[i])
      }
    }
    
    # add abline
    abline(0, 1, lty = 2, col = "gray")
    # add legend
    legend("bottomright",
           legend = selected_tests,
           col = palette,
           lwd = 2,
           title = "Selected Normality Tests")
  })
  
  # add caption
  output$roc_caption <- renderUI({
    req(roc_results())
    HTML(paste("ROC analysis run using sample size", input$roc_sample_size,
               "and", input$roc_Nsim, "simulations for tests:",
               paste(input$roc_tests, collapse = ", ")))
  })
  
  # =============================================================================
  # -----------------------------------------------------------------------------
  # Function to compute p-values for normality test & downstream tests
  # -----------------------------------------------------------------------------
  generate_pval <- function(N, test = "SW", dist = NULL, ...) {
    pval_t.test_H0 <- pval_wilcox.test_H0 <- numeric(N)
    pval_t.test_H1 <- pval_wilcox.test_H1 <- numeric(N)
    norm_pvals_H0 <- norm_pvals_H1 <- vector("list", N)
    
    # 1. Parse the H1 means the user typed:
    effect_size <- as.numeric(trimws(strsplit(input$effect_size, ",")[[1]]))
    if (any(is.na(effect_size))) {
      stop("Invalid H₁ means: please enter comma-separated numbers, e.g. 0.2,0.5")
    }
    
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    
    for (i in 1:N) {
      # Under H0
      paras_H0 <- user_env$get_parameters(n = input$roc_sample_size)
      data_H0 <- do.call(user_env$gen_data, paras_H0)
      
      # Under H1
      paras_H1 <- user_env$get_parameters(n = input$roc_sample_size, means = effect_size)
      data_H1 <- do.call(user_env$gen_data, paras_H1)
      
      # Normality objects
      norm_obj_H0 <- user_env$fn_get_norm_obj(data_H0)
      norm_obj_H1 <- user_env$fn_get_norm_obj(data_H1)
      
      norm_H0 <- normality_test(norm_obj_H0, test = test, alpha = 0.05)
      norm_H1 <- normality_test(norm_obj_H1, test = test, alpha = 0.05)
      
      norm_pvals_H0[[i]] <- norm_H0$p_values
      norm_pvals_H1[[i]] <- norm_H1$p_values
      
      pval_t.test_H0[i] <- user_env$fn_for_ds_test_1(data_H0)$p.value
      pval_wilcox.test_H0[i] <- user_env$fn_for_ds_test_2(data_H0)$p.value
      pval_t.test_H1[i] <- user_env$fn_for_ds_test_1(data_H1)$p.value
      pval_wilcox.test_H1[i] <- user_env$fn_for_ds_test_2(data_H1)$p.value
      
      setTxtProgressBar(pb, i)
    }
    
    close(pb)
    
    return(list(
      pval_t.test_H0 = pval_t.test_H0,
      pval_wilcox.test_H0 = pval_wilcox.test_H0,
      pval_t.test_H1 = pval_t.test_H1,
      pval_wilcox.test_H1 = pval_wilcox.test_H1,
      norm_pvals_H0 = norm_pvals_H0,
      norm_pvals_H1 = norm_pvals_H1
    ))
  }
  
  
  # --------------- # Function to analyze trade-offs across alphas --------------
  perform_analysis <- function(N, distributions, alpha_pretest, test_alpha) {
    distributions <- c(select_dist()$dist1, select_dist()$dist2)
    ds_test_results <- list()
    error_ds_test <- list()
    power_ds_test <- list()
    
    for(dist in distributions) {
      cat("Processing distribution:", dist, "\n")
      ds_test_results[[dist]] <- generate_pval(N, test = "SW")
      
      # Calculate Type I error rates (under H0)
      error_ds_test[[dist]] <- list(
        error_t.test = mean(ds_test_results[[dist]]$pval_t.test_H0 < test_alpha),
        error_wilcox.test = mean(ds_test_results[[dist]]$pval_wilcox.test_H0 < test_alpha),
        adaptive_wilcox_error = numeric(length(alpha_pretest)))
      
      # Calculate Power (under H1)
      power_ds_test[[dist]] <- list(
        power_t.test = mean(ds_test_results[[dist]]$pval_t.test_H1 < test_alpha),
        power_wilcox.test = mean(ds_test_results[[dist]]$pval_wilcox.test_H1 < test_alpha),
        adaptive_wilcox_power = numeric(length(alpha_pretest)))
      
      # Pre-compute decisions for each alpha level
      for(j in seq_along(alpha_pretest)) {
        alpha <- alpha_pretest[j]
        
        # For Type I error (H0)
        use_t_test_H0 <- sapply(ds_test_results[[dist]]$norm_pvals_H0, function(x) all(x > alpha))
        adaptive_pvals_H0 <- ifelse(use_t_test_H0,
                                    ds_test_results[[dist]]$pval_t.test_H0,
                                    ds_test_results[[dist]]$pval_wilcox.test_H0)
        error_ds_test[[dist]]$adaptive_wilcox_error[j] <- mean(adaptive_pvals_H0 < test_alpha)
        
        # For Power (H1)
        use_t_test_H1 <- sapply(ds_test_results[[dist]]$norm_pvals_H1, function(x) all(x > alpha))
        adaptive_pvals_H1 <- ifelse(use_t_test_H1,
                                    ds_test_results[[dist]]$pval_t.test_H1,
                                    ds_test_results[[dist]]$pval_wilcox.test_H1)
        power_ds_test[[dist]]$adaptive_wilcox_power[j] <- mean(adaptive_pvals_H1 < test_alpha)
      }
    }
    
    return(list(
      error_ds_test = error_ds_test,
      power_ds_test = power_ds_test
    ))
  }
  
  # run analysis to get power and error 
  # run analysis to get power and error 
  analysis_ds_tests <- reactive({
    withProgress(message = "Performing trade-off analysis...", value = 0, {
      res <- perform_analysis(
        N            = input$roc_Nsim,
        distributions = c(select_dist()$dist1, select_dist()$dist2),
        alpha_pretest = seq(input$sig_range[1], input$sig_range[2], by = input$sig_by),
        test_alpha    = input$test_alpha
      )
      incProgress(1)
      res
    })
  })
  
  # ------------------------------------------------------------------------------
  #                         Function to compute ROC-like metrics
  # ------------------------------------------------------------------------------
  compute_roc_metrics <- function(error_ds_test, power_ds_test, test_alpha) {
    
    # Non-normal case
    power_t_test_nonnormal <- power_ds_test[[ 2 ]]$power_t.test
    adaptive_power_nonnormal <- power_ds_test[[ 2 ]]$adaptive_wilcox_power
    adaptive_error_nonnormal <- error_ds_test[[ 2 ]]$adaptive_wilcox_error
    EPG <- adaptive_power_nonnormal - power_t_test_nonnormal
    EDE <- adaptive_error_nonnormal - test_alpha
    
    # normal case
    power_t_test_normal <- power_ds_test[[ 1 ]]$power_t.test
    adaptive_power_normal <- power_ds_test[[ 1 ]]$adaptive_wilcox_power
    adaptive_error_normal <- error_ds_test[[ 1 ]]$adaptive_wilcox_error
    EPL <- power_t_test_normal - adaptive_power_normal
    EIE <- adaptive_error_normal - test_alpha
    
    # Point estimates for benchmark comparison
    power_gain <- power_ds_test[[ 2 ]]$power_wilcox.test - power_ds_test[[ 2 ]]$power_t.test
    power_loss <- power_ds_test[[ 1 ]]$power_t.test - power_ds_test[[ 1 ]]$power_wilcox.test
    
    return(list(EPL = EPL,
                EPG = EPG,
                EIE = EIE,
                EDE = EDE,
                power_gain = power_gain,
                power_loss = power_loss)
    )
  }
  
  # Compute metrics after running analysis
  metrics <- reactive({
    tryCatch({
      req(analysis_ds_tests())
      compute_roc_metrics(
        error_ds_test = analysis_ds_tests()$error_ds_test,
        power_ds_test = analysis_ds_tests()$power_ds_test,
        test_alpha    = input$test_alpha
      )
    }, error = function(e) {
      showNotification(paste("Error computing metrics:", e$message), type = "error")
      return(NULL)
    })
  })
  
  # --------------------------------------------------------------------
  # Power and Type I error trade-off plots
  plot_power_error_tradeoff <- function(alpha_pretest, metrics, file_name = NULL) {
    # if (!is.null(file_name)) pdf(file_name, width = 12, height = 10)
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
         metrics$EIE, 
         type = "l", 
         col = "orange", 
         lwd = 2,
         ylab = "Type I Error Inflation", 
         xlab = expression(alpha),
         main = "Inflation (Normal)")
    
    plot(alpha_pretest, 
         metrics$EDE, 
         type = "l", 
         col = "green", 
         lwd = 2,
         ylab = "Type I Error Inflation", 
         xlab = expression(alpha),
         main = "Inflation (LogNormal)")
    if (!is.null(file_name)) dev.off()
  }
  
  # Generate all plots
  alpha_pretest <- reactive(seq(from = input$alpha_range[1], to = input$alpha_range[2], by = input$alpha_by))
  output$power_error_plot <- renderPlot({
    req(metrics())
    plot_power_error_tradeoff(alpha_pretest = alpha_pretest(), metrics = metrics())
  })
  
  
  # --------------------------------------------------------------------
  # ----------- Function to ROC curve data for downstream test ---------
  generate_roc_tables <- function(N, distributions, sig_levels) {
    
    withProgress(message = 'Preparing ROC data...', value = 0, {
    
    distributions <- c(select_dist()$dist1, select_dist()$dist2)
    # data frames for Type I error and Power
    roc_data <- data.frame()
    
    total_steps <- length(distributions) * length(sig_levels)
    current_step <- 0
    
    for (dist in distributions) {
      cat("Generating ROC data for:", dist, "\n")
      res <- generate_pval(N, test = "SW", dist)
      
      for (alpha in sig_levels) {
        
        # Update progress
        current_step <- current_step + 1
        incProgress(1/total_steps, detail = paste(dist, "at α =", alpha))
        
        # Compute Type I error (under H0)
        error_t <- mean(res$pval_t.test_H0 < alpha)
        error_w <- mean(res$pval_wilcox.test_H0 < alpha)
        
        # Compute Power (under H1)
        power_t <- mean(res$pval_t.test_H1 < alpha)
        power_w <- mean(res$pval_wilcox.test_H1 < alpha)
        
        # Append results
        roc_data <- rbind(roc_data, 
                          data.frame(
                            Distribution = dist,
                            Method = "test 1",
                            Alpha = alpha,
                            Power = power_t,
                            TypeIError = error_t
                          ),
                          data.frame(
                            Distribution = dist,
                            Method = "test 2",
                            Alpha = alpha,
                            Power = power_w,
                            TypeIError = error_w
                          ))
      }
    }
    
    return(roc_data)
    })
  }
  
  # Run the roc data
  roc_data <- reactive({
  
    generate_roc_tables(N = input$roc_Nsim, 
                        distributions = distributions, 
                        sig_levels = seq(from = input$sig_range[1], to = input$sig_range[2], by = input$sig_by))
  })
  
  # --------------------- create the plots -------------------------------------
  plot_power_vs_error_ROC_curve <- function(roc_data, file_name) {
    # Add progress bar
    withProgress(message = 'Generating ROC curves...', value = 0, {
    # pdf(file_name, width = 10, height = 6)
    par(mfrow = c(1, length(unique(roc_data$Distribution))))
    
    methods <- unique(roc_data$Method)
    colors <- c("blue", "red")
    pch_vals <- c(19, 17)
    # Get unique distributions
    distributions <- unique(roc_data$Distribution)
    
    for (i in seq_along(distributions)) {
      dist <- distributions[i]
      # Update progress
      incProgress(1/length(distributions), detail = paste("Processing", dist))
      
      plot(NA, xlim = c(0, 1), 
           ylim = c(0, 1), 
           xlab = "Type I Error", 
           ylab = "Power",
           main = paste("ROC-like Curve (", dist, ")", sep = ""))
      
      for (m in seq_along(methods)) {
        method <- methods[m]
        data_subset <- subset(roc_data, Distribution == dist & Method == method)
        lines(data_subset$TypeIError, data_subset$Power, type = "l",col = colors[m], lwd = 2, pch = pch_vals[m])
      }
      legend("bottomright", legend = methods, col = colors, lwd = 2, pch = pch_vals, title = "Method")
    }
    
    #  dev.off()
    })
  }
  
  # create the plot
  output$roc_like_metrics_plot <- renderPlot({
    req(roc_data())
    plot_power_vs_error_ROC_curve(roc_data())
  })
  
  output$roc_like_caption <- renderUI({
    HTML("ROC-like plot showing the tradeoff of Power against Type I Error as a results of pretesting for normality.")
  })
  
  
}

shinyApp(ui, server)
