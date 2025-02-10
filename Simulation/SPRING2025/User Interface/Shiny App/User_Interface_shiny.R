if(!require("pacman")) install.packages("pacman")
pacman::p_load(shiny, ggplot2, dplyr, purrr, tidyr, e1071, tseries, nortest)
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

ui <- fluidPage(
  titlePanel("Statistical Test Analysis"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_source", "Data Source:",
                   choices = c("Generate Data" = "gen", "Upload Data" = "upload")),
      
      conditionalPanel(
        condition = "input.data_source == 'gen'",
        selectInput("dist", "Distribution:",
                    choices = c("Standard Normal", "Chi-Square", "Gamma", "Exponential",
                                "t", "Uniform", "Laplace", "Weibull", "LogNormal",
                                "Contaminated", "Pareto")),
        numericInput("n", "Sample Size:", value = 30, min = 5, max = 1000)
      ),
      
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("file", "Choose CSV File",
                  accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
      ),
      
      radioButtons("test_type", "Test Type:",
                   choices = c("Normality Test" = "normality",
                               "Location Test" = "location")),
      
      conditionalPanel(
        condition = "input.test_type == 'normality'",
        selectInput("norm_test", "Select Normality Test:",
                    choices = c("Kolmogorov-Smirnov (KS)" = "KS",
                                "Shapiro-Wilk (SW)" = "SW",
                                "Jarque-Bera (JB)" = "JB",
                                "D'Agostino-Pearson (DAP)" = "DAP",
                                "Anderson-Darling (AD)" = "AD",
                                "Shapiro-Francia (SF)" = "SF",
                                "Cramér-von Mises (CVM)" = "CVM"))
      ),
      
      conditionalPanel(
        condition = "input.test_type == 'location'",
        selectInput("loc_test", "Select Location Test:",
                    choices = c("t-test" = "t",
                                "Wilcoxon" = "Wilcox",
                                "Adaptive t/Wilcox" = "t_Wilcox",
                                "Permutation" = "perm",
                                "Bootstrap" = "boot")),
        radioButtons("samples", "Number of Samples:",
                     choices = c("One-Sample" = "one", "Two-Sample" = "two")),
        numericInput("effect_size", "Effect Size:", value = 0.5, step = 0.1)
      ),
      
      numericInput("alpha", "Significance Level:", value = 0.05, min = 0.01, max = 0.2, step = 0.01),
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Summary",
                 plotOutput("data_plot"),
                 tableOutput("summary_stats")),
        
        tabPanel("Test Results",
                 verbatimTextOutput("test_output")),
        
        tabPanel("Power Analysis",
                 plotOutput("power_plot"),
                 tableOutput("power_table"))
      )
    )
  )
)

server <- function(input, output) {
  
  data <- reactive({
    if(input$data_source == "gen"){
      Generate_data(datagen.type = 1, 
                    n = input$n, 
                    dist = input$dist,
                    two_samples = (input$samples == "two"))
    } else {
      req(input$file)
      Generate_data(datagen.type = 2,
                    two_samples = (input$samples == "two"))
    }
  })
  
  output$data_plot <- renderPlot({
    df <- data()
    if(input$test_type == "normality" || input$samples == "one"){
      ggplot(data.frame(x = df$x), aes(x)) + 
        geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue") +
        geom_density(color = "darkblue") +
        ggtitle("Data Distribution") +
        theme_minimal()
    } else {
      df_long <- data.frame(
        value = c(df$x, df$y),
        group = rep(c("Sample 1", "Sample 2"), each = length(df$x))
      )
      ggplot(df_long, aes(value, fill = group)) +
        geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
        ggtitle("Two-Sample Distributions") +
        theme_minimal()
    }
  })
  
  output$summary_stats <- renderTable({
    df <- data()
    if(input$test_type == "normality" || input$samples == "one"){
      data.frame(
        Statistic = c("Mean", "SD", "Skewness", "Kurtosis"),
        Value = c(mean(df$x), sd(df$x), skewness(df$x), kurtosis(df$x))
      )
    } else {
      data.frame(
        Sample = c("Sample 1", "Sample 2"),
        Mean = c(mean(df$x), mean(df$y)),
        SD = c(sd(df$x), sd(df$y)),
        Skewness = c(skewness(df$x), skewness(df$y)),
        Kurtosis = c(kurtosis(df$x), kurtosis(df$y))
      )
    }
  })
  
  test_result <- eventReactive(input$run, {
    df <- data()
    if(input$test_type == "normality") {
      generate_tests(df$x, input$norm_test)
    } else {
      if(input$samples == "one") {
        OneSample.test(df$x, input$loc_test, input$alpha, B = 1000)
      } else {
        TwoSample.test(df$x, df$y, input$loc_test, input$alpha, B = 1000)
      }
    }
  })
  
  output$test_output <- renderPrint({
    result <- test_result()
    if(input$test_type == "normality") {
      cat("Normality Test:", input$norm_test, "\n")
      cat("Test Statistic:", result$statistic, "\n")
      cat("P-value:", result$p.value, "\n")
    } else {
      cat("Location Test:", input$loc_test, "\n")
      cat("P-value:", result, "\n")
    }
    cat("Conclusion:", ifelse(result$p.value < input$alpha, 
                              "Reject Null Hypothesis", 
                              "Fail to Reject Null Hypothesis"))
  })
  
  output$power_plot <- renderPlot({
    req(input$test_type == "location")
    
    df <- data()
    sample_sizes <- seq(10, 100, by = 10)
    effect_size <- input$effect_size
    test_names <- c("t", "Wilcox", "t_Wilcox", "perm", "boot")
    
    if(input$data_source == "upload") {
      # Bootstrap analysis for uploaded data
      results <- map_dfr(test_names, ~{
        if(input$samples == "one") {
          pwr <- bootstrap_one_sample(df$x, effect_size, input$alpha, 
                                            n_bootstrap = 100, sample_sizes)
        } else {
          pwr <- bootstrap_two_sample_power(df$x, df$y, effect_size, input$alpha,
                                            n_bootstrap = 100, sample_sizes)
        }
        data.frame(
          SampleSize = sample_sizes,
          Value = pwr,
          Test = .x
        )
      })
    } else {
      # Simulation-based analysis for generated data
      results <- map_dfr(test_names, ~{
        pwr <- Calculate_power(alpha = input$alpha, N = 100, 
                               twosamples = (input$samples == "two"),
                               dist = input$dist, 
                               sample_size = sample_sizes, 
                               test = .x,
                               effect_size = effect_size, B = 50)
        data.frame(
          SampleSize = sample_sizes,
          Value = pwr,
          Test = .x
        )
      })
    }
    
    # Determine analysis type label
    analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
    
    ggplot(results, aes(SampleSize, Value, color = Test)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Comparison of Location Tests -", analysis_type),
           x = "Sample Size",
           y = analysis_type) +
      theme_minimal() +
      scale_color_brewer(palette = "Set1")
  })
}

shinyApp(ui = ui, server = server)
