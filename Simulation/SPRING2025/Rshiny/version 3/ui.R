# ui.R
library(shiny)

fluidPage(
  titlePanel("Statistical Test Analysis"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("data_source", "Data Source:",
                   choices = c("Generate Data" = "gen", "Upload Data" = "upload")),
      
      conditionalPanel(
        condition = "input.data_source == 'gen'",
        selectInput("dist", "Distribution:",
                    choices = c("Normal", "Chi-Square", "Gamma", "Exponential",
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
        checkboxGroupInput("norm_tests", "Select Normality Tests:",
                           choices = c("Kolmogorov-Smirnov (KS)" = "KS",
                                       "Shapiro-Wilk (SW)" = "SW",
                                       "Jarque-Bera (JB)" = "JB",
                                       "D'Agostino-Pearson (DAP)" = "DAP",
                                       "Anderson-Darling (AD)" = "AD",
                                       "Shapiro-Francia (SF)" = "SF",
                                       "Cramér-von Mises (CVM)" = "CVM"),
                           selected = "SW"),
        conditionalPanel(
          condition = "input.data_source == 'gen'",
          numericInput("N_norm", "Number of simulations (N):", value = 1000, min = 1)
        )
      ),
      
      conditionalPanel(
        condition = "input.test_type == 'location'",
        selectInput("loc_test", "Select Location Test:",
                    choices = c("t-test" = "t",
                                "Wilcoxon" = "Wilcox",
                                "Adaptive t/Wilcox" = "t_Wilcox",
                                "Permutation" = "perm",
                                "Custom" = "custom")),
        conditionalPanel(
          condition = "input.loc_test == 'custom'",
          fileInput("custom_test_file", "Upload Custom Test Function (R script)",
                    accept = c(".R"))
        ),
        radioButtons("samples", "Number of Samples:",
                     choices = c("One-Sample" = "one", "Two-Sample" = "two")),
        numericInput("effect_size", "Effect Size:", value = 0.5, step = 0.1),
        conditionalPanel(
          condition = "input.data_source == 'gen'",
          numericInput("N", "Number of iterations (N):", value = 1000, min = 1),
          numericInput("B", "Number of permutations (B):", value = 1000, min = 1)
        ),
        conditionalPanel(
          condition = "input.data_source == 'upload'",
          numericInput("n_bootstrap", "Number of bootstraps:", value = 1000, min = 1)
        )
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
