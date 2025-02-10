library(shiny)
library(ggplot2)
library(dplyr)

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
                   choices = c("One-Sample" = "one", "Two-Sample" = "two")),
      
      selectInput("test", "Select Test:",
                  choices = c("Kolmogorov-Smirnov (KS)" = "KS",
                              "Shapiro-Wilk (SW)" = "SW",
                              "Jarque-Bera (JB)" = "JB",
                              "D'Agostino-Pearson (DAP)" = "DAP",
                              "Anderson-Darling (AD)" = "AD",
                              "Shapiro-Francia (SF)" = "SF",
                              "Cramér-von Mises (CVM)" = "CVM")),
      
      numericInput("alpha", "Significance Level:", value = 0.05, min = 0.01, max = 0.2, step = 0.01),
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data Summary",
                 plotOutput("data_plot"),
                 tableOutput("summary_stats")),
        
        tabPanel("Test Results",
                 verbatimTextOutput("test_output"),
                 plotOutput("power_plot")),
        
        tabPanel("Bootstrap Analysis",
                 plotOutput("bootstrap_plot"),
                 tableOutput("bootstrap_table"))
      )
    )
  )
)
