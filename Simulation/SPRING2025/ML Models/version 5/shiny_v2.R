library(shiny)
library(ggplot2)
library(pROC)
library(PRROC)
library(dplyr)
library(DT)

# Simulate evaluation results externally and load here
load("eval_results.RData")

ui <- fluidPage(
  titlePanel("Model Evaluation Dashboard"),
  sidebarLayout(
    sidebarPanel(
      selectInput("model_select", "Choose Model:", choices = names(eval_results)),
      actionButton("refresh", "Refresh Evaluation Results"),
      br(),
      downloadButton("download_metrics", "Download Metrics Summary"),
      downloadButton("download_roc", "Download ROC Plot"),
      downloadButton("download_pr", "Download PR Plot"),
      br(),
      helpText("Metrics Definitions:"),
      tags$ul(
        tags$li("Accuracy: Overall correctness"),
        tags$li("Sensitivity: True Positive Rate"),
        tags$li("Specificity: True Negative Rate"),
        tags$li("Precision: Positive Predictive Value"),
        tags$li("F1 Score: Harmonic mean of Precision and Sensitivity")
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Metrics Summary", DTOutput("metrics_table")),
        tabPanel("ROC Curve", plotOutput("roc_plot")),
        #tabPanel("Precision-Recall Curve", plotOutput("pr_plot")),
        tabPanel("Confusion Matrix", verbatimTextOutput("conf_matrix"))
      )
    )
  )
)

server <- function(input, output, session) {
  observeEvent(input$refresh, {
    load("eval_results.RData", envir = globalenv())
    updateSelectInput(session, "model_select", choices = names(eval_results))
  })
  
  output$metrics_table <- renderDT({
    get_metrics_summary(eval_results)
  })
  
  output$roc_plot <- renderPlot({
    plot_combined_roc_ggplot(eval_results)
  })
  
  output$pr_plot <- renderPlot({
    pr_data <- data.frame()
    for (model in names(eval_results)) {
      preds <- eval_results[[model]]$Predictions
      probs <- preds$Prob_Non_Normal
      actuals <- ifelse(preds$True_Class == "Non_Normal", 1, 0)
      pr_obj <- pr.curve(scores.class0 = probs[actuals == 1],
                         scores.class1 = probs[actuals == 0], curve = TRUE)
      pr_points <- data.frame(
        Recall = pr_obj$curve[, 1],
        Precision = pr_obj$curve[, 2],
        Model = model
      )
      pr_data <- rbind(pr_data, pr_points)
    }
    # uncomment to plot precision graph
    # ggplot(pr_data, aes(x = Recall, y = Precision, color = Model)) +
    #   geom_line(size = 1.2) +
    #   labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
    #   theme_minimal(base_size = 14) +
    #   theme(legend.title = element_blank()) +
    #   scale_color_brewer(palette = "Dark2")
  })
  
  output$conf_matrix <- renderPrint({
    eval_results[[input$model_select]]$ConfusionMatrix
  })
  
  output$download_metrics <- downloadHandler(
    filename = function() { "metrics_summary.csv" },
    content = function(file) {
      write.csv(get_metrics_summary(eval_results), file, row.names = FALSE)
    }
  )
  
  output$download_roc <- downloadHandler(
    filename = function() { "roc_curve.png" },
    content = function(file) {
      png(file, width = 800, height = 600)
      print(plot_combined_roc_ggplot(eval_results))
      dev.off()
    }
  )
  
  output$download_pr <- downloadHandler(
    filename = function() { "pr_curve.png" },
    content = function(file) {
      png(file, width = 800, height = 600)
      pr_data <- data.frame()
      for (model in names(eval_results)) {
        preds <- eval_results[[model]]$Predictions
        probs <- preds$Prob_Non_Normal
        actuals <- ifelse(preds$True_Class == "Non_Normal", 1, 0)
        pr_obj <- pr.curve(scores.class0 = probs[actuals == 1],
                           scores.class1 = probs[actuals == 0], curve = TRUE)
        pr_points <- data.frame(
          Recall = pr_obj$curve[, 1],
          Precision = pr_obj$curve[, 2],
          Model = model
        )
        pr_data <- rbind(pr_data, pr_points)
      }
      print(
        ggplot(pr_data, aes(x = Recall, y = Precision, color = Model)) +
          geom_line(size = 1.2) +
          labs(title = "Precision-Recall Curve", x = "Recall", y = "Precision") +
          theme_minimal(base_size = 14) +
          theme(legend.title = element_blank()) +
          scale_color_brewer(palette = "Dark2")
      )
      dev.off()
    }
  )
}

shinyApp(ui, server)
