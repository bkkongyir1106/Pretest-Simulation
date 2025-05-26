# ============================
# MODEL EVALUATION DASHBOARD
# ============================

# Load required libraries
library(shiny)
library(ggplot2)
library(pROC)
library(dplyr)

# === UI ===
ui <- fluidPage(
  titlePanel("Model Evaluation Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_model", "Select Model:", choices = NULL),
      actionButton("refresh", "Refresh Models")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("ROC Curve", plotOutput("rocPlot")),
        #tabPanel("Precision-Recall Curve", plotOutput("prPlot")),
        tabPanel("Metrics Summary", tableOutput("metricsTable"))
      )
    )
  )
)

# === SERVER ===
server <- function(input, output, session) {
  
  # Load evaluation results (assumes pre-computed in environment)
  eval_results <- reactiveVal()
  
  observeEvent(input$refresh, {
    # If eval_results already exists in environment, load it
    if (exists("eval_results", envir = .GlobalEnv)) {
      eval_results(get("eval_results", envir = .GlobalEnv))
      updateSelectInput(session, "selected_model", choices = names(eval_results()))
    }
  })
  
  output$rocPlot <- renderPlot({
    req(eval_results())
    
    roc_data <- data.frame()
    
    for (model_name in names(eval_results())) {
      pred_df <- eval_results()[[model_name]]
      
      if (!"Prob_Non_Normal" %in% colnames(pred_df$Predictions)) next
      
      actual_bin <- ifelse(pred_df$Predictions$True_Class == "Non_Normal", 1, 0)
      probs <- pred_df$Predictions$Prob_Non_Normal
      
      roc_obj <- roc(actual_bin, probs, levels = c(0, 1), direction = "<", quiet = TRUE)
      auc_value <- round(auc(roc_obj), 3)
      
      roc_points <- data.frame(
        FPR = 1 - roc_obj$specificities,
        TPR = roc_obj$sensitivities,
        Model = paste0(model_name, " (AUC=", auc_value, ")")
      )
      
      roc_data <- rbind(roc_data, roc_points)
    }
    
    ggplot(roc_data, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      geom_abline(linetype = "dashed", color = "gray") +
      labs(
        title = "ROC Curves (Positive Class: Non_Normal)",
        x = "False Positive Rate",
        y = "True Positive Rate"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.title = element_blank())
  })
  
  output$prPlot <- renderPlot({
    req(eval_results())
    
    pr_data <- data.frame()
    
    for (model_name in names(eval_results())) {
      pred_df <- eval_results()[[model_name]]
      
      if (!"Prob_Non_Normal" %in% colnames(pred_df$Predictions)) next
      
      actual_bin <- ifelse(pred_df$Predictions$True_Class == "Non_Normal", 1, 0)
      probs <- pred_df$Predictions$Prob_Non_Normal
      
      precision <- c()
      recall <- c()
      thresholds <- sort(unique(probs), decreasing = TRUE)
      
      for (t in thresholds) {
        pred_class <- ifelse(probs >= t, 1, 0)
        tp <- sum(pred_class == 1 & actual_bin == 1)
        fp <- sum(pred_class == 1 & actual_bin == 0)
        fn <- sum(pred_class == 0 & actual_bin == 1)
        
        prec <- ifelse(tp + fp == 0, 1, tp / (tp + fp))
        rec <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
        
        precision <- c(precision, prec)
        recall <- c(recall, rec)
      }
      
      pr_points <- data.frame(
        Precision = precision,
        Recall = recall,
        Model = model_name
      )
      
      pr_data <- rbind(pr_data, pr_points)
    }
    # uncomment to plot precision grap
    # ggplot(pr_data, aes(x = Recall, y = Precision, color = Model, linetype = Model)) +
    #   geom_line(size = 1.2) +
    #   labs(
    #     title = "Precision-Recall Curves",
    #     x = "Recall",
    #     y = "Precision"
    #   ) +
    #   theme_minimal(base_size = 14) +
    #   theme(legend.title = element_blank())
  })
  
  output$metricsTable <- renderTable({
    req(eval_results())
    summary_df <- data.frame()
    
    for (model_name in names(eval_results())) {
      metrics <- eval_results()[[model_name]]
      
      summary_df <- rbind(summary_df, data.frame(
        Model = model_name,
        Accuracy = round(metrics$Accuracy, 4),
        Sensitivity = round(metrics$Sensitivity, 4),
        Specificity = round(metrics$Specificity, 4),
        Precision = round(metrics$Precision, 4),
        F1 = round(metrics$F1, 4)
      ))
    }
    summary_df
  })
}

# === RUN APP ===
shinyApp(ui = ui, server = server)

