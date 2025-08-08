library(shiny)
library(DT)
library(caret)
library(pROC)

source("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/gen_data_fun.R")

if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest,
               evd, discretization, nnet, ggplot2, mlbench, infotheo, dplyr)

# Load pre-trained models and objects
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/trained_models.RData")
# ---------------------------
# Define Feature Extraction Functions
# ---------------------------
calculate_zero_crossing_rate <- function(samples) {
  signs <- samples > 0
  zero_crossings <- sum(abs(diff(signs)))
  return(zero_crossings / (length(samples) - 1))
}

calculate_gini_coefficient <- function(samples) {
  samples_abs <- abs(samples - min(samples))  
  return(ineq(samples_abs, type = "Gini"))
}

calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  return(sum(samples < (qnt[1] - H) | samples > (qnt[2] + H)))
}

calculate_entropy <- function(samples) {
  return(infotheo::entropy(discretize(samples, nbins = 10)))
}

calculate_peak_to_trough <- function(samples) max(samples) / abs(min(samples))

calculate_box_test  <- function(samples) {
  as.numeric(Box.test(samples, lag = 1, type = "Ljung-Box")$statistic)
}


calculate_spectral_entropy <- function(samples) {
  spec <- stats::spec.pgram(samples, plot = FALSE)
  p <- spec$spec / sum(spec$spec)
  -sum(p * log(p))
}

calculate_spectral_centroid <- function(samples) {
  spec <- stats::spec.pgram(samples, plot = FALSE)
  p <- spec$spec / sum(spec$spec)
  freq <- spec$freq
  sum(freq * p) / sum(p)
}

# Fractal & complexity measures
calculate_hjorth <- function(samples) {
  activity   <- var(samples)
  mobility   <- sqrt(var(diff(samples)) / activity)
  complexity <- (sqrt(var(diff(diff(samples))) / var(diff(samples)))) / mobility
  c(Activity = activity, Mobility = mobility, Complexity = complexity)
}


calculate_fractal_dimension <- function(samples) {
  res <- fractaldim::fd.estimate(samples, methods = "madogram")
  as.numeric(res$fd)
}

calculate_energy <- function(samples) {
  return(sum(samples^2))
}

# Dispersion & variation
calculate_cv             <- function(samples) sd(samples) / mean(samples)
calculate_range          <- function(samples) max(samples) - min(samples)
calculate_entropy        <- function(samples) {
  infotheo::entropy(discretize(samples, nbins = 10))
}

# ---------------------------
# extract all features
# ---------------------------
calculate_features <- function(samples) {
  # Central & spread
  mean_val        <- mean(samples)
  median_val      <- median(samples)
  var_val         <- var(samples)
  iqr_val         <- IQR(samples)
  mad_val         <- mad(samples)
  range_val       <- max(samples) - min(samples)
  cv_val          <- calculate_cv(samples)
  rms_val         <- sqrt(mean(samples^2))
  
  # Normality & shape
  skewness        <- e1071::skewness(samples)
  kurtosis        <- e1071::kurtosis(samples)
  jb_stat         <- as.numeric(tseries::jarque.bera.test(samples)$statistic)
  ad_stat         <- as.numeric(nortest::ad.test(samples)$statistic)
  sw_stat         <- shapiro.test(samples)$statistic
  sf_stat         <- nortest::sf.test(samples)$statistic
  lf_stat         <- nortest::lillie.test(samples)$statistic
  cvm_stat        <- nortest::cvm.test(samples)$statistic
  
  # Time & distributional
  zcr             <- calculate_zero_crossing_rate(samples)
  gini            <- calculate_gini_coefficient(samples)
  outliers        <- calculate_outliers(samples)
  
  # Others
  entropy_val     <- calculate_entropy(samples)
  pt_ratio        <- calculate_peak_to_trough(samples)
  box_val         <- calculate_box_test(samples)
  spec_entropy    <- calculate_spectral_entropy(samples)
  spec_centroid   <- calculate_spectral_centroid(samples)
  fd_val          <- calculate_fractal_dimension(samples)
  hjorth_vals     <- calculate_hjorth(samples)
  energy          <- calculate_energy(samples)
  
  # create a dataframe
  features <- data.frame(
    Skewness                = skewness,
    Kurtosis                = kurtosis,
    Jarque_Bera             = jb_stat,
    Anderson_Darling        = ad_stat,
    Shapiro_Wilk            = sw_stat,
    Shapiro_Francia         = sf_stat,
    Lilliefors              = lf_stat,
    Cramer_Von_Misse        = cvm_stat,
    Zero_Cross_Rate         = zcr,
    Gini_Coefficient        = gini,
    #Outliers                = outliers,
    Mean                     = mean_val,
    Median                   = median_val,
    Variance                 = var_val,
    IQR                      = iqr_val,
    MAD                      = mad_val,
    Range                    = range_val,
    CV                       = cv_val,
    Root_Mean_Square         = rms_val,
    energy                   = energy,
    Peak_to_Trough           = pt_ratio,
    Enropy                   = entropy_val,
    Spectral_Entropy         = spec_entropy,
    Spectral_Centroid        = spec_centroid,
    Box_Ljung_Stat           = box_val,
    Fractal_Dimension        = fd_val,
    Hjorth_Activity          = hjorth_vals["Activity"],
    Hjorth_Mobility          = hjorth_vals["Mobility"],
    Hjorth_Complexity        = hjorth_vals["Complexity"]
  )
  return(features)
}


# ---------------------------
# Standardization & Normalization Function
# ---------------------------
preprocess_data <- function(train_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  
  # Step 1: Standardization (Z-score scaling: mean = 0, std = 1)
  preProcStandard <- preProcess(numeric_train, method = c("center", "scale"))
  train_std <- predict(preProcStandard, numeric_train)
  
  # Step 2: Normalization (Rescale to [0,1])
  preProcNorm <- preProcess(train_std, method = "range")
  train_norm <- predict(preProcNorm, train_std)
  
  train_norm$Label <- train_data$Label
  
  return(list(train = train_norm, preProcStandard = preProcStandard, preProcNorm = preProcNorm))
}

# ---------------------------
# Data Generation Function 
# ---------------------------
generate_data <- function(sample_size, N, dist = "normal", label) {
  data <- do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(sample_size, dist) 
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}

# ---------------------------
# Classification Function
# ---------------------------
classify_sample <- function(sample_data, trained_models, preProcStandard, preProcNorm) {
  sample_features <- calculate_features(sample_data)
  sample_std <- predict(preProcStandard, sample_features)
  sample_norm <- predict(preProcNorm, sample_std)
  
  results <- lapply(names(trained_models), function(model_name) {
    model <- trained_models[[model_name]]
    pred_class <- as.character(predict(model, newdata = sample_norm))
    pred_prob <- predict(model, newdata = sample_norm, type = "prob")$Non_Normal
    
    data.frame(
      Model = model_name,
      Prediction = pred_class,
      Non_Normal_Probability = round(pred_prob, 4),
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, results)
}

# ---------------------------
# Base R Plot Functions
# ---------------------------
plot_var_importance <- function(rf_model, top_n = 30) {
  rf_var_imp <- varImp(rf_model)
  imp_df <- data.frame(Feature = rownames(rf_var_imp$importance), 
                       Importance = rf_var_imp$importance$Overall)
  imp_df <- imp_df[order(imp_df$Importance), ]
  if (nrow(imp_df) > top_n) imp_df <- tail(imp_df, top_n)
  
  par(mar = c(5, 8, 4, 2))
  barplot(imp_df$Importance, 
          names.arg = imp_df$Feature,
          horiz = TRUE,
          las = 1,
          col = "steelblue",
          main = "Top Feature Importance (Random Forest)",
          xlab = "Importance Score")
}

plot_combined_roc_base <- function(eval_results) {
  colors <- rainbow(length(eval_results))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = "Combined ROC Curves")
  abline(a = 0, b = 1, col = "gray", lty = 2)
  
  legend_labels <- c()
  for (i in seq_along(names(eval_results))) {
    model_name <- names(eval_results)[i]
    pred_df <- eval_results[[model_name]]$Predictions
    
    actual_bin <- ifelse(pred_df$True_Class == "Non_Normal", 1, 0)
    probs <- pred_df$Prob_Non_Normal
    
    roc_obj <- roc(actual_bin, probs, quiet = TRUE)
    lines(1 - roc_obj$specificities, roc_obj$sensitivities, 
          col = colors[i], lwd = 2)
    
    auc_val <- round(auc(roc_obj), 3)
    legend_labels <- c(legend_labels, paste0(model_name, " (AUC=", auc_val, ")"))
  }
  
  legend("bottomright", legend = legend_labels, 
         col = colors, lty = 1, lwd = 2, cex = 0.8)
  grid()
}

plot_input_data <- function(samples) {
  par(mfrow = c(ceiling(sqrt(length(samples))), ceiling(sqrt(length(samples)))))
  for (i in seq_along(samples)) {
    plot(samples[[i]], type = "b", col = "blue", pch = 19,
         main = ifelse(length(samples) > 1, paste("Sample", i), "Input Data"),
         xlab = "Index", ylab = "Value")
  }
}

# ---------------------------
# Shiny App UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Machine Learning Approach to Normality Test"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Input Data"),
      radioButtons("data_source", "Data Source:",
                   choices = c("Upload CSV" = "upload", 
                               "Generate Data" = "generate"),
                   selected = "upload"),
      
      conditionalPanel(
        condition = "input.data_source == 'upload'",
        fileInput("file", "Choose CSV File",
                  accept = c("text/csv", ".csv")),
        helpText("File format: Single column with numeric values")
      ),
      
      conditionalPanel(
        condition = "input.data_source == 'generate'",
        selectInput("dist", "Distribution:",
                    choices = c("Normal" = "normal",
                                "LogNormal" = "LogNormal",
                                "Chi_Square" = "Chi_Square",
                                "Laplace"    = "Laplace",
                                "beta"       = "beta",
                                "Exponential" = "Exponential",
                                "Weibull" = "Weibull",
                                "Uniform" = "Uniform")),
        numericInput("sample_size", "Sample Size:", value = 8, min = 8),
        numericInput("n_samples", "Number of Samples:", value = 1, min = 1)
      ),
      
      actionButton("run", "Run Classification", class = "btn-primary"),
      br(), br(),
      downloadButton("download", "Download Results")
    ),
    
    mainPanel(
      tabsetPanel(
        # tabPanel("Input Data", 
        #          plotOutput("data_plot"),
        #          dataTableOutput("data_table")),
        
        tabPanel("Predictions", 
                 dataTableOutput("pred_table")),
        
        tabPanel("Feature Importance", 
                 plotOutput("varimp_plot" , width = "700px", height = "800px")),
        
        tabPanel("Model Performance",
                 plotOutput("roc_plot" , width = "800px", height = "600px"),
                 dataTableOutput("metrics_table") )
      )
    )
  )
)

# ---------------------------
# Shiny App Server
# ---------------------------
server <- function(input, output) {
  
  # Reactive data input
  sample_data <- reactive({
    if (input$data_source == "upload") {
      req(input$file)
      df <- read.csv(input$file$datapath, header = TRUE)
      if (ncol(df) != 1) {
        showNotification("Error: File should have exactly one column", type = "error")
        return(NULL)
      }
      return(df[, 1])
    } else {
      set.seed(123)
      do.call(c, lapply(1:input$n_samples, function(i) {
        generate_samples(input$sample_size, input$dist)
      }))
    }
  })
  
  # Processed samples
  processed_samples <- reactive({
    data <- sample_data()
    req(data)
    samples <- split(data, ceiling(seq_along(data)/8))
    samples <- samples[lengths(samples) == 8]  # Keep only complete samples
    samples
  })
  
  # Classification results
  classification_results <- eventReactive(input$run, {
    req(processed_samples())
    samples <- processed_samples()
    
    results <- lapply(seq_along(samples), function(i) {
      res <- classify_sample(
        samples[[i]],
        models_list,
        norm_result$preProcStandard,
        norm_result$preProcNorm
      )
      res$SampleID <- i
      res
    })
    
    do.call(rbind, results)
  })
  
  # Output: Input data plot (base R)
  output$data_plot <- renderPlot({
    samples <- processed_samples()
    req(samples)
    plot_input_data(samples)
  })
  
  # Output: Input data table
  output$data_table <- renderDataTable({
    samples <- processed_samples()
    req(samples)
    sample_matrix <- do.call(rbind, samples)
    colnames(sample_matrix) <- paste0("X", 1:8)
    datatable(sample_matrix, 
              options = list(pageLength = 5, scrollX = TRUE),
              caption = "Input Data Samples")
  })
  
  # Output: Prediction table
  output$pred_table <- renderDataTable({
    results <- classification_results()
    req(results)
    datatable(results, 
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Classification Results") %>%
      formatStyle(
        "Prediction",
        backgroundColor = styleEqual(c("Normal", "Non_Normal"), c("#E6F5E6", "#FFE6E6"))
      )
  })
  
  # Output: Variable importance plot 
  output$varimp_plot <- renderPlot({
    plot_var_importance(models_list$RF)
  })
  
  # Output: ROC curve 
  output$roc_plot <- renderPlot({
    plot_combined_roc_base(eval_results)
  })
  
  # Output: Performance metrics table
  output$metrics_table <- renderDataTable({
    metrics <- lapply(names(eval_results), function(model) {
      data.frame(
        Model = model,
        Accuracy = round(eval_results[[model]]$Accuracy, 4),
        Sensitivity = round(eval_results[[model]]$Sensitivity, 4),
        Specificity = round(eval_results[[model]]$Specificity, 4),
        F1_Score = round(eval_results[[model]]$F1, 4)
      )
    })
    
    datatable(do.call(rbind, metrics), 
              options = list(dom = 't', pageLength = 6),
              caption = "Model Performance Metrics",
              rownames = FALSE)
  })
  
  # Download handler
  output$download <- downloadHandler(
    filename = function() {
      paste("classification-results-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(classification_results(), file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)