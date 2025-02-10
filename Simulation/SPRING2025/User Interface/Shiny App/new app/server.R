
function(input, output) {
  
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
    analysis_type <- if(input$test_type == "normality"){
      cat("Analysis Type:", analysis_type, "\n")
      conclusion <- ifelse(
        if(input$test_type == "normality") result$p.value else result < input$alpha,
        "Reject Null Hypothesis", 
        "Fail to Reject Null Hypothesis"
      )
      cat("Conclusion:", conclusion)
    }else{
      if(input$test_type == "loc_test"){
      analysis_type <- ifelse(input$effect_size == 0, "Type I Error", "Power")
      cat("Analysis Type:", analysis_type, "\n")
      conclusion <- ifelse(
        if(input$test_type == "loc_test") result$p.value else result < input$alpha,
        "Reject Null Hypothesis", 
        "Fail to Reject Null Hypothesis"
      )
      cat("Conclusion:", conclusion)
      }
    }
  })
  
  output$power_plot <- renderPlot({
    req(input$test_type == "location")
    
    df <- data()
    effect_size <- input$effect_size
    analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
    
    if(input$data_source == "upload") {
      # Bootstrap analysis for uploaded data
      results <- map_dfr(test_names, ~{
        if(input$samples == "one") {
          pwr <- bootstrap_one_sample(df$x, effect_size, input$alpha, 
                                      n_bootstrap = 100, sample_size = sample_sizes_power)
        } else {
          pwr <- bootstrap_two_sample(df$x, df$y, effect_size, input$alpha,
                                      n_bootstrap = 100, sample_size = sample_sizes_power)
        }
        data.frame(
          SampleSize = sample_sizes_power,
          Value = pwr,
          Test = .x
        )
      })
    } else {
      # Simulation-based analysis for generated data
      results <- map_dfr(test_names, ~{
        pwr <- Calculate_power(
          alpha = input$alpha, 
          N = 100, 
          twosamples = (input$samples == "two"),
          dist = input$dist, 
          sample_size = sample_sizes_power, 
          test = .x,
          effect_size = effect_size, 
          B = 50
        )
        data.frame(
          SampleSize = sample_sizes_power,
          Value = pwr,
          Test = .x
        )
      })
    }
    
    ggplot(results, aes(SampleSize, Value, color = Test)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Comparison of Location Tests -", analysis_type),
           x = "Sample Size",
           y = analysis_type) +
      theme_minimal() +
      scale_color_brewer(palette = "Set1") +
      if(effect_size == 0){
        ylim(0, 0.5)
      }else{
        ylim(0, 1)
      }
      
  })
}