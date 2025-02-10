server <- function(input, output) {
  
  data <- reactive({
    if(input$data_source == "gen"){
      Generate_data(datagen.type = 1, 
                    n = input$n, 
                    dist = input$dist,
                    two_samples = (input$test_type == "two"))
    } else {
      req(input$file)
      Generate_data(datagen.type = 2,
                    two_samples = (input$test_type == "two"))
    }
  })
  
  output$data_plot <- renderPlot({
    df <- data()
    if(input$test_type == "one"){
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
    if(input$test_type == "one"){
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
    if(input$test_type == "one"){
      generate_tests(df$x, input$test)
    } else {
      generate_tests(df$x - df$y, input$test)
    }
  })
  
  output$test_output <- renderPrint({
    result <- test_result()
    cat("Test:", input$test, "\n")
    cat("Test Statistic:", result$statistic, "\n")
    cat("P-value:", result$p.value, "\n")
    cat("Conclusion:", ifelse(result$p.value < input$alpha, 
                              "Reject Null Hypothesis", 
                              "Fail to Reject Null Hypothesis"))
  })
  
  output$power_plot <- renderPlot({
    df <- data()
    sample_sizes <- seq(10, 100, by = 10)
    power_values <- Calculate_power(alpha = input$alpha, N = 1000, 
                                    twosamples = (input$test_type == "two"),
                                    dist = input$dist, 
                                    sample_size = sample_sizes, 
                                    test = input$test,
                                    effect_size = 0.5, B = 100)
    
    ggplot(data.frame(SampleSize = sample_sizes, Power = power_values), 
           aes(SampleSize, Power)) +
      geom_line(color = "steelblue") +
      geom_point(color = "darkblue") +
      ggtitle("Power Analysis") +
      theme_minimal()
  })
  
  output$bootstrap_plot <- renderPlot({
    df <- data()
    sample_sizes <- seq(10, 100, by = 10)
    boot_results <- if(input$test_type == "two"){
      bootstrap_two_sample_power(df$x, df$y, effect_size = 0.5, 
                                 alpha = input$alpha, n_bootstrap = 1000,
                                 sample_size = sample_sizes)
    } else {
      bootstrap_one_sample_power(df$x, effect_size = 0.5, 
                                 alpha = input$alpha, n_bootstrap = 1000,
                                 sample_size = sample_sizes)
    }
    
    ggplot(data.frame(SampleSize = sample_sizes, Power = boot_results), 
           aes(SampleSize, Power)) +
      geom_line(color = "firebrick") +
      geom_point(color = "darkred") +
      ggtitle("Bootstrap Power Analysis") +
      theme_minimal()
  })
}

shinyApp(ui = ui, server = server)
