# server.R
function(input, output, session) {
  
  data <- reactive({
    if (input$data_source == "gen") {
      Generate_data(
        datagen.type = 1,
        n = input$n,
        dist = input$dist,
        two_samples = (input$samples == "two")
      )
    } else {
      req(input$file)
      Generate_data(
        datagen.type = 2,
        file_path = input$file$datapath,
        two_samples = (input$samples == "two")
      )
    }
  })
  
  output$data_plot <- renderPlot({
    df <- data()
    if (input$test_type == "normality" || input$samples == "one") {
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
      ggplot(df_long, aes(x = value, fill = group)) +
        geom_histogram(position = "identity", alpha = 0.6, bins = 30) +
        ggtitle("Two-Sample Distributions") +
        theme_minimal()
    }
  })
  
  output$summary_stats <- renderTable({
    df <- data()
    if (input$test_type == "normality" || input$samples == "one") {
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
    if (input$test_type == "normality") {
      df <- data()
      tests <- input$norm_tests
      results <- lapply(tests, function(test) generate_tests(df$x, test))
      names(results) <- tests
      results
    } else {
      df <- data()
      effect_size <- input$effect_size
      alpha <- input$alpha
      test <- input$loc_test
      twosamples <- (input$samples == "two")
      
      if (input$data_source == "gen") {
        withProgress(message = "Running simulations...", value = 0, {
          N <- input$N
          B <- input$B
          sample_size <- if(twosamples) input$n else input$n
          rejections <- numeric(N)
          
          for (i in 1:N) {
            if (twosamples) {
              data_sim <- Generate_data(datagen.type = 1, n = sample_size, dist = input$dist, two_samples = TRUE)
              pval <- TwoSample.test(data_sim$x, data_sim$y + effect_size, test, alpha, B)
            } else {
              data_sim <- Generate_data(datagen.type = 1, n = sample_size, dist = input$dist, two_samples = FALSE)
              pval <- OneSample.test(data_sim$x + effect_size, test, alpha, B)
            }
            rejections[i] <- (pval < alpha)
            incProgress(1/N, detail = paste("Iteration", i, "of", N))
          }
          power <- mean(rejections)
          analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
          list(power = power, analysis_type = analysis_type)
        })
      } else {
        withProgress(message = "Bootstrapping...", value = 0, {
          n_bootstrap <- input$n_bootstrap
          if (twosamples) {
            boot_results <- bootstrap_two_sample(df$x, df$y, effect_size, alpha, n_bootstrap, test)
          } else {
            boot_results <- bootstrap_one_sample(df$x, effect_size, alpha, n_bootstrap, test)
          }
          power <- mean(boot_results < alpha)
          analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
          incProgress(1, detail = "Complete")
          list(power = power, analysis_type = analysis_type)
        })
      }
    }
  })
  
  output$test_output <- renderPrint({
    if (input$test_type == "normality") {
      results <- test_result()
      cat("Normality Test Results:\n\n")
      for (test in names(results)) {
        cat("Test:", test, "\n")
        cat("Statistic:", results[[test]]$statistic, "\n")
        cat("P-value:", results[[test]]$p.value, "\n")
        conclusion <- ifelse(
          results[[test]]$p.value < input$alpha,
          "Reject Null Hypothesis",
          "Fail to Reject Null Hypothesis"
        )
        cat("Conclusion:", conclusion, "\n\n")
      }
    } else {
      result <- test_result()
      cat("Location Test:", input$loc_test, "\n")
      cat("Analysis Type:", result$analysis_type, "\n")
      cat(result$analysis_type, ":", round(result$power, 4), "\n")
    }
  })
  
  output$power_plot <- renderPlot({
    req(input$test_type)
    if (input$test_type == "location") {
      effect_size <- input$effect_size
      analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
      sample_sizes <- sample_sizes_power
      
      withProgress(message = "Generating power analysis...", value = 0, {
        n_tests <- length(test_names)
        results <- purrr::map_dfr(1:n_tests, function(i) {
          test <- test_names[i]
          incProgress(1/n_tests, detail = paste("Processing", test))
          
          if (input$data_source == "upload") {
            df <- data()
            if (input$samples == "one") {
              pwr <- bootstrap_one_sample(
                df$x, effect_size, input$alpha, input$n_bootstrap, sample_sizes, test = test
              )
            } else {
              pwr <- bootstrap_two_sample(
                df$x, df$y, effect_size, input$alpha, input$n_bootstrap, sample_sizes, test = test
              )
            }
          } else {
            pwr <- Calculate_power(
              alpha = input$alpha,
              N = input$N,
              twosamples = (input$samples == "two"),
              dist = input$dist,
              sample_size = sample_sizes,
              test = test,
              effect_size = effect_size,
              B = input$B
            )
          }
          data.frame(
            SampleSize = sample_sizes,
            Value = pwr,
            Test = test
          )
        })
      })
      
      ggplot(results, aes(SampleSize, Value, color = Test)) +
        geom_line() +
        geom_point() +
        labs(title = paste("Comparison of Location Tests -", analysis_type),
             x = "Sample Size",
             y = analysis_type) +
        theme_minimal() +
        scale_color_brewer(palette = "Set1") +
        ylim(0, ifelse(effect_size == 0, 0.5, 1))
    } else if (input$test_type == "normality") {
      req(input$norm_tests, input$data_source == "gen")
      alpha <- input$alpha
      sample_sizes <- sample_sizes_power
      dist <- input$dist
      N <- input$N_norm
      tests <- input$norm_tests
      
      withProgress(message = "Calculating power for normality tests...", value = 0, {
        n_tests <- length(tests)
        results <- purrr::map_dfr(tests, function(test) {
          incProgress(1/n_tests, detail = paste("Processing", test))
          
          pwr_values <- sapply(sample_sizes, function(n) {
            pvals <- replicate(N, {
              data_sim <- Generate_data(datagen.type = 1, n = n, dist = dist, two_samples = FALSE)
              generate_tests(data_sim$x, test)$p.value
            })
            mean(pvals < alpha)
          })
          
          data.frame(
            SampleSize = sample_sizes,
            Power = pwr_values,
            Test = test
          )
        })
      })
      
      ggplot(results, aes(x = SampleSize, y = Power, color = Test)) +
        geom_line() +
        geom_point() +
        labs(title = paste("Power of Normality Tests for", dist, "Distribution"),
             x = "Sample Size",
             y = ifelse(dist == "Normal", "Type I Error Rate", "Power")) +
        theme_minimal() +
        scale_color_brewer(palette = "Set1") +
        ylim(0, 1)
    }
  })
  
  output$power_table <- renderTable({
    # Similar modifications as in power_plot for table output if needed
  })
}