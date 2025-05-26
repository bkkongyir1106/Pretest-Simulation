# server.R
function(input, output, session) {
  # Create a reactive trigger for the Run button
  run_trigger <- reactiveVal(0)
  
  observeEvent(input$run, {
    run_trigger(run_trigger() + 1)
  })
  
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
  
  # Create density plots for sample data
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
  # Generate summary statistics for sample data
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
  # function for performing location test
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
          # Source custom functions if needed
          if (test == "custom") {
            req(input$custom_test_file)
            source(input$custom_test_file$datapath, local = TRUE)
            if (twosamples) {
              if (!exists("custom_two_sample")) stop("Custom two-sample function not found.")
            } else {
              if (!exists("custom_one_sample")) stop("Custom one-sample function not found.")
            }
          }
          
          for (i in 1:N) {
            if (twosamples) {
              data_sim <- Generate_data(datagen.type = 1, n = sample_size, dist = input$dist, two_samples = TRUE)
              if (test == "custom") {
                pval <- custom_two_sample(data_sim$x, data_sim$y, effect_size)
              } else {
                pval <- TwoSample.test(data_sim$x, data_sim$y, effect_size, test, alpha, B)
              }
            } else {
              data_sim <- Generate_data(datagen.type = 1, n = sample_size, dist = input$dist, two_samples = FALSE)
              if (test == "custom") {
                pval <- custom_one_sample(data_sim$x,  effect_size)
              } else {
                pval <- OneSample.test(data_sim$x, effect_size, test, alpha, B)
              }
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
            if (test == "custom") {
              req(input$custom_test_file)
              source(input$custom_test_file$datapath, local = TRUE)
              if (!exists("custom_two_sample")) stop("Custom two-sample function not found.")
              boot_results <- replicate(n_bootstrap, {
                x_boot <- sample(df$x, replace = TRUE)
                y_boot <- sample(df$y, replace = TRUE)
                custom_two_sample(x_boot, y_boot, effect_size)
              })
            } else {
              boot_results <- bootstrap_two_sample(df$x, df$y, effect_size, alpha, n_bootstrap, test)
            }
          } else {
            if (test == "custom") {
              req(input$custom_test_file)
              source(input$custom_test_file$datapath, local = TRUE)
              if (!exists("custom_one_sample")) stop("Custom one-sample function not found.")
              boot_results <- replicate(n_bootstrap, {
                x_boot <- sample(df$x, replace = TRUE)
                custom_one_sample(x_boot, effect_size)
              })
            } else {
              boot_results <- bootstrap_one_sample(df$x, effect_size, alpha, n_bootstrap, test)
            }
          }
          power <- mean(boot_results < alpha)
          analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
          incProgress(1, detail = "Complete")
          list(power = power, analysis_type = analysis_type)
        })
      }
    }
  })
  
  # Produce Results
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
    run_trigger()  # Trigger when button is clicked
    isolate({
    req(input$test_type)
    if (input$test_type == "location") {
      effect_size <- input$effect_size
      analysis_type <- ifelse(effect_size == 0, "Type I Error", "Power")
      sample_sizes <- sample_sizes_power
      
      withProgress(message = "Generating power analysis...", value = 0, {
        tests <- test_names
        # Include custom test if selected
        if (input$loc_test == "custom" && !is.null(input$custom_test_file)) {
          tests <- c(tests, "custom")
          source(input$custom_test_file$datapath, local = TRUE)
        }
        n_tests <- length(tests)
        results <- purrr::map_dfr(1:n_tests, function(i) {
          test <- tests[i]
          incProgress(1/n_tests, detail = paste("Processing", test))
          
          # Check for custom functions
          if (test == "custom") {
            if (input$samples == "two") {
              req(exists("custom_two_sample"))
            } else {
              req(exists("custom_one_sample"))
            }
          }
          
          
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
              B = input$B,
              custom_func_path = if (test == "custom") input$custom_test_file$datapath else NULL
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
  })
  
  output$power_table <- renderTable({
    run_trigger()  # Trigger when button is clicked
    isolate({
     
    })
    
  })
}
