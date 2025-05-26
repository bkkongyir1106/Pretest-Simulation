install.packages("shiny", dependencies = TRUE)
# Load shiny package
library(shiny)

# Define UI
ui <- fluidPage(
  h1("GeeksforGeeks article on shiny package in R"),
  p(style = "font-family:Impact", "My first shiny app")
)

# Define Server
server <- function(input, output) {}

# Run the application 
shinyApp(ui = ui, server = server)
# -------------------------------

# define fluid page layout
ui <- fluidPage(
  sliderInput(inputId = "num", 
              label = "Choose a number", 
              value = 10, min = 1, max = 1000),
  plotOutput("hist")
)

server <- function(input, output)
{
  output$hist <- renderPlot({
    hist(rnorm(input$num))
  })
}

# create shiny app object 
# using shinyApp
shinyApp(ui = ui, server = server)

# -------------------------------
# reactive(): It creates a reactive expression.

# define fluid page layout
ui <- fluidPage(
  textInput(inputId = "num", label = "Choose a number", 
            value = "", width = 100, placeholder = NULL),
  plotOutput("hist"),
  verbatimTextOutput("stats")
)
server <- function(input, output) 
{
  # use reactive to create
  # a reactive expression
  data <- reactive({rnorm(input$num)})
  
  output$hist <- renderPlot({hist(data())})
  output$stats <- renderPrint({summary(data())})
}

# create shiny app object 
# using shinyApp
shinyApp(ui = ui, server = server)
# --------------------------------------