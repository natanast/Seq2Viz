# app.R
library(shiny)

# Load and run the actual app from the R/ folder
source("R/main.R")

shinyApp(ui, server)