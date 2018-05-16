library(shiny)


shinyUI(fluidPage(
  pageWithSidebar(
    titlePanel("Samples Sheet Validator"),
    sidebarPanel(width = 4,
      fileInput("samples_sheet_file", 
                 label = "Sample sheet file (csv)",
                 multiple = FALSE),
      selectInput('seq_tech',
                  label = "Sequencing technology",
                  choices = c('HiSeq' = 'HiSeq',
                              'MiSeq' = 'MiSeq'),
                  selected = 'HiSeq'),
      h5('Run the validation'),
      actionButton("runBtn", "Validate Samples Sheet"),
      hr(),
      h6('For problems, contact Nick Youngblut (nyoungblut@tuebingen.mpg.de)')
    ),
    mainPanel(
      h4('Validator script output'),
      verbatimTextOutput('script_out')
    )
  )
))