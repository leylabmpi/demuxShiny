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
      h4('Run the validation'),
      actionButton("runBtn", "Validate")
    ),
    mainPanel(
      verbatimTextOutput('script_out')
    )
  )
))