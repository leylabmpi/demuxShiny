library(shiny)
source("../utils/io.R")
source("../utils/format.R")

call_script = function(script_path, input){
  if(is.null(input$samples_sheet_file)){
    options = c('-h')
  } else {
    options = c(
      rename_tmp_file(input$samples_sheet_file),
      c('--seq-tech', add_quotes(input$seq_tech))
    )
  }
  system2(script_path, options, stdout=TRUE, stderr=TRUE)
}

shinyServer(function(input, output, session) {
  script_path = './samples_sheet_validator.py'
  
  # calling script
  script_out = eventReactive(input$runBtn, {
    # run command 
    call_script(script_path, input)
  })
  
  # adding script output to output
  output$script_out = reactive({
    paste(script_out(), collapse='\n')
  })
})