
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools")

render_report = function(AF) {
  rmarkdown::render(
    "step9b_summarize_results.Rmd", params = list(
      AF = AF
    ),
    output_file = sprintf("step9b_summarize_results_%s.html", AF)
  )
}

render_report(AF="CMC_SCZ")
render_report(AF="CMC_Control")
render_report(AF="GTExBrain")

sessionInfo()

q(save = "no")
