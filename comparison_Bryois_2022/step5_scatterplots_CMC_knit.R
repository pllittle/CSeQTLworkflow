
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools")

render_report = function(AF) {
  rmarkdown::render(
    "step5_scatterplots_CMC.Rmd", params = list(
      AF = AF
    ),
    output_file = sprintf("step5_scatterplots_CMC_%s.html", AF)
  )
}

render_report(AF="SCZ")
render_report(AF="Control")

sessionInfo()

q(save = "no")
