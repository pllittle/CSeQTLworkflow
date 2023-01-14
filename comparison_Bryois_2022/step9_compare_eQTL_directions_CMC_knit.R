
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools")

render_report = function(AF) {
  rmarkdown::render(
    "step9_compare_eQTL_directions_CMC.Rmd", params = list(
      AF = AF
    ),
    output_file = sprintf("step9_compare_eQTL_directions_CMC_%s.html", AF)
  )
}

# render_report(AF="SCZ")
render_report(AF="Control")

sessionInfo()

q(save = "no")
