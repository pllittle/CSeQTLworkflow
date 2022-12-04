
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools")

render_report = function(cs_q_cutoff, fold_cutoff, AF) {
  config = sprintf("%s_cs_q_%.0e_fold_%.1f", 
                   AF, cs_q_cutoff, fold_cutoff)
  config
  rmarkdown::render(
    "step4_overlap_CMC.Rmd", params = list(
      cs_q_cutoff = cs_q_cutoff, 
      fold_cutoff = fold_cutoff, 
      AF = AF
    ),
    output_file = sprintf("step4_overlap_CMC_%s.html", config)
  )
}

render_report(cs_q_cutoff=5e-3, fold_cutoff=1, AF="SCZ")
render_report(cs_q_cutoff=5e-3, fold_cutoff=1.5, AF="SCZ")
render_report(cs_q_cutoff=1e-3, fold_cutoff=1.5, AF="SCZ")

render_report(cs_q_cutoff=5e-3, fold_cutoff=1, AF="Control")
render_report(cs_q_cutoff=5e-3, fold_cutoff=1.5, AF="Control")
render_report(cs_q_cutoff=1e-3, fold_cutoff=1.5, AF="Control")

sessionInfo()

q(save = "no")
