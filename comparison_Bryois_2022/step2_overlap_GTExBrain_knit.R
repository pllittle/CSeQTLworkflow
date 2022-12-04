
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools")

render_report = function(cs_q_cutoff, fold_cutoff) {
  config = sprintf("cs_q_%.0e_fold_%.1f", 
                   cs_q_cutoff, fold_cutoff)
  rmarkdown::render(
    "step2_overlap_GTExBrain.Rmd", params = list(
      cs_q_cutoff = cs_q_cutoff, 
      fold_cutoff = fold_cutoff
    ),
    output_file = sprintf("step2_overlap_GTExBrain_%s.html", config)
  )
}

render_report(cs_q_cutoff=5e-3, fold_cutoff=1)
render_report(cs_q_cutoff=5e-3, fold_cutoff=1.5)
render_report(cs_q_cutoff=1e-3, fold_cutoff=1.5)

sessionInfo()

q(save = "no")
