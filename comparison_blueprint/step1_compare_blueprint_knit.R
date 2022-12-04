
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/quarto/bin/tools")

render_report = function(ot_q_cutoff, cs_q_cutoff, ot_fold_cutoff, 
                         cs_fold_cutoff) {
  config = sprintf("ot_cs_q_%.0e_%.0e_fold_%.1f_%.1f", 
                   ot_q_cutoff, cs_q_cutoff, ot_fold_cutoff, cs_fold_cutoff)
  config
  
  rmarkdown::render(
    "step1_compare_blueprint.Rmd", params = list(
      ot_q_cutoff = ot_q_cutoff, 
      cs_q_cutoff = cs_q_cutoff, 
      ot_fold_cutoff = ot_fold_cutoff,
      cs_fold_cutoff = cs_fold_cutoff
    ),
    output_file = sprintf("step1_compare_blueprint_%s.html", config)
  )
}


render_report(ot_q_cutoff=5e-5, cs_q_cutoff=5e-3, ot_fold_cutoff=1, 
              cs_fold_cutoff=1)

render_report(ot_q_cutoff=5e-5, cs_q_cutoff=5e-3, ot_fold_cutoff=1, 
              cs_fold_cutoff=1.5)

render_report(ot_q_cutoff=5e-5, cs_q_cutoff=1e-3, ot_fold_cutoff=1, 
              cs_fold_cutoff=1.5)

sessionInfo()

q(save = "no")
