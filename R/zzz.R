# Global variables to avoid R CMD check notes
# These are column names used in dplyr/tidyr/ggplot2 pipelines
#' @importFrom stats coef complete.cases embed fitted lm logLik pchisq pf pnorm quantile residuals rnorm sd var vcov
#' @importFrom ggplot2 ggplot aes geom_bar geom_line geom_point geom_ribbon geom_hline geom_abline geom_tile geom_text geom_smooth geom_histogram geom_density stat_qq stat_qq_line facet_wrap labs theme theme_minimal element_text element_rect element_blank coord_fixed position_dodge scale_color_manual scale_color_brewer scale_color_viridis_d scale_fill_manual scale_fill_gradient2 scale_shape_manual after_stat
#' @importFrom tidyr pivot_longer
NULL

utils::globalVariables(c(
  # Column names used in plotting
  "value", "type", "estimate", "se", "tau", "variable", "coefficient",
  "multiplier", "horizon", "persistence", "positive", "negative", 
  "asymmetric", "density",
  

  # Other symbols
  ".", "index"
))
