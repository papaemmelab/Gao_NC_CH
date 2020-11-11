plot_forest <- function(class_results, x = "class", y = "estimate", ymin = "conf.low", ymax = "conf.high",
                       label = "p.label", limits = NULL, breaks = waiver(), title = "", col = NULL, fill = NULL,
                        dodge_width = 0.8, outer_limit_arrows = FALSE, ps=3, eb_w=0.4, eb_s=0.4, or_s=4, OR=T, yinter = 1,
                        nudge = 0, base_plot = geom_blank(), bar_col = 'black') { 
    # suppressWarnings(ggthemr::ggthemr("fresh"))
  #' Forest plot of the coefficients in 'class_results'
  output_plot <- ggplot(class_results, aes_string(x = x, y = y, ymin = ymin, ymax = ymax, col = col, fill = fill, label = label, vjust = nudge)) + 
    base_plot +
    geom_hline(yintercept = yinter, color = "gray", linetype = "solid") +
    geom_errorbar(position = position_dodge(width = dodge_width), width = eb_w, size=eb_s) + 
    geom_point(position = position_dodge(width = dodge_width), size=ps) +
    geom_text(position = position_dodge(width = dodge_width), size = or_s, alpha = .9) +
    coord_flip() +
    ggtitle(title)
  
  if (OR==T) {
  output_plot <-  output_plot +   
    scale_y_log10(limits = limits, breaks = breaks)
  }
  
  if (outer_limit_arrows) {
    stopifnot(length(limits) > 0)
    # Check for errorbar values outside 
    class_results[, "ymin_out"] <- ifelse(class_results[, ymin] < limits[1], limits[1], NA)
    class_results[, "linestyle_min"] <- ifelse(class_results[, ymin] < limits[1], "b", NA)
    class_results[, "ymax_out"] <- ifelse(class_results[, ymax] > limits[2], limits[2], NA)
    class_results[, "linestyle_max"] <- ifelse(class_results[, ymax] > limits[2], "b", NA)
    
    output_plot <- output_plot + geom_linerange(data = class_results, 
                                                aes_string(x = x, ymin = "ymin_out", ymax = ymax), linetype = 3, position = position_dodge(width = dodge_width))
    output_plot <- output_plot + geom_linerange(data = class_results, 
                                                aes_string(x = x, ymin = ymin, ymax = "ymax_out"), linetype = 3, position = position_dodge(width = dodge_width))
    output_plot <- output_plot + geom_linerange(data = class_results, 
                                                aes_string(x = x, ymin = "ymin_out", ymax = "ymax_out"), linetype = 3, position = position_dodge(width = dodge_width))
    
  }
  return(output_plot)
}

format_variable = function(vars) {
  as.character(vars) %>%
    str_remove_all(regex('binary|scaled|bin|verbose|_b$| b$|^g_|^ind_|^pct_|_c$', ignore_case = T)) %>%
    str_replace_all('_', " ") %>%
    str_replace_all('-$', " ") %>%
    str_replace_all('^M$', "Male") %>%
    str_replace_all('male', "Male") %>%
    str_replace_all('^F$', "Female") %>%
    str_replace_all('SystemTumorType', 'Tumor Type') %>%
    str_trim() %>%
    tolower() %>%
    tools::toTitleCase() %>%
    cap_susbet(cap_set)
}

cap_set = c('DNA', '^HR$', 'HRD', 'MMR', 'CHEK2', 'CH', 'PD', 'II',
 'DNMT3A', 'TET2', 'PPM1D', 'ASXL1', 'TP53', 'ATM', 'SF3B1', 'JAK2', 'SRSF2', '^VAF',
 '^ANC$', '^ALC$', '^AMC$', '^HGB$', '^MCV$', '^PLT$', '^WBC$', '^RDW$', 'XRT', 'IDH')

# capitalize a subset of words in a vector of string
cap_susbet = function(x, words) {
    for (word in words) {
      x = str_replace(x, regex(word, ignore_case = T), str_remove_all(word, regex("\\^|\\$| ")))
    }
    return(x)
}