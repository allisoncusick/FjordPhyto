scale_discrete_manual_ext <- function(aesthetics, values, name, labels)
{
  lapply(aesthetics, function(aesthetic) {
    ggplot2::scale_discrete_manual(aesthetic, values = values[[aesthetic]], name = name, labels = labels, drop = F)
  })
}

ggnested_pattern <- function (data, mapping = aes(), ...,
                       legend_labeling = c("sub", "join", "main"),
                       join_str = " - ", legend_title = NULL,
                       main_keys = TRUE, nested_aes = c("fill", "color", "pattern"),
                       gradient_type = c("both", "shades", "tints"), min_l = 0.05, max_l = 0.95,
                       main_palette = NULL, 
                       base_clr = "#008CF0") 
{
  aes_args <- names(mapping)
  if (!"main_group" %in% aes_args) {
    stop("Error: provide the main_group in the aesthetic mapping argument. For non-nested data, use the regular ggplot2 function.")
  }
  if (!"sub_group" %in% aes_args) {
    stop("Error: provide a subgroup in the aesthetic mapping argument. For non-nested data, use the regular ggplot2 function.")
  }
  if ("fill" %in% aes_args & "fill" %in% nested_aes) {
    warning("Warning: fill aesthetics will be ignored in the main ggnested function. Please specify non-nested fill in the geom_* layer. Alternatively,\n            remove 'fill' from mapping_aes.")
    mapping$fill <- NULL
  }
  if (("colour" %in% aes_args | "color" %in% aes_args) & ("colour" %in% 
                                                          nested_aes | "color" %in% nested_aes)) {
    warning("Warning: colour aesthetics will be ignored in the main ggnested function. Please specify non-nested colour in the geom_* layer. Alternatively,\n            remove 'colour' from mapping_aes.")
    mapping$colour <- NULL
    mapping$color <- NULL
  }
  group <- quo_name(mapping$main_group)
  subgroup <- quo_name(mapping$sub_group)
  patgroup <- quo_name(mapping$pattern)
  
  pal <- nested_palette(data, group, subgroup, gradient_type, 
                        min_l, max_l, main_palette, base_clr, join_str) %>%
    left_join(., data[,c(group, subgroup, patgroup)])
  
  colours <- pal %>% rename(sublabel = !!subgroup, label = !!group, pattern = !!patgroup) %>% 
    as.data.frame() 
  
  if (main_keys) {
    colours <- colours %>% group_by(label) %>%
      group_modify(~add_row(.x,.before = 0)) %>% ungroup() %>%
      mutate(subgroup_colour = ifelse(is.na(subgroup_colour),"#FFFFFF", subgroup_colour),
             sublabel = ifelse(is.na(sublabel),
                               sprintf("**%s**", as.character(label)), as.character(sublabel)),
             group_subgroup = ifelse(is.na(group_subgroup),
                                     sprintf("**%s**",as.character(label)), group_subgroup)) %>%
      as.data.frame() %>%
      mutate(pattern = ifelse(is.na(pattern), "none", pattern))
  }
  vals <- colours$subgroup_colour
  names(vals) <- colours$group_subgroup
  df <- left_join(data, pal, by = c(group, subgroup)) %>% 
    arrange(group, subgroup) %>%
    mutate(group_subgroup = factor(
      group_subgroup, ordered = T, levels = colours$group_subgroup),
      `:=`(!!subgroup,  factor(!!sym(subgroup), ordered = T)),
      `:=`(!!group, factor(!!sym(group), ordered = T))) %>%
    ungroup() %>% 
    arrange(group_subgroup)
  if (legend_labeling[1] == "join") {
    labels <- colours$group_subgroup
    leg_title <- sprintf("%s%s%s", group, join_str, subgroup)
  }
  else if (legend_labeling[1] == "main") {
    labels <- colours$label
    leg_title <- group
  }
  else if (legend_labeling[1] == "sub") {
    labels <- colours$sublabel
    leg_title <- subgroup
  }
  else {
    stop("Invalid option for legend_labeling. Pick one of c('join', 'main', 'sub')")
  }
  if (!is.null(legend_title)) {
    leg_title <- legend_title
  }
  vals_pat <- colours$pattern
  names(vals_pat) <- colours$group_subgroup
  nested_scale <- scale_discrete_manual_ext(..., aesthetics = nested_aes, name = leg_title,
                                            values = list(
                                              fill = vals,
                                              color = vals,
                                              pattern = vals_pat
                                            ),
                                            labels = labels)
  
  if ("fill" %in% nested_aes) {
    mapping$fill <- quo(group_subgroup)
  }
  if ("colour" %in% nested_aes | "color" %in% nested_aes) {
    mapping$colour <- quo(group_subgroup)
  }
  if ("pattern" %in% nested_aes | "pattern" %in% nested_aes) {
    mapping$pattern <- quo(group_subgroup)
  }
  p <- ggplot(df, mapping, ...) + nested_scale
  if (main_keys) {
    p <- p + theme_nested(theme) + 
      guides(fill = guide_legend(override.aes = list(pattern = vals_pat)))

  }
  return(p)
}
