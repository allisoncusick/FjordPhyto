tag_facet_manual <-  function(p, open=c("(",""), close = c(")","."),
                        tag_fun_top = function(i) letters[i],
                        x = 0, y = 0.5, t_adjust = NULL,
                        hjust = 0, vjust = 0.5, 
                        fontface = c(2,2), ...){
  
  if(is.null(t_adjust)){
    t_adjust = 8
  }
  
  p <- heatmap_time
  ## specific y axis values: gb$layout$panel_scales_y[[1]]$get_labels()
  ylabel <- function(label1,label2){
    L1 <- nchar(label1)
    L2 <- nchar(label2)
    scaler <- ifelse(L1 + L2 > 8, 4, 0)
    space1 = paste0(rep("",27 - (L1/2)),collapse = " ")
    space2 = paste0(rep("",44 - (L1/2 + L2/2) - scaler), collapse = " ")
    space3 = paste0(rep("",22 - (L2/2)), collapse = " ")
    paste0(space1,label1,space2,label2,space3)
  }
  
  p + ylab(ylabel("Testing","Testing")) + 
                        
  gb <- ggplot_build(p)
  lay <- gb$layout$layout
  nm <- names(gb$layout$facet$params$rows)
  nm_lvl <- levels(gb$layout$layout$phytogroups)
  
  
  
  tl <- lapply(tags_top, grid::textGrob, x=x[1], y=y[1],
               hjust=hjust[1], vjust=vjust[1], gp=grid::gpar(fontface=fontface[1]))

  facet_move <- lapply(nm_lvl, grid::textGrob, x=x[1], y=y[1],
                       hjust=hjust[1], vjust=vjust[1], gp=grid::gpar(fontface=fontface[1]))
  
  g <- ggplot_gtable(gb)
  
  which.ylab = grep('ylab-l', g$layout$name)
  which.axes = grep('axis-l', g$layout$name)
  axis.rows  = g$layout$t[which.axes]-1
  label.col  = g$layout$l[which.ylab]+1
  which.facet =  grep('strip-r', g$layout$name)
  
  g = gtable::gtable_add_grob(g, g$grobs[which.facet], axis.rows, label.col)
  g = gtable::gtable_add_row_space(g, grid::unit(0.2,"line"))
  #g = gtable::gtable_filter(g, 'ylab-l') 
  grid::grid.draw(g)
  
  g <- gtable::gtable_add_rows(g, grid::unit(0.1,"line"), pos = 1)
  l <- unique(g$layout[grepl("panel",g$layout$name), "l"])
  g <- gtable::gtable_add_grob(g, grobs = tl, t=t_adjust, l=l,clip = "off", b = 10)
  

  
  grid::grid.newpage()
  grid::grid.draw(g)
}