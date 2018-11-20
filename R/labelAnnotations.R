#' label ComplexHeatmap annotations
#' 
#' @param annoTop     top annotations
#' @param annoBot     bottom annotations
#' @param fontsize    font size (default 11) 
#' 
#' @import ComplexHeatmap 
#' 
#' @export
labelAnnotations <- function(annoTop=NULL, annoBot=NULL, fontsize=11) { 
  if (!is.null(annoTop) && 
      is(annoTop, "list") && 
      all(c("annoTop", "annoBot") %in% names(annoTop))) {
    annoBot <- annoTop$annoBot
    annoTop <- annoTop$annoTop  
  }
  message("Labeling bottom annotations...")
  if (!is.null(annoBot) && "anno_list" %in% slotNames(annoBot)) {
    for(an in sapply(annoBot@anno_list, slot, "name")) {
      decorate_annotation(an, {
        grid.text(an, unit(1,"npc") + unit(2,"mm"), 0.5,
                  default.units="npc", just="left", gp=gpar(fontsize=fontsize))
      })
    }
  }
  if (!is.null(annoTop) && "anno_list" %in% slotNames(annoTop)) {
    message("Labeling top annotations...")
    for(an in sapply(annoTop@anno_list, slot, "name")) {
      decorate_annotation(an, {
        grid.text(an, unit(1,"npc") + unit(2,"mm"), 0.5,
                  default.units="npc", just="left", gp=gpar(fontsize=fontsize))
      })
    }
  }
}

# helper fn
.labelBottomAnnotations <- function(annoBot, fontsize=11) { 
  if (!is.null(annoBot) && "anno_list" %in% slotNames(annoBot)) {
    message("Labeling bottom annotations...")
    for(an in sapply(annoBot@anno_list, slot, "name")) {
      decorate_annotation(an, {
        grid.text(an, unit(1,"npc") + unit(2,"mm"), 0.5,
                  default.units="npc", just="left", gp=gpar(fontsize=fontsize))
      })
    }
  }
}
