#' bead array control processing/QC: get metrics from control probe signals
#' 
#' mostly lifted from Sean Maden's recountmethylation supplement
#' https://github.com/metamaden/recountmethylationManuscriptSupplement
#' modified slightly to take advantage of rgSet annotations and to Go Big,
#' i.e. this works reasonably well on fairly large rgSets (like all of GEO)
#' 
#' @param rgSet           an rgSet with control probe intensities and annotation
#' @param dft             optional single-row data.frame of thresholds (default)
#' @param baseline        offset for Extension green based background (3000)
#' @param biotin.baseline baseline offset for biotin staining probes (1)
#'
#' @return                a matrix (cols = metrics, rows = samples)
#' 
#' @examples 
#'   
#' if (exists("MPAL_rgSet")) {
#'
#'   # TARGET MPAL data, from Alexander et al, Nature 2018
#'   cm <- control_metrics(MPAL_rgSet)
#'
#'   library(ComplexHeatmap)
#'   Heatmap(t(cm), name="metric", column_names_side="top",
#'           column_names_gp=gpar(fontsize=6), row_names_side="left")
#'
#'   library(circlize)
#'   flagcols <- colorRamp2(c(0, 1), c("white", "darkred"))
#'   flagged <- flag_control_failures(cm, platform="epic")
#'   all_samples_passed <- ifelse(colSums(flagged) == 0, 1, 0)
#'   all_probes_passed <- ifelse(rowSums(flagged) == 0, 1, 0)
#'   Heatmap(t(flagged), name="failed", col=flagcols, column_names_side="top",
#'           column_split=all_probes_passed, column_names_gp=gpar(fontsize=6), 
#'           row_split=all_samples_passed, row_names_side="left")
#'
#' } 
#' 
#' @seealso               flag_control_failures
#' 
#' @import                minfi
#'
#' @export 
control_metrics <- function(rgSet, dft=NULL, baseline=3000, biotin.baseline=1) {

  # make this easier on RAM by NOT grabbing hundreds of unneeded control probes 
  #
  cdf <- subset(getProbeInfo(rgSet, "Control"), !grepl("(NORM|NEGATIVE)", Type))

  message("Retrieving control probe red intensities...") 
  rs <- as.matrix(t(getRed(rgSet[cdf$Address,])))
  
  message("Retrieving control probe green intensities...") 
  gs <- as.matrix(t(getGreen(rgSet[cdf$Address,])))

  message("Loading thresholds...")
  if (is.null(dft)) dft <- .get_dft() 
  cnl <- names(dft) 

  # just dimnames
  rmat <- rs[, 0]


  #---------------------
  # Background addresses
  #---------------------
  # 1 metric, grn channel
  #
  ci <- .getsub(cdf, "EXTENSION")
  addr.bkg <- subset(ci, grepl("\\([AT]\\)", ExtendedType))$Address 
 

  #---------------------
  # FFPE Restoration kit
  #---------------------
  # 1 metric, grn channel
  #
  message("calculating restore metric...")
  ci = .getsub(cdf, "RESTORATION")
  which.rest = which(colnames(gs) %in% ci$Address)
  m1 = apply(gs, 1, function(x) x[which.rest]/(max(x[addr.bkg]) + baseline))
  mi = 1
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #----------------
  # BIOTIN STAINING
  #----------------
  # 2 metrics, 1 per chan
  #
  message("calculating biotin staining metrics...")
  ci = cdf[grepl("Biotin|DNP", cdf$ExtendedType),]
  
  # red 
  which.stain = which(colnames(rs) %in% .getetaddr(ci, "DNP \\(High"))
  which.bkg = which(colnames(rs) %in% .getetaddr(ci, "DNP \\(Bkg"))
  m1 = apply(rs, 1, function(x) x[which.stain]/(x[which.bkg] + biotin.baseline))
  mi = 2
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]

  # grn
  which.stain = which(colnames(gs) %in% .getetaddr(ci, "Biotin \\(High"))
  which.bkg = which(colnames(gs) %in% .getetaddr(ci, "Biotin \\(Bkg"))
  m2 = apply(gs, 1, function(x) x[which.stain]/(x[which.bkg] + biotin.baseline))
  mi = 3
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #--------------
  # SPECIFICITY I
  #--------------
  # 2 metrics, 1 per channel
  #
  message("calculating specificity I...")
  ci1 <- .getsub(cdf, "SPECIFICITY I")

  # red
  ci = ci1[grepl("Mismatch (4|5|6)", ci1$ExtendedType),]
  addr.mm.index = which(colnames(rs) %in% .getetaddr(ci, "MM"))
  addr.pm.index = which(colnames(rs) %in% .getetaddr(ci, "PM"))
  m1 = apply(rs, 1, function(x) min(x[addr.pm.index])/max(x[addr.mm.index]))
  mi = 4
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]

  # grn
  ci = ci1[grepl("Mismatch (1|2|3)", ci1$ExtendedType),]
  addr.mm.index = which(colnames(gs) %in% .getetaddr(ci, "MM"))
  addr.pm.index = which(colnames(gs) %in% .getetaddr(ci, "PM"))
  m2 = apply(gs, 1, function(x) min(x[addr.pm.index])/max(x[addr.mm.index]))
  mi = 5
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #----------------
  # SPECIFICITY II
  #----------------
  # 1 metric, uses both channels
  # 
  message("calculating specificityII metric...")
  ci = .getsub(cdf, "SPECIFICITY II")
  which.addr.red = which(colnames(rs) %in% ci$Address)
  which.addr.grn = which(colnames(gs) %in% ci$Address)
  m0.1 = apply(rs, 1, function(x) min(x[which.addr.red]))
  m0.2 = apply(gs, 1, function(x) max(x[which.addr.grn]))
  m1 = m0.1/m0.2
  mi = 6
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #-----------
  # EXTENSION
  #-----------
  # 2 metrics, 1 per channel
  #
  message("calculating extension metrics...")
  ci = .getsub(cdf, "EXTENSION")
  addr.ext.cg = ci$Address[grepl("(C)|(G)", ci$ExtendedType)]
  addr.ext.at = ci$Address[grepl("(A)|(T)", ci$ExtendedType)]

  # red
  which.cg = which(colnames(rs) %in% addr.ext.cg)
  which.at = which(colnames(rs) %in% addr.ext.at)
  m1 = apply(rs, 1, function(x) min(x[which.at])/max(x[which.cg]))
  mi = 7
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]

  # grn
  which.cg = which(colnames(gs) %in% addr.ext.cg)
  which.at = which(colnames(gs) %in% addr.ext.at)
  m2 = apply(gs, 1, function(x) min(x[which.cg])/max(x[which.at]))
  mi = 8
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #---------------
  # HYBRIDIZATION
  #---------------
  # 2 metrics, 1 per channel
  #
  message("calculating hybridization metrics...")
  ci = .getsub(cdf, "HYBRIDIZATION")
  which.hi = which(colnames(gs) %in% .getetaddr(ci, "High"))
  which.med = which(colnames(gs) %in% .getetaddr(ci, "Med"))
  which.low = which(colnames(gs) %in% .getetaddr(ci, "Low"))
  m1 = apply(gs, 1, function(x) x[which.hi]/x[which.med])
  m2 = apply(gs, 1, function(x) x[which.med]/x[which.low])
  mi = 9
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
  mi = 10
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #----------------
  # TARGET REMOVAL
  #----------------
  # 2 metrics, both from grn channel
  #
  message("calculating target removal metics...")
  ci = .getsub(cdf, "EXTENSION")
  which.bkg = which(colnames(gs) %in% .getetaddr(ci, "(A)|(T)"))

  ci = .getsub(cdf, "TARGET REMOVAL")
  which.t1 = which(colnames(gs) %in% .getetaddr(ci, "Removal 1"))
  which.t2 = which(colnames(gs) %in% .getetaddr(ci, "Removal 2"))

  # rem 1
  m1 = apply(gs, 1, function(x) (max(x[which.bkg]) + baseline)/x[which.t1])
  mi = 11
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]

  # rem 2
  m2 = apply(gs, 1, function(x) (max(x[which.bkg]) + baseline)/x[which.t2]) 
  mi = 12
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]


  #------------------------
  # BISULFITE CONVERSION I
  #------------------------
  # 2 metrics, 1 per channel
  #
  # Somewhat useless (see Zhou et al, NAR 2017)
  #
  message("calculating bisulfite conversion I metrics...")
  ci = .getsub(cdf, "BISULFITE CONVERSION I")
  which.c123 = which(colnames(gs) %in% .getetaddr(ci, "C[123]"))
  which.u123 = which(colnames(gs) %in% .getetaddr(ci, "U[123]"))
  which.c456 = which(colnames(gs) %in% .getetaddr(ci, "C[456]"))
  which.u456 = which(colnames(gs) %in% .getetaddr(ci, "U[456]"))

  # red
  m1 = apply(rs, 1, function(x) min(x[which.c456])/max(x[which.u456]))
  mi = 13
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]

  # grn
  m2 = apply(gs, 1, function(x) min(x[which.c123])/max(x[which.u123]))
  mi = 14
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 

  #-------------------------
  # BISULFITE CONVERSION II
  #-------------------------
  # 2 metrics, 1 per channel
  #
  # Somewhat useless (see Zhou et al, NAR 2017)
  #
  message("calculating bisulfite conversion II metric...")
  ci = .getsub(cdf, "BISULFITE CONVERSION II")
  which.ci.red = which(colnames(rs) %in% ci$Address)
  which.ci.grn = which(colnames(gs) %in% ci$Address)
  m0.1 = apply(rs, 1, function(x) min(x[which.ci.red]))
  m0.2 = apply(gs, 1, function(x) max(x[which.ci.grn]))
  m1 = m0.1/m0.2
  mi = 15
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]


  #----------------
  # NON-POLYMORPHIC
  #----------------
  # 2 metrics, 1 per channel
  #
  message("calculating non-polymorphic metrics...")
  ci = .getsub(cdf, "NON-POLYMORPHIC")

  # red
  which.cg = which(colnames(rs) %in% .getetaddr(ci, "(C|G)"))
  which.at = which(colnames(rs) %in% .getetaddr(ci, "(A|T)"))
  m1 = apply(rs, 1, function(x) min(x[which.at])/max(x[which.cg]))
  mi = 16
  rmat = cbind(rmat, m1)
  colnames(rmat)[ncol(rmat)] = cnl[mi]

  # grn
  which.cg = which(colnames(rs) %in% .getetaddr(ci, "(C|G)"))
  which.at = which(colnames(rs) %in% .getetaddr(ci, "(A|T)"))
  m2 = apply(gs, 1, function(x) min(x[which.cg])/max(x[which.at]))
  mi = 17
  rmat = cbind(rmat, m2)
  colnames(rmat)[ncol(rmat)] = cnl[mi]
 
  #-----
  # DONE
  #----- 
  return(rmat)

}


# helper 
.getsub <- function(cdf, ct) subset(cdf, Type == ct)


# helper 
.getetaddr <- function(cdf, patt) subset(cdf, grepl(patt, ExtendedType))$Address


# helper
.get_dft <- function() { 

  # should go into data(dft) or similar
  dft <- data.frame(restoration.grn = 0, 
                    biotin.stain.red = 5, 
                    biotin.stain.grn = 5,
                    specificityI.red = 1, 
                    specificityI.grn = 1, 
                    specificityII = 1,
                    extension.red = 5, 
                    extension.grn = 5, 
                    hyb.hi.med = 1, 
                    hyb.med.low = 1,
                    target.removal.1 = 1, 
                    target.removal.2 = 1, 
                    bisulfite.conv.I.red = 1,
                    bisulfite.conv.I.grn = 1, 
                    bisulfite.conv.II = 1, 
                    nonpolymorphic.red = 5, 
                    nonpolymorphic.grn = 5) 

  return(dft)

}
