ksp_to_clust_frame <- function(ksp) {
  ptid <- vapply(ksp$curves, "[[", numeric(1), "id")
  data.frame(ptid=ptid, cluster=ksp$clusters)
}

#' @export
merge_clusters <- function(ksp1, ksp2, ksp3, ksp_names) {
  d1 <- ksp_to_clust_frame(ksp1)
  names(d1)[2] <- ksp_names[1]
  d2 <- ksp_to_clust_frame(ksp2)
  names(d2)[2] <- ksp_names[2]
  d3 <- ksp_to_clust_frame(ksp3)
  names(d3)[2] <- ksp_names[3]
  dat <- merge(d1, d2, by="ptid", all=TRUE)
  dat <- merge(dat, d3, by="ptid", all=TRUE)
  dat
}

#' @import plyr
#' @export
make_ksplines_counts <- function(skin, pfvc, rvsp) {
  ksp_names <- c("skin", "pfvc", "rvsp")
  dat <- merge_clusters(skin, pfvc, rvsp, ksp_names)
  dat <- ddply(dat, ksp_names, summarize, n=length(ptid))
  arrange(dat, desc(n))
}

#' @import gridExtra
#' @export
plot_triple <- function(triple, skin_ksp, pfvc_ksp, rvsp_ksp) {
  all_clusters <- merge_clusters(skin_ksp, pfvc_ksp, rvsp_ksp, c("skin", "pfvc", "rvsp"))
  patients <- subset(all_clusters, skin == triple[1] & pfvc == triple[2] & rvsp == triple[3])
  patients <- patients$ptid
  subcurves <- list()
  subcurves$skin <- Filter(function(crv) crv$id %in% patients, skin_ksp$curves)
  subcurves$pfvc <- Filter(function(crv) crv$id %in% patients, pfvc_ksp$curves)
  subcurves$rvsp <- Filter(function(crv) crv$id %in% patients, rvsp_ksp$curves)

  blank <- grid.rect(gp=gpar(col="white"))
  p <- ggplot() + theme_bw()
  
  p1 <- p + geom_curves(subcurves$skin) +
      ylim(0, 75) + labs(x="years since diagnosis", y="total skin score")
  p2 <- p + geom_curves(subcurves$pfvc) +
      ylim(0, 120) + labs(x="years since diagnosis", y="pfvc")
  p3 <- p + geom_curves(subcurves$rvsp) +
      ylim(0, 120) + labs(x="years since diagnosis", y="rvsp")
  
  grid.arrange(p1, p2, p3, nrow=1, main=paste(triple, collapse=", "))
}
