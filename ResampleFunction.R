######################################################################
###   Distance & similarity based resampling of vegetation plots   ###
######################################################################

#Author: Jan Divíšek
#Version 2025-11-27

resample <- function(coord, spec, longlat = FALSE, dist.threshold = 1000, sim.threshold = 0.8,
                     sim.method = c("simpson", "sorensen", "jaccard", "bray"),
                     remove = c("random", "less diverse", "more diverse", "lower var.value", "higher var.value"),
                     var.value = NULL, strata = NULL, seed = 1234) {
  
  start_time <- Sys.time()
  
  # --- 1. Load packages ---
  pkgs <- c("spdep", "data.table", "Matrix", "vegan", "igraph")
  lapply(pkgs, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
  })
  
  spec[, PlotObservationID := as.character(PlotObservationID)]
  coord[, PlotObservationID := as.character(PlotObservationID)]
  
  # --- 2. Check data ---
  sim.method <- match.arg(sim.method); remove <- match.arg(remove)
  if(!("PlotObservationID" %chin% colnames(coord))) stop("Error: Column 'PlotObservationID' not found in 'coord'")
  if(ncol(coord) < 3) stop("Error: 'coord' must contain at least three columns (PlotObservationID, X, Y)")
  if(!all(sapply(coord[, 2:3], is.numeric))) stop("Error: Coordinates (X, Y or Lon, Lat) must be numeric")
  if(any(is.na(coord[, 2:3]))) stop("Error: Coordinates contain NA values")
  if(remove %chin% c("lower var.value", "higher var.value") && (is.null(var.value) || !var.value %chin% colnames(coord))) stop("Error: Invalid or missing var.value")
  if(!is.null(strata) && !strata %chin% colnames(coord)) stop("Error: Invalid strata column name")
  if(!all(c("PlotObservationID", "Taxon_name", "cover") %chin% colnames(spec))) stop("Error: Check column names in 'spec'")
  if(any(spec$cover == 0) && sim.method == "bray") stop("Error: Zero covers not allowed")
  if(any(is.na(spec$cover))) stop("Error: NAs in cover")
  if(!all(spec$PlotObservationID %chin% coord$PlotObservationID)) stop("Error: Some PlotObservationID in 'spec' not found in 'coord'")
  if(!all(coord$PlotObservationID %chin% spec$PlotObservationID)) stop("Error: Some PlotObservationID in 'coord' not found in 'spec'")
  
  # --- 3. Prepare data ---
  set.seed(seed) 
  coord <- coord[sample(1:nrow(coord)), ]
  
  if(remove == "less diverse"){
    coord <- coord[spec[, .(.N), by = .(PlotObservationID)], on = c("PlotObservationID")][order(N)]}
  if(remove == "more diverse"){
    coord <- coord[spec[, .(.N), by = .(PlotObservationID)], on = c("PlotObservationID")][order(-N)]}
  if(remove == "lower var.value"){
    coord <- setorderv(coord, var.value, 1, na.last=FALSE)}
  if(remove == "higher var.value"){
    coord <- setorderv(coord, var.value, -1, na.last=FALSE)}
  
  if (sim.method != "bray") { spec[, cover := 1] }
  
  ##set new ids based on ordered coord
  spec[coord[, .(.I, PlotObservationID)], on = "PlotObservationID", id := I]
  
  # --- 4. Indentification of neighbouring plots ---
  if (!is.null(strata)) {
    cat("Searching for neighbouring plots within strata. Please wait...\n")
    coord[, (strata) := as.character(get(strata))]
    d <- vector("list", nrow(coord))
    original_indices <- coord[,.I]
    
    for (s in unique(coord[[strata]])) {
      stratum_mask <- coord[[strata]] == s
      coord_sub <- coord[stratum_mask,]
      d_sub <- dnearneigh(as.matrix(coord_sub[, 2:3]), d1 = 0, d2 = dist.threshold, bounds = c("GE", "LT"), longlat = longlat)
      
      sub_indices <- original_indices[stratum_mask]
      for (i in 1:length(d_sub)) {
        d[[sub_indices[i]]] <- as.integer(sub_indices[d_sub[[i]]])
      }
    }
    class(d) <- "nb"
    rm(original_indices, stratum_mask, coord_sub, d_sub, sub_indices)
  } else { 
    cat("Searching for neighbouring plots. Please wait...\n")
    d <- dnearneigh(as.matrix(coord[, 2:3]), d1 = 0, d2 = dist.threshold, bounds = c("GE", "LT"), longlat = longlat) 
  }
  
  # --- 5. Split plots to groups ---
  g <- igraph::components(igraph::graph_from_adj_list(d, mode = "all"))
  
  spec[data.table(id = coord[,.I], group = g$membership), on = "id", grp := group]
  spec <- spec[grp %in% which(g$csize > 1), ]
  
  d <- mapply(append, seq_along(d), d, SIMPLIFY = FALSE, USE.NAMES = FALSE)
  names(d) <- coord[,.I]
  d <- d[g$membership %in% which(g$csize > 1)]
  
  filtering_task <- function(spec_sub, d, sim.threshold, sim.method){
    
    d_sub <- d[names(d) %chin% as.character(unique(spec_sub$id))]
    
    spec_sub <- Matrix::sparseMatrix(i = as.integer(factor(spec_sub$id), labels = names(d_sub)),
                                     j = as.integer(factor(spec_sub$Taxon_name)),
                                     x = spec_sub$cover,
                                     dimnames = list(names(d_sub), levels(factor(spec_sub$Taxon_name))))
    
    ds <- lapply(d_sub, FUN = calc_sim, spec_sub, sim.threshold, sim.method)
    
    pairs <- rbindlist(lapply(names(ds), function(p1) {
      p2 <- ds[[p1]][ds[[p1]] > sim.threshold & as.integer(names(ds[[p1]])) > as.integer(p1)]
      if (length(p2) == 0) return(NULL)
      data.table(p1 = p1, p2 = names(p2), sim = p2) }))
    
    if(nrow(pairs) > 0){
      pairs <- setorderv(pairs, c("sim"), c(-1))
      
      black <- NULL
      repeat
      {
        p <- pairs$p1[1]
        black <- c(black, as.integer(p))
        pairs <- pairs[!(pairs$p1 == p | pairs$p2 == p)]
        
        if(nrow(pairs) == 0){ break}
      }
      return(black)
    } else {
      return(NULL)
    }
    
  }
  
  calc_sim <- function(plots, spec_sub, sim.threshold, sim.method){
    
    n <- length(plots)-1
    
    sim <- switch(sim.method,
                  "simpson" = 1 - vegan::betadiver(as.matrix(spec_sub[as.character(plots), ]), method = "sim")[1:n],
                  "sorensen" = ,
                  "bray" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "bray")[1:n],
                  "jaccard" = 1 - vegan::vegdist(spec_sub[as.character(plots), ], method = "jaccard")[1:n])
    names(sim) <- plots[-1]
    return(sim)
    
  }
  
  cat("Similarity-based resampling:\n")
  pb <- txtProgressBar(min = 0, max = uniqueN(spec$grp), style = 3)
  blacklist <-  spec[, {setTxtProgressBar(pb, .GRP);
    filtering_task(.SD, d, sim.threshold, sim.method)}, 
    by = .(grp),
    .SDcols = PlotObservationID:grp]
  
  close(pb)
  
  # --- 8. Return selected plots ---
  
  coord.filtered <- coord[-blacklist$V1, 1:3]
  
  cat(sprintf("Removed %.1f%% of plots (%d out of %d)\n", nrow(blacklist) / nrow(coord) * 100, nrow(blacklist), nrow(coord)))
  end_time <- Sys.time(); elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat(sprintf("Elapsed time %02d:%02d:%02d\n", floor(elapsed / 3600), floor((elapsed %% 3600) / 60), round(elapsed %% 60)))
  
  return(coord.filtered)
}