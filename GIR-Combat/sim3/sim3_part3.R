

nbatches <- length(unlist(merge.tree))
nmerges <- nbatches - 1L
mnn.pairings <- left.set <- right.set <- vector("list", 
                                                nmerges)
for (mdx in seq_len(nmerges)) {
  # mdx = 1
  # mdx = 2
  next.step <- NEXT(merge.tree)
  left <- next.step$left
  right <- next.step$right
  
  .get_node_data <- function(node) node@data
  
  left.data <- .get_node_data(left)
  
  .get_node_restrict <- function(node) node@restrict
  
  left.restrict <- .get_node_restrict(left)
  
  .get_node_index <- function(node) node@index
  
  left.index <- .get_node_index(left)
  
  .get_node_origin <- function(node) node@origin
  
  left.origin <- .get_node_origin(left)
  
  .get_node_extras <- function(node) node@extras
  
  left.extras <- .get_node_extras(left)[[1]]
  right.data <- .get_node_data(right)
  right.restrict <- .get_node_restrict(right)
  right.index <- .get_node_index(right)
  right.origin <- .get_node_origin(right)
  right.extras <- .get_node_extras(right)[[1]]
  left.data <- as.matrix(left.data)
  right.data <- as.matrix(right.data)
  
  .restricted_mnn  <- function (left.data, left.restrict, right.data, right.restrict, 
                                k, prop.k = NULL, ...) 
  {
    if (!is.null(left.restrict)) {
      left.data <- left.data[left.restrict, , drop = FALSE]
    }else {
      left.data <- left.data
    }
    if (!is.null(right.restrict)) {
      right.data <- right.data[right.restrict, , drop = FALSE]
    }else {
      right.data <- right.data
    }
    
    .choose_k <-function (k, prop.k, N) 
    {
      if (is.null(prop.k)) {
        k
      }else {
        min(N, max(k, round(prop.k * N)))
      }
    }
    
    k1 <- .choose_k(k, prop.k, nrow(left.data))
    k2 <- .choose_k(k, prop.k, nrow(right.data))
    
    findMutualNN <- function (data1, data2, k1, k2 = k1, BNINDEX1 = NULL, BNINDEX2 = NULL) 
    {
      data1 <- as.matrix(data1)
      data2 <- as.matrix(data2)
      library(BiocNeighbors)
    
      W21 <- queryAnnoy(data2, query=data1, k=k2, get.distance=FALSE) # W21 <- do.call(queryKNN, args)
    
      W12 <- queryAnnoy(data1, query=data2, k=k1, get.distance=FALSE) # W12 <- do.call(queryKNN, args)
      
      find_mutual_nns <- function(left, right) {
        .Call('_batchelor_find_mutual_nns', PACKAGE = 'batchelor', left, right)
      }
      
      out <- find_mutual_nns(W21$index, W12$index)
      names(out) <- c("first", "second")
      return(out)
    }
    
    # pairs <- findMutualNN(left.data, right.data, k1 = k1, k2 = k2, ...)
    pairs <- findMutualNN(left.data, right.data, k1 = k1, k2 = k2)
    
    .unrestrict_indices<- function (index, restrict) 
    {
      if (!is.null(restrict)) 
        index <- restrict[index]
      index
    }
    
    pairs$first <- .unrestrict_indices(pairs$first, left.restrict)
    pairs$second <- .unrestrict_indices(pairs$second, right.restrict)
    pairs
  }
  
  mnn.sets <- .restricted_mnn(left.data, left.restrict, 
                              right.data, right.restrict, k = k, prop.k = prop.k)
  s1 <- mnn.sets$first
  s2 <- mnn.sets$second
  mnn.pairings[[mdx]] <- DataFrame(left = s1, right = s2)
  left.set[[mdx]] <- left.index
  right.set[[mdx]] <- right.index
  trans.right <- t(right.data)
  
  .compute_correction_vectors <- function (data1, data2, mnn1, mnn2, tdata2, sigma) 
  {
 
    
    vect <- data1[mnn1, , drop = FALSE] - data2[mnn2, , drop = FALSE]
    averaged <- sumCountsAcrossCells(t(vect), DataFrame(ID = mnn2), 
                                     average = TRUE)
    
    
    library(Rcpp)

    Rcpp::sourceCpp('D:/Users/lenovo/Desktop/小论文/dusmooth_gaussian_kernel.cpp')
    
    cell.vect <- dusmooth_gaussian_kernel(assay(averaged, withDimnames = FALSE),
                                          averaged$ID - 1L, tdata2, sigma)
    
    # cell.vect <- smooth_gaussian_kernel(assay(averaged, withDimnames = FALSE), 
    #                                     averaged$ID - 1L, tdata2, sigma)
    t(cell.vect)
  }
  
  correction.in <- .compute_correction_vectors(left.data, 
                                               right.data, s1, s2, trans.right, sigma)
  if (!same.set) {
    correction.out <- .compute_correction_vectors(left.extras, 
                                                  right.extras, s1, s2, trans.right, sigma)
  }
  if (svd.dim > 0) {
    u1 <- unique(s1)
    u2 <- unique(s2)
    in.span1 <- .get_bio_span(t(left.data[u1, , drop = FALSE]), 
                              ndim = svd.dim, BSPARAM = BSPARAM, BPPARAM = BPPARAM)
    in.span2 <- .get_bio_span(t(right.data[u2, , drop = FALSE]), 
                              ndim = svd.dim, BSPARAM = BSPARAM, BPPARAM = BPPARAM)
    correction.in <- .subtract_bio(correction.in, in.span1, 
                                   in.span2)
    if (!same.set) {
      out.span1 <- .get_bio_span(t(left.extras[u1, 
                                               , drop = FALSE]), subset.row = subset.row, 
                                 ndim = svd.dim, BSPARAM = BSPARAM, BPPARAM = BPPARAM)
      out.span2 <- .get_bio_span(t(right.extras[u2, 
                                                , drop = FALSE]), subset.row = subset.row, 
                                 ndim = svd.dim, BSPARAM = BSPARAM, BPPARAM = BPPARAM)
      correction.out <- .subtract_bio(correction.out, 
                                      out.span1, out.span2, subset.row = subset.row)
    }
  }
  if (var.adj) {
    args <- list(sigma = sigma, restrict1 = left.restrict, 
                 restrict2 = right.restrict)
    
    .adjust_shift_variance <- function(data1, data2, correction, sigma, subset.row=NULL, restrict1=NULL, restrict2=NULL) 
      # Performs variance adjustment to avoid kissing effects.    
    {
      cell.vect <- correction 
      
      if (!is.null(subset.row)) { 
        # Only using subsetted genes to compute locations, consistent with our policy in SVD. 
        cell.vect <- cell.vect[,subset.row,drop=FALSE]
        data1 <- data1[subset.row,,drop=FALSE]
        data2 <- data2[subset.row,,drop=FALSE]
      }
      
      .col_subset_to_index <- function(x, index) {
        if (is.null(index)) {
          seq_len(ncol(x))
        } else {
          i <- seq_len(ncol(x))
          names(i) <- colnames(x)
          unname(i[index])
        }
      }
      
      restrict1 <- .col_subset_to_index(data1, restrict1) - 1L
      restrict2 <- .col_subset_to_index(data2, restrict2) - 1L
      
      adjust_shift_variance <- function(data1, data2, vect, sigma2, restrict1, restrict2) {
        .Call('_batchelor_adjust_shift_variance', PACKAGE = 'batchelor', data1, data2, vect, sigma2, restrict1, restrict2)
      }
      
      scaling <- adjust_shift_variance(data1, data2, cell.vect, sigma, restrict1, restrict2)
      
      scaling <- pmax(scaling, 1)
      scaling * correction
    }
    
    correction.in <- do.call(.adjust_shift_variance, 
                             c(list(t(left.data), t(right.data), correction.in), 
                               args))
    if (!same.set) {
      correction.out <- do.call(.adjust_shift_variance, 
                                c(list(t(left.extras), t(right.extras), correction.out, 
                                       subset.row = subset.row), args))
    }
  }
  right.data <- right.data + correction.in
  if (!same.set) {
    right.extras <- right.extras + correction.out
  }
  
  MNN_treenode <- function(index, data, restrict, origin=rep(index, nrow(data)), extras=list()) {
    new("MNN_treenode", index=index, data=data, restrict=restrict, origin=origin, extras=extras)
  }
  
  .combine_restrict <- function(left.data, left.restrict, right.data, right.restrict) {
    if (is.null(left.restrict) && is.null(right.restrict)) {
      NULL
    } else {
      if (is.null(left.restrict)) {
        left.restrict <- seq_len(nrow(left.data))
      }
      if (is.null(right.restrict)) {
        right.restrict <- seq_len(nrow(right.data))
      }
      c(left.restrict, right.restrict + nrow(left.data))
    }
  }
  
  merge.tree <- UPDATE(merge.tree, next.step$chosen, data = rbind(left.data, 
                                                                  right.data), index = c(left.index, right.index), 
                       restrict = .combine_restrict(left.data, left.restrict, 
                                                    right.data, right.restrict), origin = c(left.origin, 
                                                                                            right.origin), extras = list(rbind(left.extras, 
                                                                                                                               right.extras)))
}
full.order <- .get_node_index(merge.tree)
full.origin <- .get_node_origin(merge.tree)
if (same.set) {
  full.data <- .get_node_data(merge.tree)
}else {
  full.data <- .get_node_extras(merge.tree)[[1]]
}
for (mdx in seq_along(mnn.pairings)) {
  bonus1 <- match(left.set[[mdx]][1], full.origin) - 1L
  mnn.pairings[[mdx]]$left <- mnn.pairings[[mdx]]$left + 
    bonus1
  bonus2 <- match(right.set[[mdx]][1], full.origin) - 
    1L
  mnn.pairings[[mdx]]$right <- mnn.pairings[[mdx]]$right + 
    bonus2
}
if (is.unsorted(full.order)) {
  ncells.per.batch <- tabulate(full.origin)
  ordering <- .restore_original_order(full.order, ncells.per.batch)
  full.data <- full.data[ordering, , drop = FALSE]
  full.origin <- full.origin[ordering]
  mnn.pairings <- .reindex_pairings(mnn.pairings, ordering)
}


output <-SingleCellExperiment(list(corrected = t(full.data)), colData = DataFrame(batch = full.origin), 
                              metadata = list(merge.info = DataFrame(left = I(as(left.set, 
                                                                                 "List")), right = I(as(right.set, "List")), pairs = I(as(mnn.pairings, 
                                                                                                                                          "List")))))





nms <- names(batches)
if (!is.null(nms)) {
  output$batch <- nms[output$batch]
  
  .fix_names_in_merge_info <- function (output, names) 
  {
    L1 <- metadata(output)$merge.info$left
    R1 <- metadata(output)$merge.info$right
    L2 <- R2 <- List()
    for (i in seq_along(L1)) {
      L2[[i]] <- names[L1[[i]]]
      R2[[i]] <- names[R1[[i]]]
    }
    metadata(output)$merge.info$left <- L2
    metadata(output)$merge.info$right <- R2
    output
  }
  
  output <- .fix_names_in_merge_info(output, nms)
}

output

if (do.split) {
  d.reo <- divided$reorder
  output <- output[, d.reo, drop = FALSE]
  
  .reindex_pairings <- function(pairings, new.order) 
    # Restores the MNN pairing indices to the original positions,
    # i.e., after cells have been reordered by 'new.order'.
  {
    rev.order <- integer(length(new.order))
    rev.order[new.order] <- seq_along(new.order)
    for (x in seq_along(pairings)) {
      current <- pairings[[x]]
      current$left <- rev.order[current$left]
      current$right <- rev.order[current$right]
      pairings[[x]] <- current
    }
    pairings
  }
  
  metadata(output)$merge.info$pairs <- .reindex_pairings(metadata(output)$merge.info$pairs, 
                                                         d.reo)
}

.rename_output<- function (output, batches, subset.row = NULL, correct.all = FALSE, 
                           cells.in.columns = TRUE) 
{
  GENERATE_NAMES <- function(batches, OTHERDIMFUN, OTHERDIMNAMEFUN) {
    collected <- lapply(batches, OTHERDIMNAMEFUN)
    nulled <- vapply(collected, is.null, FUN.VALUE = TRUE)
    if (any(nulled) && !all(nulled)) {
      collected[nulled] <- lapply(batches[nulled], FUN = function(x) character(OTHERDIMFUN(x)))
    }
    unlist(collected)
  }
  if (!cells.in.columns) {
    cell.names <- GENERATE_NAMES(batches, nrow, rownames)
    rownames(output) <- cell.names
    colnames(output) <- colnames(batches[[1]])
  }
  else {
    cell.names <- GENERATE_NAMES(batches, ncol, colnames)
    colnames(output) <- cell.names
    feat.names <- rownames(batches[[1]])
    if (!is.null(feat.names) && !is.null(subset.row) && 
        !correct.all) {
      feat.names <- feat.names[.row_subset_to_index(batches[[1]], 
                                                    subset.row)]
    }
    rownames(output) <- feat.names
  }
  output
}

.rename_output(output, original, subset.row = subset.row, 
               correct.all = correct.all)


Xmnn <- output
