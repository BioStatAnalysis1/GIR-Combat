
batch = old_pheno$measure
restrict = NULL
k=20
prop.k = NULL
sigma=1
var.adj=TRUE; subset.row = NULL; correct.all = FALSE;
merge.order = NULL; auto.merge = FALSE;assay.type = "logcounts";

cos.norm.in = TRUE 
cos.norm.out = TRUE
svd.dim = 0L 

library("batchelor")


.unpackLists<- function (data) 
{
  objects <- list(data)
  for (i in seq_along(objects)) {
    current <- objects[[i]]
    if (!is.list(current)) {
      if (is(current, "List") && !is(current, "DataFrame")) {
        current <- as.list(current)
      }
      else {
        current <- list(current)
      }
      objects[[i]] <- current
    }
  }
  do.call(c, objects)
}

original <- batches <- .unpackLists(data)

checkBatchConsistency <- function (batches, cells.in.columns = TRUE) 
{
  if (length(batches) == 0L) {
    return(invisible(NULL))
  }
  if (cells.in.columns) {
    DIMFUN <- nrow
    DIMNAMEFUN <- rownames
    DIM <- "row"
  }
  else {
    DIMFUN <- ncol
    DIMNAMEFUN <- colnames
    DIM <- "column"
  }
  first <- batches[[1]]
  ref.n <- DIMFUN(first)
  ref.names <- DIMNAMEFUN(first)
  for (b in seq_along(batches)[-1]) {
    current <- batches[[b]]
    if (!identical(DIMFUN(current), ref.n)) {
      stop(sprintf("number of %ss is not the same across batches (see batch %s)", 
                   DIM, .identify_failed_batch(b, names(batches))))
    }
    cur.names <- DIMNAMEFUN(current)
    if (!identical(cur.names, ref.names)) {
      stop(sprintf("%s names are not the same across batches (see batch %s)", 
                   DIM, .identify_failed_batch(b, names(batches))))
    }
  }
  invisible(NULL)
}

checkBatchConsistency(batches)

checkRestrictions <- function (batches, restrictions, cells.in.columns = TRUE) 
{
  if (is.null(restrictions)) {
    return(NULL)
  }
  if (length(batches) != length(restrictions)) {
    stop("'restrictions' must of length equal to the number of batches")
  }
  if (!identical(names(batches), names(restrictions))) {
    stop("'restrictions' must have the same names as the batches")
  }
  for (b in seq_along(batches)) {
    if (is.null(restrictions[[b]])) {
      next
    }
    FUN <- if (!cells.in.columns) 
      .row_subset_to_index
    else .col_subset_to_index
    restrictions[[b]] <- FUN(batches[[b]], restrictions[[b]])
    if (length(restrictions[[b]]) == 0L) {
      stop("no cells remaining in a batch after restriction")
    }
  }
  restrictions
}

restrict <- checkRestrictions(batches, restrict)

checkIfSCE <-  function (batches) 
{
  vapply(batches, is, class2 = "SingleCellExperiment", FUN.VALUE = TRUE)
}

is.sce <- checkIfSCE(batches)
if (any(is.sce)) {
  batches[is.sce] <- lapply(batches[is.sce], assay, i = assay.type,
                            withDimnames = FALSE)
}

do.split <- length(batches) == 1L
if (do.split) {
  divided <- divideIntoBatches(batches[[1]], batch = batch,
                               restrict = restrict[[1]]) 
  batches <- divided$batches 
  restrict <- divided$restrict 
}

  nbatches <- length(batches)
  if (nbatches < 2L) {
    stop("at least two batches must be specified")
  }
  
  library(scuttle)
  .apply_cosine_norm <- function(x, l2) {
    l2 <- pmax(1e-8, l2) # protect against zero-L2.
    normalizeCounts(x, size_factors=l2, center_size_factors=FALSE, log=FALSE)
  }
  
  .prepare_input_data <- function (batches, cos.norm.in, cos.norm.out, subset.row, correct.all) 
  {
    nbatches <- length(batches)
    in.batches <- out.batches <- batches
    same.set <- TRUE
    if (!is.null(subset.row)) {
      subset.row <- .row_subset_to_index(batches[[1]], subset.row)
      if (identical(subset.row, seq_len(nrow(batches[[1]])))) {
        subset.row <- NULL
      }
      else {
        in.batches <- lapply(in.batches, "[", i = subset.row, 
                             , drop = FALSE)
        if (correct.all) {
          same.set <- FALSE
        }
        else {
          out.batches <- in.batches
        }
      }
    }
    if (cos.norm.in) {
      norm.scaling <- vector("list", nbatches)
      for (b in seq_len(nbatches)) {
        current.in <- in.batches[[b]]
        cos.out <- cosineNorm(current.in, mode = "all")
        in.batches[[b]] <- cos.out$matrix
        norm.scaling[[b]] <- cos.out$l2norm
      }
    }
    if (cos.norm.out) {
      if (!cos.norm.in) {
        norm.scaling <- lapply(in.batches, cosineNorm, mode = "l2norm")
      }
      out.batches <- mapply(FUN = .apply_cosine_norm, out.batches, 
                            norm.scaling, SIMPLIFY = FALSE)
    }
    if (cos.norm.out != cos.norm.in) {
      same.set <- FALSE
    }
    list(In = in.batches, Out = out.batches, Subset = subset.row, 
         Same = same.set)
  }
  
  prep.out <- .prepare_input_data(batches, cos.norm.in = cos.norm.in,
                                  cos.norm.out = cos.norm.out, subset.row = subset.row,
                                  correct.all = correct.all)
  in.batches <- prep.out$In # 3个100(基因)*150(细胞)
  out.batches <- prep.out$Out
  subset.row <- prep.out$Subset
  same.set <- prep.out$Same
  in.batches <- lapply(in.batches, t) # 3个150(细胞)*100(基因)
  if (!same.set) {
    out.batches <- lapply(out.batches, t)
  }
  if (!auto.merge) {
    
   
    .create_tree_predefined <-  function (batches, restrict, merge.order) 
    {
      if (is.null(merge.order)) {
        merge.order <- seq_along(batches)
      }
      if (!is.list(merge.order) && length(merge.order) > 1L) {
        merge.tree <- list(merge.order[1], merge.order[2])
        for (i in tail(merge.order, -2L)) {
          merge.tree <- list(merge.tree, i)
        }
      }
      else {
        merge.tree <- merge.order
      }
      
      .binarize_tree <- function (merge.tree) 
      {
        if (!is.list(merge.tree) && length(merge.tree) == 1L) {
          if (is.factor(merge.tree)) {
            merge.tree <- as.character(merge.tree)
          }
          return(merge.tree)
        }
        N <- length(merge.tree)
        if (N == 1L) {
          .binarize_tree(merge.tree[[1]])
        }
        else if (N == 2L) {
          list(.binarize_tree(merge.tree[[1]]), .binarize_tree(merge.tree[[2]]))
        }
        else if (N > 2L) {
          current <- list(.binarize_tree(merge.tree[[1]]), .binarize_tree(merge.tree[[2]]))
          for (i in 3:N) {
            current <- list(current, .binarize_tree(merge.tree[[i]]))
          }
          current
        }
        else {
          stop("merge tree contains a node with no children")
        }
      }
      
      merge.tree <- .binarize_tree(merge.tree)
      leaves <- unlist(merge.tree)
      if (!is.numeric(leaves)) {
        leaves <- match(as.character(leaves), names(batches))
      }
      else {
        leaves <- as.integer(leaves)
      }
      if (any(is.na(leaves)) || anyDuplicated(leaves) || any(leaves < 
                                                             1) || any(leaves > length(batches))) {
        stop("invalid leaf nodes specified in 'merge.order'")
      }
      merge.tree <- relist(leaves, merge.tree)
      
      .fill_tree <- function (merge.tree, batches, restrict) 
      {
        if (!is.list(merge.tree)) {
          val <- batches[[merge.tree]]
          
          MNN_treenode <-  function (index, data, restrict, origin = rep(index, nrow(data)), 
                    extras = list()) 
          {
            new("MNN_treenode", index = index, data = data, restrict = restrict, 
                origin = origin, extras = extras)
          }
          
          return(MNN_treenode(index = merge.tree, data = val, 
                              restrict = restrict[[merge.tree]]))
        }
        if (length(merge.tree) != 2L) {
          stop("merge tree structure should contain two children per node")
        }
        merge.tree[[1]] <- .fill_tree(merge.tree[[1]], batches, 
                                      restrict)
        merge.tree[[2]] <- .fill_tree(merge.tree[[2]], batches, 
                                      restrict)
        merge.tree
      }
      
      .fill_tree(merge.tree, batches, restrict)
    }
    
    merge.tree <- .create_tree_predefined(in.batches, restrict,
                                          merge.order)
    
    .get_next_merge <- function(merge.tree, path=NULL) {
      if (!is.list(merge.tree[[1]]) && !is.list(merge.tree[[2]])) {
        list(left=merge.tree[[1]], right=merge.tree[[2]], chosen=path)
      } else if (is.list(merge.tree[[2]])) {
        .get_next_merge(merge.tree[[2]], c(path, 2L))
      } else {
        .get_next_merge(merge.tree[[1]], c(path, 1L))
      }
    }
    
    
    .update_tree <- function(merge.tree, path, ...) {
      if (length(path)==0L) {
        return(MNN_treenode(...))
      }
      merge.tree[[path[1]]] <- .update_tree(merge.tree[[path[[1]]]], path[-1], ...)
      merge.tree
    }
    
    NEXT <- .get_next_merge
    UPDATE <- .update_tree
    
    .add_out_batches_to_tree <- function (merge.tree, out.batches) 
    {
      .get_node_index <- function(node) node@index
      if (!is.list(merge.tree)) {
        merge.tree@extras <- list(out.batches[[.get_node_index(merge.tree)]])
        return(merge.tree)
      }
      merge.tree[[1]] <- .add_out_batches_to_tree(merge.tree[[1]], 
                                                  out.batches)
      merge.tree[[2]] <- .add_out_batches_to_tree(merge.tree[[2]], 
                                                  out.batches)
      merge.tree
    }
    
    merge.tree <- .add_out_batches_to_tree(merge.tree, if (same.set)
      NULL
      else out.batches)
  }

