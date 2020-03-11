


##------------------------------------------------------------------------
## block-diagonal matrix addition, multiplication, and trace functions
##------------------------------------------------------------------------

# turn matrix into a list of sub-matrices

sub_f <- function(x, fac, dim) {
  function(f) switch(dim,
                     row = x[fac==f, ,drop=FALSE],
                     col = x[ ,fac==f, drop=FALSE],
                     both = x[fac==f, fac==f, drop=FALSE])
}

matrix_list <- function(x, fac, dim) {
  if (is.vector(x)) {
    if (dim != "both") stop(paste0("Object must be a matrix in order to subset by ",dim,"."))
    x_list <- split(x, fac)
    lapply(x_list, function(x) diag(x, nrow = length(x)))
  } else {
    lapply(levels(fac), sub_f(x, fac, dim))
  }
}

# turn block-diagonal into regular matrix

unblock <- function(A, block=attr(A, "groups")) {

  if (is.null(block)) block <- factor(rep(names(A), times = sapply(A, function(x) dim(x)[1])))
  n <- length(block)
  mat <- matrix(0, n, n)
  for (i in levels(block)) {
    index <- i == block
    mat[index,index] <- A[[i]]
  }
  return(mat)
}


# sum of two conformable block-diagonal matrices

sum_blockblock <- function(A, B)
  mapply(function(a,b) a + b, a = A, b = B, SIMPLIFY = FALSE)


# generic matrix minus block-diagonal

matrix_minus_block <- function(A, B, block=NULL) {
  if (is.null(block)) block <- rep(names(B), times = sapply(B, function(x) dim(x)[1]))

  mat <- A
  for (i in unique(block)) {
    index <- i == block
    mat[index,index] <- mat[index, index] - B[[i]]
  }
  return(mat)
}


# block-diagonal minus generic matrix

block_minus_matrix <- function(A, B, block=NULL) {
  if (is.null(block))
    block <- rep(names(A), times = sapply(A, function(x) dim(x)[1]))

  mat <- -B
  for (i in unique(block)) {
    index <- i == block
    mat[index,index] <- mat[index, index] + A[[i]]
  }
  return(mat)
}

add_submatrices <- function(indices, small_mat, big_mat) {
  levs <- levels(indices)
  if (nlevels(indices) != length(small_mat)) stop("Levels of indices do not match entries of small_mat.")
  for (i in 1:length(levs)) {
    ind <- levs[i] == indices
    big_mat[ind,ind] <- big_mat[ind,ind] + small_mat[[i]]
  }
  big_mat
}

add_bdiag <- function(small_mats, big_mats, crosswalk) {
  small_indices <- lapply(split(crosswalk[[1]], crosswalk[[2]]), droplevels)
  big_indices <- unique(crosswalk)
  big_indices <- big_indices[[2]][order(big_indices[[1]])]
  small_mats_list <- split(small_mats, big_indices)
  Map(add_submatrices, indices = small_indices, small_mat = small_mats_list, big_mat = big_mats)
}

# sum of conformable diagonal matrix and block-diagonal matrix

add_diag <- function(d, M) {
  diag(M) <- diag(M) + d
  M
}

add_diag_bdiag <- function(diag_mats, big_mats) {
  Map(add_diag, d = diag_mats, M = big_mats)
}

# product of two block-diagonal matrices

prod_blockblock <- function(A, B, crosswalk = NULL) {

  if (is.null(crosswalk)) {
    A_groups <- attr(A, "groups")
    B_groups <- attr(B, "groups")
    if (is.null(A_groups) | is.null(B_groups)) {
      stop("Must specify a crosswalk or use matrices with groups attribute.")
    }
  } else {
    A_groups <- crosswalk[[1]]
    B_groups <- crosswalk[[2]]
  }

  B_in_A <- all(tapply(A_groups, B_groups, function(x) length(unique(x)) == 1))
  A_in_B <- all(tapply(B_groups, A_groups, function(x) length(unique(x)) == 1))

  if (B_in_A & A_in_B) {
    res <- mapply(function(a, b) a %*% b, a = A, b = B, SIMPLIFY = FALSE)
  } else if (B_in_A) {
    # B is nested in A
    if (is.null(names(A))) names(A) <- levels(A_groups)
    block_list <- split(B_groups, A_groups)
    B_map <- tapply(levels(A_groups)[A_groups], B_groups, unique)[names(B)]
    B_list <- split(B, B_map)[names(A)]
    res <- mapply(prod_matrixblock, A = A, B = B_list, block = block_list)
  } else if (A_in_B) {
    # A is nested in B
    if (is.null(names(B))) names(B) <- levels(B_groups)
    block_list <- split(A_groups, B_groups)
    A_map <- tapply(levels(B_groups)[B_groups], A_groups, unique)[names(A)]
    A_list <- split(A, A_map)[names(B)]
    res <- mapply(prod_blockmatrix, A = A_list, B = B, block = block_list)
  } else {
    stop("The A and B matrices are not nested.")
  }

  if (!is.null(attr(A, "groups"))) {
    attr(res, "groups") <- A_groups
  }

  return(res)
}



# product of a block-diagonal matrix and a generic matrix

prod_blockmatrix <- function(A, B, block = attr(A, "groups")) {

  if (is.null(names(A))) names(A) <- 1:length(A)
  A_names <- names(A)

  if (is.null(block)) block <- rep(A_names, times = sapply(A, function(x) dim(x)[1]))

  C <- matrix(0, length(block), dim(B)[2])

  for (b in A_names) {
    ind <- block == b
    C[ind, ] <- A[[b]] %*% B[ind,]
  }
  return(C)
}


# product of a generic matrix and a block-diagonal matrix

prod_matrixblock <- function(A, B, block = attr(B, "groups")) {

  if (is.null(names(B))) names(B) <- 1:length(B)
  B_names <- names(B)

  if (is.null(block)) block <- rep(B_names, times = sapply(B, function(x) dim(x)[2]))

  C <- matrix(0, dim(A)[1], length(block))

  for (b in B_names) {
    ind <- block == b
    C[,ind] <- A[,ind] %*% B[[b]]
  }
  return(C)
}


# trace of the product of two generic matrices

product_trace <- function(A,B) sum(as.vector(t(A)) * as.vector(B))


# trace of the product of two conformable block-diagonal matrices

product_trace_blockblock <- function(A, B)
  sum(mapply(function(a, b) product_trace(a,b), a = A, b = B))


