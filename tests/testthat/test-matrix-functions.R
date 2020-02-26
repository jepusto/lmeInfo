set.seed(20200225)
N <- 200
fac <- factor(sample(LETTERS, size = N, replace = TRUE))
gac <- sample(LETTERS[3:8], size = nlevels(fac), replace = TRUE)
gac <- factor(gac[fac], levels = LETTERS)
gac_drop <- droplevels(gac)

X <- matrix(rnorm(N * 4), nrow = N)
Z <- matrix(rnorm(N * 2), nrow = N)

random_sym_matrix <- function(dim, df)
  rWishart(1, df = df, Sigma = diag(1, nrow = dim))[,,1]

U_list <- lapply(table(fac), random_sym_matrix, df = max(table(fac)))
V_list <- lapply(table(fac), random_sym_matrix, df = max(table(fac)))
W_list <- lapply(table(gac_drop), random_sym_matrix, df = max(table(gac)))
H_list <- lapply(table(gac_drop), random_sym_matrix, df = max(table(gac)))

U_full <- unblock(U_list, block = fac)
V_full <- unblock(V_list, block = fac)
W_full <- unblock(W_list, block = gac)
H_full <- unblock(H_list, block = gac)

X_list_fac <- matrix_list(X, fac = fac, dim = "row")
X_list_gac <- matrix_list(X, fac = gac_drop, dim = "row")
Z_list_fac <- matrix_list(Z, fac = fac, dim = "row")
Z_list_gac <- matrix_list(Z, fac = gac_drop, dim = "row")

test_that("unblock() and matrix_list() work.", {

  for (i in levels(fac)) {
    expect_identical(V_full[fac == i, fac == i,drop=FALSE], V_list[[i]])
  }
  V_relist <- matrix_list(V_full, fac = fac, dim = "both")
  names(V_relist) <- levels(fac)
  expect_identical(V_list, V_relist)

  for (i in levels(gac_drop)) {
    expect_identical(W_full[gac == i, gac == i,drop=FALSE], W_list[[i]])
  }
  W_relist <- matrix_list(W_full, fac = gac_drop, dim = "both")
  names(W_relist) <- levels(gac_drop)
  expect_identical(W_list, W_relist)

  XVX <- Map(function(x, v) t(x) %*% v %*% x, x = X_list_fac, v = V_list)
  XVX <- Reduce("+", XVX)
  XVX_full <- t(X) %*% V_full %*% X
  expect_equal(XVX, XVX_full)

  XWX <- Map(function(x, v) t(x) %*% v %*% x, x = X_list_gac, v = W_list)
  XWX <- Reduce("+", XWX)
  XWX_full <- t(X) %*% W_full %*% X
  expect_equal(XWX, XWX_full)

})

test_that("sum_blockblock(), add_submatrices(), and add_bdiag() work.", {

  U_V_list <- sum_blockblock(U_list, V_list)
  names(U_V_list) <- NULL
  U_V_full <- U_full + V_full
  expect_identical(U_V_list, matrix_list(U_V_full, fac = fac, dim = "both"))

  W_H_list <- sum_blockblock(W_list, H_list)
  names(W_H_list) <- NULL
  W_H_full <- W_full + H_full
  expect_identical(W_H_list, matrix_list(W_H_full, fac = gac_drop, dim = "both"))


})

test_that("add_diag() and add_diag_bdiag() work.", {

})


test_that("prod_blockmatrix(), prod_matrixblock(), and prod_blockblock() work.", {

})

test_that("product_trace() and product_trace_blockblock() work.", {

})
