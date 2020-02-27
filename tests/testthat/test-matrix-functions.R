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

  XXt <- tcrossprod(X)
  XXt_U <- add_submatrices(indices = fac, small_mat = U_list, big_mat = XXt)
  expect_identical(XXt_U, XXt + U_full)
  XXt_V <- add_submatrices(indices = fac, small_mat = V_list, big_mat = XXt)
  expect_identical(XXt_V, XXt + V_full)
  XXt_W <- add_submatrices(indices = gac_drop, small_mat = W_list, big_mat = XXt)
  expect_identical(XXt_W, XXt + W_full)
  XXt_H <- add_submatrices(indices = gac_drop, small_mat = H_list, big_mat = XXt)
  expect_identical(XXt_H, XXt + H_full)

  expect_error(add_submatrices(indices = gac, small_mat = W_list, big_mat = XXt))
  expect_error(add_submatrices(indices = gac, small_mat = H_list, big_mat = XXt))

  U_W_list <- add_bdiag(small_mats = U_list, big_mats = W_list,
                        crosswalk = data.frame(fac, gac_drop))
  expect_identical(unblock(U_W_list, block = gac_drop), U_full + W_full)

  expect_error(add_bdiag(small_mats = U_list, big_mats = W_list,
                         crosswalk = data.frame(fac, gac)))

})

test_that("add_diag_bdiag() work.", {
  D <- rnorm(N)
  D_fac <- split(D, fac)
  D_gac <- split(D, gac)
  D_gac_drop <- split(D, gac, drop = TRUE)

  D_U_list <- add_diag_bdiag(diag_mats = D_fac, big_mats = U_list)
  expect_identical(diag(D) + U_full, unblock(D_U_list, block = fac))
  D_W_list <- add_diag_bdiag(diag_mats = D_gac_drop, big_mats = W_list)
  expect_identical(diag(D) + W_full, unblock(D_W_list, block = gac_drop))

  expect_error(add_diag_bdiag(D_gac, W_list))
})


test_that("prod_blockmatrix(), prod_matrixblock(), and prod_blockblock() work.", {

  VZ <- prod_blockmatrix(A = V_list, B = Z, block = fac)
  expect_identical(VZ, V_full %*% Z)

  HZ <- prod_blockmatrix(A = H_list, B = Z, block = gac)
  HZ_drop <- prod_blockmatrix(A = H_list, B = Z, block = gac_drop)
  expect_identical(HZ, H_full %*% Z)
  expect_identical(HZ, HZ_drop)

  XtU <- prod_matrixblock(A = t(X), B = U_list, block = fac)
  expect_identical(XtU, t(X) %*% U_full)

  XtW <- prod_matrixblock(A = t(X), B = W_list, block = gac)
  XtW_drop <- prod_matrixblock(A = t(X), B = W_list, block = gac_drop)
  expect_identical(XtW, t(X) %*% W_full)
  expect_identical(XtW, XtW_drop)

  Zt_list_fac <- lapply(Z_list_fac, t)
  Zt_list_gac <- lapply(Z_list_gac, t)

  UX_fac <- prod_blockblock(A = U_list, B = X_list_fac, crosswalk = data.frame(fac, fac))
  UX_full <- Reduce(rbind, UX_fac)[order(order(fac)),]
  expect_identical(t(XtU), UX_full)

  ZtUX <- prod_blockblock(A = Zt_list_fac, B = UX_fac, crosswalk = data.frame(fac, fac))
  ZtUX <- Reduce("+", ZtUX)
  expect_equal(ZtUX, t(Z) %*% U_full %*% X)

  UX_gac <- prod_blockblock(A = U_list, B = X_list_gac, crosswalk = data.frame(fac, gac_drop))
  UX_full <- Reduce(rbind, UX_gac)[order(order(gac)),]
  expect_identical(t(XtU), UX_full)

  ZtUX <- prod_blockblock(A = Zt_list_gac, B = UX_gac, crosswalk = data.frame(gac_drop, gac_drop))
  ZtUX <- Reduce("+", ZtUX)
  expect_equal(ZtUX, t(Z) %*% U_full %*% X)

  WX_gac <- prod_blockblock(A = W_list, B = X_list_gac, crosswalk = data.frame(gac_drop, gac_drop))
  WX_full <- Reduce(rbind, WX_gac)[order(order(gac)),]
  expect_identical(t(XtW), WX_full)

  ZtWX <- prod_blockblock(A = Zt_list_gac, B = WX_gac,
                          crosswalk = data.frame(gac_drop, gac_drop))
  ZtWX <- Reduce("+", ZtWX)
  expect_equal(ZtWX, t(Z) %*% W_full %*% X)

  UV_list <- prod_blockblock(U_list, V_list, crosswalk = data.frame(fac,fac))
  expect_identical(unblock(UV_list, fac), U_full %*% V_full)

  UW_list <- prod_blockblock(U_list, W_list, crosswalk = data.frame(fac, gac_drop))
  UW_full <- U_full %*% W_full
  expect_identical(unblock(UW_list, gac_drop), UW_full)

  UWZ <- prod_blockmatrix(UW_list, Z, gac_drop)
  expect_identical(UWZ, UW_full %*% Z)

  HV_list <- prod_blockblock(H_list, V_list, crosswalk = data.frame(gac_drop, fac))
  HV_full <- H_full %*% V_full
  expect_identical(unblock(HV_list, gac_drop), HV_full)

  XtHV <- prod_matrixblock(t(X), HV_list, gac_drop)
  expect_identical(XtHV, t(X) %*% HV_full)

})

test_that("product_trace() and product_trace_blockblock() work.", {

  XtX_tr <- product_trace(t(X), X)
  expect_equal(XtX_tr, sum(diag(t(X) %*% X)))

  UZ <- prod_blockmatrix(U_list, Z, fac)
  VZ <- prod_blockmatrix(V_list, Z, fac)
  ZtUVZ_tr <- product_trace(t(UZ), VZ)
  expect_equal(ZtUVZ_tr, sum(diag(t(UZ) %*% VZ)))

  UZ_list <- prod_blockblock(U_list, Z_list_fac, crosswalk = data.frame(fac, fac))
  ZtUt_list <- lapply(UZ_list, t)
  VZ_list <- prod_blockblock(V_list, Z_list_fac, crosswalk = data.frame(fac, fac))
  ZtUVZ_tr_list <- product_trace_blockblock(ZtUt_list, VZ_list)
  expect_equal(ZtUVZ_tr, ZtUVZ_tr_list)
})
