#------------------------------------------------------------------------------
# Build list of derivative matrices wrt model parameters
#------------------------------------------------------------------------------

build_dV_list <- function(mod) UseMethod("build_dV_list")

build_dV_list.default <- function(mod) {
  mod_class <- paste(class(mod), collapse = "-")
  stop(paste0("Derivatives not available for models of class ", mod_class, "."))

}

build_dV_list.gls <- function(mod) {
  cor_params <- dV_dcorStruct(mod)                          # correlation structure
  var_params <- dV_dvarStruct(mod)                          # variance structure
  sigma_sq <- dV_dsigmasq(mod)                              # sigma_sq

  # Create a list of derivative matrices
  c(cor_params, var_params, sigma_sq)
}

build_dV_list.lme <- function(mod) {
  Tau_params <- dV_dreStruct(mod)                           # random effects structure(s)
  cor_params <- dV_dcorStruct(mod)                          # correlation structure
  var_params <- dV_dvarStruct(mod)                          # variance structure
  sigma_sq <- dV_dsigmasq(mod)                              # sigma_sq

  # Create a list of derivative matrices
  c(unlist(Tau_params, recursive = FALSE), cor_params, var_params, sigma_sq)
}



#------------------------------------------------------------------------------
# First derivative matrices wrt unstructured random effects covariances
#------------------------------------------------------------------------------

dV_dTau_index <- function(tau_index, Z_blocks, block) {
  dV_dTau <- lapply(Z_blocks, function(Z) Z[,tau_index, drop=FALSE] %*% t(Z)[rev(tau_index),,drop=FALSE])
  attr(dV_dTau, "groups") <- block
  return(dV_dTau)
}

dV_dTau_unstruct <- function(block, pdMat_class, Z_design) {
  Tau_q <- dim(Z_design)[2]
  Z_blocks <- by(Z_design, block, as.matrix, simplify = FALSE)

  if ("pdDiag" %in% pdMat_class) {
    tau_index <- cbind(seq(1,Tau_q), seq(1,Tau_q))
  } else if ("pdSymm" %in% pdMat_class) {
    tau_index <- cbind(unlist(sapply(1:Tau_q, function(x) seq(1,x))),
                       unlist(sapply(1:Tau_q, function(x) rep(x,x))))
  } else {
    stop("Tau_index only available for pdMat structures of class pdDiag and pdSymm.")
  }

  apply(tau_index, 1, function(t) dV_dTau_index(unique(t), Z_blocks = Z_blocks, block = block))
}

dV_dreStruct <- function(mod) {
  blocks <- mod$groups
  blocks_names <- names(blocks)
  b <- lapply(blocks_names, function(x) class(mod$modelStruct$reStruct[[x]])) # pdClass
  data <- mod$data
  Z_design <- model.matrix(mod$modelStruct$reStruct, data = data[complete.cases(data), ])

  if (length(blocks) == 1L) {
    Z_list <- list(Z_design)
  } else {
    Z_list <- sapply(names(blocks),
                     function(x) Z_design[,grep(x, colnames(Z_design)), drop = FALSE],
                     simplify = FALSE, USE.NAMES = TRUE)

  }

  mapply(dV_dTau_unstruct,
         block = blocks, pdMat_class = b, Z_design = Z_list,
         SIMPLIFY = FALSE)
}

#------------------------------------------------------------------------------
# First derivative matrices wrt correlation structures
#------------------------------------------------------------------------------

# Derivative of Sigma with respect to correlation structure

dV_dcorStruct <- function(mod) {

  # No derivatives if there's no correlation structure
  if (is.null(mod$modelStruct$corStruct)) return(NULL)

  dR_dcor <- dR_dcorStruct(mod$modelStruct$corStruct)

  grps <- get_cor_grouping(mod)

  dR_dcor <- lapply(dR_dcor, function(x) {
    x_grps <- factor(grps, levels = names(x))
    attr(x, "groups") <- x_grps
    x
  })

  lapply(dR_dcor, build_var_cor_mats, mod = mod, sigma_scale = TRUE)
}

# Methods for different corStruct classes

dR_dcorStruct <- function(struct) UseMethod("dR_dcorStruct")

dR_dcorStruct.default <- function(struct) {
  cor_class <- class(struct)[[1]]
  stop(paste0("Derivatives not available for correlation structures of class ", cor_class, "."))
}

# corAR1

dR_dcorStruct.corAR1 <- function(struct) {
  cor_AR1 <- as.double(coef(struct, FALSE))
  covariate <- attr(struct, "covariate")
  if (!is.list(covariate)) covariate <- list(A = covariate)
  dR <- lapply(covariate, function(x) as.matrix(dist(x)) * cor_AR1^(as.matrix(dist(x)) - 1L))
  list(dR)
}

# corCAR1

dR_dcorStruct.corCAR1 <- function(struct) {
  cor_CAR1 <- as.double(coef(struct, FALSE))
  covariate <- attr(struct, "covariate")
  if (!is.list(covariate)) covariate <- list(A = covariate)
  dR <- lapply(covariate, function(x) as.matrix(dist(x)) * cor_CAR1^(as.matrix(dist(x)) - 1L))
  list(dR)
}

# corARMA

dR_dcorMA1 <- function(covariate, cor) {
  dist_mat <- as.matrix(dist(covariate))
  dist_mat[dist_mat != 1] <- 0
  dist_mat[dist_mat == 1] <- (1 - cor^2) / ((1 + cor^2)^2)
  return(dist_mat)
}

dR_dcorStruct.corMA1 <- function(struct) {
  cor_MA1 <- as.double(coef(struct, FALSE))
  covariate <- attr(struct, "covariate")
  if (!is.list(covariate)) covariate <- list(A = covariate)
  dR <- lapply(covariate, dR_dcorMA1, cor = cor_MA1)
  list(dR)
}

dR_dcorStruct.corARMA <- function(struct) {
  cor_ARMA <- coef(struct, FALSE)
  cor_name <- names(cor_ARMA)
  p <- length(grep("Phi", cor_name, value = TRUE))
  q <- length(grep("Theta", cor_name, value = TRUE))

  if (p == 1 & q == 0) {
    dR_dcorStruct.corAR1(struct)
  } else if (p == 0 & q == 1) {
    dR_dcorStruct.corMA1(struct)
  } else {
    stop("Derivatives not available for correlation structures of class corARMA, except for the simple cases of AR(1) or MA(1).")
  }
}

# corCompSymm

dR_dcorStruct.corCompSymm <- function(struct) {
  covariate <- attr(struct, "covariate")
  if (!is.list(covariate)) covariate <- list(A = covariate)
  dR <- lapply(covariate, function(x) 1L - diag(1L, nrow = length(x)))
  list(dR)
}

# corSymm

# return list of derivative matrices for one cor parameter

replace <- function(x, y, row, col) {
  row_index <- y == row
  col_index <- y == col
  x[row_index, col_index] <- x[col_index, row_index] <- 1L
  x
}

dR_dcor_index <- function(row, col, covariate) {
  R_null <- lapply(covariate, function(x) matrix(0L, nrow = length(x), ncol = length(x)))
  dR <- mapply(replace, x = R_null, y = covariate, MoreArgs = list(row = row, col = col), SIMPLIFY = FALSE)
  dR
}

# return a list of derivative matrices for all cor parameters

dR_dcorStruct.corSymm <- function(struct) {
  cor_Symm <- as.double(coef(struct, FALSE)) # parameters
  cor_q <- (1 + sqrt(1 + 8 * length(cor_Symm))) / 2 # number of cols/rows of R mat
  cor_index <- as.matrix(t(utils::combn(1:cor_q, 2))) # sort cor_index appropriately
  groups <- attr(struct, "groups")
  covariate <- attr(struct, "covariate")
  covariate <- lapply(covariate, function(x) as.integer(x + 1)) # add 1 to covariate to align with cor_index
  apply(cor_index, 1, function(t) dR_dcor_index(t[1], t[2], covariate = covariate))
}

#------------------------------------------------------------------------------
# First derivative matrices wrt sigma^2
#------------------------------------------------------------------------------

dV_dsigmasq <- function(mod) {

  fixed_sigma <- attr(mod$modelStruct, "fixedSigma")

  if (fixed_sigma) {
    NULL
  } else {
    list(build_var_cor_mats(mod, sigma_scale = FALSE))
  }

}

#------------------------------------------------------------------------------
# First derivative matrices wrt variance structures
#------------------------------------------------------------------------------

# Derivative of Sigma with respect to variance structure parameters

sdRds <- function(dsd_list, sd_list, R_list) {
  dsd_list <- Map(function(s, d) tcrossprod(s, d),
             s = sd_list, d = dsd_list)
  dV_list <- Map(function(dsd, R) (dsd + t(dsd)) * R,
                 dsd = dsd_list, R = R_list)
  attr(dV_list, "groups") <- attr(R_list, "groups")
  dV_list
}

dV_dvarStruct <- function(mod) {

  # No derivatives if there's no variance structure or only varFixed structure
  if (is.null(mod$modelStruct$varStruct) | inherits(mod$modelStruct$varStruct,"varFixed")) return(NULL)

  # all_groups <- rev(mod$groups)
  # sort_order <- order(do.call(order, all_groups))
  sort_order <- get_sort_order(mod)

  dsd_dvar <- dsd_dvarStruct(mod$modelStruct$varStruct)

  dsd_dvar <- lapply(dsd_dvar, function(x) x[sort_order]) # reorder based on input data

  R_list <- build_corr_mats(mod)

  wts <- nlme::varWeights(mod$modelStruct$varStruct)[sort_order]

  if (is.null(R_list)) {

    dV_dvar <- lapply(dsd_dvar, function(d) 2 * d * mod$sigma^2 / wts)

    groups <- if (is.null(mod$groups)) factor(1:mod$dims$N) else mod$groups[[1]]

    dV_list <- lapply(dV_dvar, function(v) {
      v_list <- tapply(v, groups, function(x) diag(x, nrow = length(x)), simplify = FALSE)
      attr(v_list, "groups") <- groups
      v_list
    })

  } else {
    grps <- attr(R_list, "groups")
    sigmasq_S_list <- split(mod$sigma^2 / wts, grps)
    dsd_list <- lapply(dsd_dvar, split, f = grps)
    dV_list <- lapply(dsd_list, sdRds, sd_list = sigmasq_S_list, R_list = R_list)

  }

  dV_list

}

# Methods for different varStruct classes

dsd_dvarStruct <- function(struct) UseMethod("dsd_dvarStruct")

dsd_dvarStruct.default <- function(struct) {
  var_class <- class(struct)[[1]]
  stop(paste0("Derivatives not available for variance structures of class ", var_class, "."))
}

# varIdent

dsd_dvarStruct.varIdent <- function(struct) {
  grps <- attr(struct, "groups")
  par_name <- names(coef(struct, FALSE)) # sd parameter name
  lapply(par_name, function(x) as.integer(x == grps))
}

# varExp

dsd_dvarExp <- function(val, grp, groups, covariate) {
  exp(covariate * val) * covariate * as.integer(grp == groups)
}

dsd_dvarStruct.varExp <- function(struct) {
  var_Exp <- coef(struct, FALSE)
  par_val <- as.double(var_Exp)
  par_name <- attr(var_Exp, "names")
  groups <- attr(struct, "groups")
  covariate <- attr(struct, "covariate")

  if (length(var_Exp) == 1) {

    covariate <- attr(struct, "covariate")
    list(exp(covariate * par_val) * covariate)

  } else{

    mapply(dsd_dvarExp, val = par_val, grp = par_name,
           MoreArgs = list(groups = groups, covariate = covariate),
           SIMPLIFY = FALSE)
  }
}

# varPower

dsd_dvarPower <- function(val, grp, groups, covariate) {
  abs_covariate <- abs(covariate)
  abs_covariate^val * log(abs_covariate) * as.integer(grp == groups)
}

dsd_dvarStruct.varPower <- function(struct) {
  var_Power <- coef(struct, FALSE)
  par_val <- as.double(var_Power)
  par_name <- attr(var_Power, "names")
  groups <- attr(struct, "groups")
  covariate <- attr(struct, "covariate")

  if (length(var_Power) == 1) {

    abs_covariate <- abs(covariate)
    list(abs_covariate^par_val * log(abs_covariate))

  } else {

    mapply(dsd_dvarPower, val = par_val, grp = par_name,
           MoreArgs = list(groups = groups, covariate = covariate),
           SIMPLIFY = FALSE)

  }
}

# varConstPower

# for one stratum (two parameters: const and power)

dsd_dConstPower1 <- function(x, val, covariate) {
  abs_covariate <- abs(covariate)
  if (x == "const") {
    rep(1, length(covariate))
  } else {
    abs_covariate^val[2] * log(abs_covariate) # [2] is the power par
  }
}

# for two or more strata (multiple const and power parameters)

dsd_dConstPower2 <- function(val, type, grp, groups, covariate) {
  abs_covariate <- abs(covariate)
  if (type == "const") {
    as.integer(grp == groups)
  } else {
    abs_covariate^val * log(abs_covariate) * as.integer(grp == groups)
  }
}

dsd_dvarStruct.varConstPower <- function(struct) {
  groups <- attr(struct, "groups")
  covariate <- attr(struct, "covariate")
  var_ConstPower <- coef(struct, FALSE)
  par_name <- names(var_ConstPower)
  par_val <- as.double(var_ConstPower)

  # indicates whether the par is a const or power, can be used for one stratum scenario
  par_type <- substr(par_name, 1, 5)

  # indicates the stratum of the par, cannot be used for one stratum scenario
  par_grp <- substring(par_name, 7)

  if (length(var_ConstPower) == 2) {

    lapply(par_name, dsd_dConstPower1, val = par_val, covariate = covariate)

  } else {

    mapply(dsd_dConstPower2, val = par_val, type = par_type, grp = par_grp,
           MoreArgs = list(groups = groups, covariate = covariate),
           SIMPLIFY = FALSE)

  }
}
