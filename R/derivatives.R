#------------------------------------------------------------------------------
# First derivative matrices wrt unstructured random effects covariances
#------------------------------------------------------------------------------

dV_dTau_index <- function(tau_index, Z_blocks) {
  lapply(Z_blocks, function(Z)
    Z[,tau_index, drop=FALSE] %*% t(Z)[rev(tau_index),,drop=FALSE])
}

dV_dTau_unstruct <- function(block, Z_design) {
  Tau_q <- dim(Z_design)[2]
  Z_blocks <- by(Z_design, block, as.matrix)
  tau_index <- cbind(unlist(sapply(1:Tau_q, function(x) seq(1,x))),
                     unlist(sapply(1:Tau_q, function(x) rep(x,x))))
  apply(tau_index, 1, function(t) dV_dTau_index(unique(t), Z_blocks = Z_blocks))
}

dV_dreStruct <- function(mod) {

  blocks <- mod$groups
  Z_design <- model.matrix(mod$modelStruct$reStruct, data = mod$data)

  if (length(blocks) == 1L) {
    Z_list <- list(Z_design)
  } else {
    Z_list <- sapply(names(blocks),
                     function(x) Z_design[,grep(x, colnames(Z_design)), drop = FALSE],
                     simplify = FALSE, USE.NAMES = TRUE)
  }

  mapply(dV_dTau_unstruct,
         block = blocks, Z_design = Z_list,
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
  cor_AR1 <- as.double(coef(mod$modelStruct$corStruct, FALSE))
  covariate <- attr(struct, "covariate")
  lapply(covariate, function(x) as.matrix(dist(x)) * cor_AR1^(as.matrix(dist(x)) - 1L))
}

# corCAR1

dR_dcorStruct.corCAR1 <- function(struct) {
  cor_CAR1 <- as.double(coef(mod$modelStruct$corStruct, FALSE))
  covariate <- attr(struct, "covariate")
  lapply(covariate, function(x) as.matrix(dist(x)) * cor_CAR1^(as.matrix(dist(x)) - 1L))
}

# corCompSymm

dR_dcorStruct.corCompSymm <- function(struct) {
  covariate <- attr(struct, "covariate")
  lapply(covariate, function(x) 1L - diag(1L, nrow = length(x)))
}

# corSymm

# return list of matrices with derivative for one parameter

dR_dcor_index <- function(row, col, dim_vec, q) {
  R_mat <- matrix(0, q, q)
  R_mat[row, col] <- 1
  R_mat[col, row] <- 1
  replicate(dim_vec, R_mat, simplify=FALSE)
}

# return a list of block diagonal matrix for all parameters

dR_dcorStruct.corSymm <- function(struct) {
  cor_Symm <- as.double(coef(struct, FALSE)) # parameters

  # max dimension of the correlation matrices
  cor_q <- (1 + sqrt(1 + 8 * length(cor_Symm))) / 2 # number of cols of R matrix

  # the index for which we are taking derivatives wrt
  cor_index <- as.matrix(t(combn(1:cor_q, 2))) # sort cor_index appropriately

  # get dimensions of block-diagonal matrices
  grps <- attr(struct, "groups")
  dim_vec <- length(unique(grps)) # dim of the block diagonal matrix

  apply(cor_index, 1, function(t) dR_dcor_index(t[1], t[2], dim_vec, cor_q))
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
}

dV_dvarStruct <- function(mod) {

  # No derivatives if there's no variance structure or only varFixed structure
  if (is.null(mod$modelStruct$varStruct) | inherits(mod$modelStruct$varStruct,"varFixed")) return(NULL)

  dsd_dvar <- dsd_dvarStruct(mod$modelStruct$varStruct)

  R_list <- build_corr_mats(mod)

  all_groups <- rev(mod$groups)
  wts <- nlme::varWeights(mod$modelStruct$varStruct)[order(do.call(order, all_groups))]

  if (is.null(R_list)) {

    dV_dvar <- lapply(dsd_dvar, function(d) 2 * d * mod$sigma / wts)

    dV_list <- lapply(dV_dvar, function(v) {
      v_list <- split(v, f = mod$groups[[1]])
      attr(v_list, "groups") <- "diagonal"
      v_list
    })

  } else {

    sd_list <- split(mod$sigma / wts, attr(R_list, "groups"))
    dsd_list <- lapply(dsd_dvar, split, f = mod$groups[[1]])
    dV_list <- lapply(dsd_list, sdRds, sd_list = sd_list, R_list = R_list)

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

dsd_dvarExp <- function(val, grp, strt = struct) {
  grps <- attr(strt, "groups")
  covariate <- as.numeric(attr(strt, "covariate"))
  exp(covariate * val) * covariate * as.integer(grp == grps)
}

dsd_dvarStruct.varExp <- function(struct) {
  var_Exp <- coef(struct, FALSE) # get the parameter
  par_val <- as.double(var_Exp)
  par_name <- attr(var_Exp, "names")

  if (length(var_Exp) == 1) {
    covariate <- attr(struct, "covariate")
    exp(covariate * par_val) * covariate
  } else{
    Map(dsd_dvarExp, val = par_val, grp = par_name)
  }
}

# varPower

dsd_dvarPower <- function(val, grp, strt = struct) {
  grps <- attr(strt, "groups")
  covariate <- as.numeric(attr(strt, "covariate"))
  abs_covariate <- abs(covariate)
  abs_covariate^val * log(abs_covariate) * as.integer(grp == grps)
}

dsd_dvarStruct.varPower <- function(struct) {
  var_Power <- as.list(coef(struct, FALSE))
  par_val <- as.double(var_Power)
  par_name <- attr(var_Power, "names")

  if (length(var_Power) == 1) {
    covariate <- attr(struct, "covariate")
    abs_covariate <- abs(covariate)
    abs_covariate^par_val * log(abs_covariate)
  } else{
    Map(dsd_dvarPower, val = par_val, grp = par_name)
  }
}

# varConstPower

# for one stratum (two parameters: const and power)

dsd_dConstPower1 <- function(x, strt = struct) {
  var_ConstPower <- coef(strt, FALSE)
  par_val <- as.double(var_ConstPower)
  covariate <- attr(strt, "covariate")
  abs_covariate <- abs(covariate)

  if (x == "const") {
    res <- rep(1, length(covariate))
  } else {
    res <- abs_covariate^par_val[2] * log(abs_covariate)
  }
  res
}

# for two or more strata (multiple const and power parameters)

dsd_dConstPower2 <- function(val, type, grp, strt = struct) {
  grps <- attr(strt, "groups")
  covariate <- as.numeric(attr(strt, "covariate"))
  abs_covariate <- abs(covariate)

  if (type == "const") {
    as.integer(grp == grps)
  } else {
    abs_covariate^val * log(abs_covariate) * as.integer(grp == grps)
  }
}

dsd_dvarStruct.varConstPower <- function(struct) {

  # get the var struct
  var_ConstPower <- coef(struct, FALSE)

  # get the var struct names
  par_name <- names(var_ConstPower)

  # get the par values
  par_val <- as.double(var_ConstPower)

  # indicates whether the par is a const or power, can be used for one stratum scenario
  par_type <- substr(par_name, 1, 5)

  # indicates the stratum of the par, cannot be used for one stratum scenario
  par_grp <- substring(par_name, 7)

  if (length(var_ConstPower) == 2) {

    lapply(par_name, dsd_dConstPower1)

  } else {

    Map(dsd_dConstPower2, val = par_val, type = par_type, grp = par_grp)

  }

}
