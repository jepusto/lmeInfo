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
  Z_list <- sapply(names(blocks),
                   function(x) Z_design[,grep(x, colnames(Z_design)), drop = FALSE],
                   simplify = FALSE, USE.NAMES = TRUE)
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

}

# corCAR1

dR_dcorStruct.corAR1 <- function(struct) {

}

# corCompSymm

dR_dcorStruct.corCompSymm <- function(struct) {

}

# corCompSymm

dR_dcorStruct.corSymm <- function(struct) {

}

#------------------------------------------------------------------------------
# First derivative matrices wrt variance structures
#------------------------------------------------------------------------------

# Derivative of Sigma with respect to variance structure parameters

dV_dvarStruct <- function(mod) {

  # No derivatives if there's no variance structure or only varFixed structure
  if (is.null(mod$modelStruct$varStruct) | inherits(mod$modelStruct$varStruct,"varFixed")) return(NULL)

  dsd_dvar <- dsd_dvarStruct(mod$modelStruct$varStruct)

  R_list <- build_corr_mats(mod)

  all_groups <- rev(mod$groups)
  wts <- nlme::varWeights(mod$modelStruct$varStruct)[order(do.call(order, all_groups))]

  if (is.null(R_list)) {

    sd_list <- split(mod$sigma / wts, mod$groups[[1]])
    dV_list <- Map(function(s, d) 2 * s * d,
                   s = sd_list, d = dsd_dvar)
    attr(dV_list, "groups") <- "diagonal"

  } else {

    sd_list <- split(mod$sigma / wts, attr(R_list, "groups"))
    dV_list <- Map(function(R, s, d) (tcrossprod(s, d) + tcrossprod(d, s)) * R,
                  R = R_list, s = sd_list, d = dsd_dvar)
    attr(dV_list, "groups") <- attr(R_list, "groups")
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

}

# varExp

dsd_dvarStruct.varExp <- function(struct) {

}

# varPower

dsd_dvarStruct.varPower <- function(struct) {

}

# varConstPower

dsd_dvarStruct.varConstPower <- function(struct) {

}
