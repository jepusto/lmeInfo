
# Construct list of block-diagonal correlation matrices

build_corr_mats <- function(mod) {

  if (is.null(mod$modelStruct$corStruct)) {
    return(NULL)
  } else {
    R_list <- nlme::corMatrix(mod$modelStruct$corStruct)
    # grps <- nlme::getGroups(mod, form = nlme::getGroupsFormula(mod$modelStruct$corStruct))
    grps <- stats::model.frame(nlme::getGroupsFormula(mod$modelStruct$corStruct), data = nlme::getData(mod))
    grps <- factor(apply(grps, 1, paste, collapse = "/"), levels = names(R_list))
    attr(R_list, "groups") <- grps
    return(R_list)
  }
}

# Construct list of block-diagonal lowest-level var-cov matrices

get_sort_order <- function(mod) {
  groups <- mod$groups
  if (is.data.frame(groups)) {
    order(do.call(order, groups))
  } else if (!is.null(groups)) {
    order(order(groups))
  } else {
    1:mod$dims$N
  }
}

build_var_cor_mats <- function(mod, R_list = build_corr_mats(mod), sigma_scale = FALSE) {

  sigma <- if (sigma_scale) mod$sigma else 1

  if (is.null(R_list)) {

    # if there is no correlation structure,
    # then build block-diagonals with first available grouping variable

    if (is.null(mod$groups)) {

      # if there are no groups then make diagonal matrix-lists

      if (is.null(mod$modelStruct$varStruct)) {
        V_list <- as.list(rep(sigma^2, mod$dims$N))
      } else {
        sd_vec <- sigma / as.numeric(nlme::varWeights(mod$modelStruct$varStruct))
        V_list <- as.list(sd_vec^2)
      }
      grps <- factor(1:mod$dims$N)
      attr(V_list, "groups") <- grps
      names(V_list) <- levels(grps)

    } else {

      # if there are groups then make block-diagonal matrix-lists

      if (is.null(mod$modelStruct$varStruct)) {
        V_list <- tapply(rep(sigma^2, length(mod$groups[[1]])),  mod$groups[[1]], diag)
      } else {
        sort_order <- get_sort_order(mod)
        sd_vec <- sigma / as.numeric(nlme::varWeights(mod$modelStruct$varStruct))[sort_order]
        V_list <- tapply(sd_vec^2, mod$groups[[1]], diag)
      }
      attr(V_list, "groups") <- mod$groups[[1]]
    }

  } else {

    # if there is a correlation structure,
    # build block-diagonals according to its grouping structure

    if (is.null(mod$modelStruct$varStruct)) {
      V_list <- if (sigma_scale) lapply(R_list, function(x) x * mod$sigma^2) else R_list
    } else {
      sort_order <- get_sort_order(mod)
      sd_vec <- sigma / as.numeric(nlme::varWeights(mod$modelStruct$varStruct))[sort_order]
      sd_list <- split(sd_vec, attr(R_list, "groups"))
      V_list <- Map(function(R, s) tcrossprod(s) * R, R = R_list, s = sd_list)
    }

    attr(V_list, "groups") <- attr(R_list, "groups")
  }

  return(V_list)
}

# Create block-diagonal covariance structure from Z-design and Tau matrices

ZDZt <- function(D, Z_list) {
  lapply(Z_list, function(z) z %*% D %*% t(z))
}

# Construct list of block-diagonal matrices for each random effects grouping structure

build_RE_mats <- function(mod, sigma_scale = FALSE) {

  # Get random effects structure
  all_groups <- rev(mod$groups)

  if (length(all_groups) == 1) {

    D_mat <- as.matrix(mod$modelStruct$reStruct[[1]])
    if (sigma_scale) D_mat <- mod$sigma^2 * D_mat
    Z_mat <- model.matrix(mod$modelStruct$reStruc, nlme::getData(mod))
    row.names(Z_mat) <- NULL
    Z_list <- matrix_list(Z_mat, all_groups[[1]], "row")
    ZDZ_list <- ZDZt(D_mat, Z_list)

    attr(ZDZ_list, "groups") <- all_groups[[1]]

  } else {
    if (sigma_scale) {
      D_list <- lapply(mod$modelStruct$reStruct, function(x) mod$sigma^2 * as.matrix(x))
    } else {
      D_list <- lapply(mod$modelStruct$reStruct, as.matrix)
    }
    Z_mat <- model.matrix(mod$modelStruct$reStruc, nlme::getData(mod))
    Z_names <- sapply(strsplit(colnames(Z_mat), ".", fixed=TRUE), function(x) x[1])
    row.names(Z_mat) <- NULL
    Z_levels <- lapply(names(all_groups), function(x) Z_mat[,x==Z_names,drop=FALSE])
    Z_levels <- Map(matrix_list, x = Z_levels, fac = all_groups, dim = "row")

    ZDZ_lists <- Map(ZDZt, D = D_list, Z_list = Z_levels)
    # ZDZ_lists <- Map(function(x,fac) x[order(fac)], x = ZDZ_lists, fac = all_groups)

    for (i in 2:length(all_groups)) {
      ZDZ_lists[[i]] <- add_bdiag(small_mats = ZDZ_lists[[i-1]],
                                  big_mats = ZDZ_lists[[i]],
                                  crosswalk = all_groups[c(i-1,i)])
    }

    ZDZ_list <- ZDZ_lists[[i]]

    attr(ZDZ_list, "groups") <- all_groups[[i]]

  }

  ZDZ_list

}

build_Sigma_mats <- function(mod, invert = FALSE, sigma_scale = FALSE) UseMethod("build_Sigma_mats")

build_Sigma_mats.default <- function(mod, invert = FALSE, sigma_scale = FALSE) {
  mod_class <- paste(class(mod), collapse = "-")
  stop(paste0("Sigma matrices not available for models of class ", mod_class, "."))
}

build_Sigma_mats.gls <- function(mod, invert = FALSE, sigma_scale = FALSE) {

  # lowest-level covariance structure
  V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)
  V_grps <- attr(V_list, "groups")

  if (invert) {
    V_list <- lapply(V_list, function(x) chol2inv(chol(x)))
  }

  attr(V_list, "groups") <- V_grps

  return(V_list)
}

build_Sigma_mats.lme <- function(mod, invert = FALSE, sigma_scale = FALSE) {

  if (inherits(mod, "nlme")) stop("not implemented for \"nlme\" objects")

  # lowest-level covariance structure
  V_list <- build_var_cor_mats(mod, sigma_scale = sigma_scale)

  # random effects covariance structure
  ZDZ_list <- build_RE_mats(mod, sigma_scale = sigma_scale)

  V_grps <- attr(V_list, "groups")

  # Check if lowest-level covariance structure is nested within RE structure
  ZDZ_grps <- attr(ZDZ_list, "groups")
  group_mapping <- tapply(ZDZ_grps, V_grps, function(x) length(unique(x)))
  nested <- all(group_mapping == 1L)

  if (nested) {
    Sigma_list <- add_bdiag(V_list, ZDZ_list, data.frame(V_grps, ZDZ_grps))
    Sigma_grps <- attr(ZDZ_list, "groups")
  } else {
    V_mat <- unblock(V_list, block = V_grps)
    ZDZ_mat <- unblock(ZDZ_list, block = ZDZ_grps)
    Sigma_list <- V_mat + ZDZ_mat
    Sigma_grps <- factor(rep("A", nrow(Sigma_list)))
  }

  if (invert) {
    Sigma_list <- lapply(Sigma_list, function(x) chol2inv(chol(x)))
  }

  attr(Sigma_list, "groups") <- Sigma_grps

  return(Sigma_list)
}

