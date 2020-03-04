#### Overview ####
# Partial (weighted) RDA functions
# 
# Motivation: easyCODA performs weighted RDA, but cannot partial out additional effects. 
# Pre-partialing them doesn't play well with either easyCODA or vegan. vegan doesn't 
# perform weighted RDA. Also, easyCODA is s-l-o-w. vegan's code is much faster.
#
# Algorithm:
#   - X is weighted by rows and columns at the beginning
#   - Y and Z are weighted by rows and columns as called
#   - The fitted values give centeroid scores (CCA)
#   - The residuals are regressed or decomposed in the next step.
#
# A number of steps are repeated, and so are packaged into hidden functions:
#   1. .make_weights(): compute row- and column-weights for a matrix
#   2. .make_dummy(): make a dummy matrix without an intercept
#   3. .prepare_matrix(): Center and weight a matrix
#   4. .regress(): Perform linear regression via qr decomposition
#   5. .remove_zeros(): Remove ordination axes that are (machine) zero.
# Most of these hidden functions are also called by permute_ord().

.make_weights = function(P, choice, weight=FALSE, as_vegan=FALSE) {
  # Make vectors of matrix weights for either rows or columns: One vector for the response data
  # (wt_init), and another for the predictors. 
  # 
  # as_vegan refers to vegan-style row- and column-weights applied to the response matrix
  # for unweighted ordination, which are 1/(nrow-1) for rows and 1 for columns
  # Otherwise, both rows and columns have easyCODA-style weights of 1/nrow or 1/ncol, respectively.
  # 
  # if weight==TRUE, then weights by marginals (1/rowSums or 1/colSums)
  # 
  # if weight is specified, then those are used
  
  choice = substr(tolower(choice), 1, 1)
  
  n = switch(choice,
             'r'=nrow(P),
             'c'=ncol(P))
  
  if (!weight[1]) {
    wt = rep(1/n, n)
    
    if (as_vegan) {          # vegan style
      wt = switch(choice,
                       'r' = rep(1/(n-1), n), 
                       'c' = rep(1, n))
      
    }
  } else if (is.logical(weight)) {
    wt = switch(choice,
                'r'=1/rowSums(P),
                'c'=1/colSums(P))
  } else {
    wt = weight
  }
  return(wt)
}

.make_dummy = function(X=NULL, nr=NULL) {
  # turn X into a dummy variable without an intercept. Or return rep(1, nr) as a matrix.
  if (is.null(X)) X = matrix(rep(1, nr))
  
  X = as.data.frame(X)
  X = model.matrix(~., data=X)[, -1, drop=FALSE]
  
  return(X)
}

.prepare_matrix = function(X, rw=NA, cw=NA, as_CA=FALSE) {
  # center and weight a matrix
  
  nr = nrow(X)
  nc = ncol(X)
  
  if (is.na(rw[1])) rw = rep(1, nr)
  if (is.na(cw[1])) cw = rep(1, nc)
  
  if (as_CA) {
    X = X/sum(X)
  } else {
    colmeans = apply(X, 2, weighted.mean, w=rw)
    X = scale(X, center=colmeans, scale=FALSE) # weighted column means?
  }
  
  X = sweep(X, 1, sqrt(rw), '*') # scale rows
  X = sweep(X, 2, sqrt(cw), '*') # scale columns
  
  return(X)
}

.regress = function(Y, X) {
  # via QR decomposition
  # Y is response, X is predictor
  
  Q = qr(X)
  out = list(
    fitted = qr.fitted(Q, Y),
    resids = qr.resid(Q, Y)
  )
  return(out)
}

.remove_zeros = function(S) {
  # remove columns/eigenvectors from SVD that are (machine) zero
  
  not_zero = S$d > sqrt(.Machine$double.eps)  # from vegan
  
  S$d = S$d[not_zero]
  S$u = S$u[, not_zero, drop=FALSE]
  S$v = S$v[, not_zero, drop=FALSE]
  
  return(S)
}

pcwOrd = function(Y, 
                  X=NULL, 
                  Z=NULL, 
                  weight_rows=FALSE, 
                  weight_columns=FALSE, 
                  as_CA=FALSE, 
                  as_vegan=FALSE) {
  # Conduct a (partial) (constrained) (weighted) ordination. Currently only runs principle component analysis
  # and derivatives. May include correspondence analaysis in the near future.
  #
  # Y: Community matrix. Columns are features, rows are samples. Alternatively, an object from easyCODA::CLR.
  # X: Constraining data.frame. Treatments or environmental variables.
  # Z: Partialing data.frame. Treatments or environmental variables.
  # weight_rows: logical or numeric. If TRUE, will weight by rowSums(Y).
  # weight_columns: logical or nueric. If TRUE, will weight by colSums(Y).
  # as_CA: NOT IMPLEMENTED. run as correspondence analysis? If TRUE, weights of rows and columns are overrridden as TRUE
  #   and Y is double-standardized. If FALSE, Y is centered and weighted as indicated.
  # as_vegan: weight the data matrix Y as in vegan?
  
  method = 'RDA'
  if (as_CA) {
    stop('as_CA is not yet implemented. You might try easyCODA::CCA or vegan::cca.')
    weight_rows = TRUE
    weight_columns = TRUE
    method = 'CCA'
  }
  
  # for easyCODA compatibility
  if (hasName(Y, 'LR.wt') & hasName(Y, 'LR')) {
    weight_columns = Y$LR.wt
    Y = Y$LR
  }
  
  nr = nrow(Y)
  nc = ncol(Y)
  
  # Begin preparing Y
  Y_orig = Y
  Y = as.matrix(Y)
  
  # Define weights
  rw = .make_weights(Y, 'row', weight_rows, as_vegan)
  cw = .make_weights(Y, 'col', weight_columns, as_vegan)
  
  # scale and weight Y
  Y = .prepare_matrix(Y, rw=rw, cw=cw, as_CA)
  
  # Partial Z from Y
  if (is.null(Z)) {
    pY = Y
    partial = NULL
    X_dummy = X
  } else {
    Z_dummy = .make_dummy(Z, nr)
    Z_act = .prepare_matrix(Z_dummy, rw)
    out = .regress(Y, Z_act)
    
    pY = out$resids
    partial = out$fitted
    if (!is.null(X)) {
      X_dummy = cbind(Z, X)
    }
  }
  
  # Constrain pY by X
  if (is.null(X)) {
    rY = pY
    fitted = NULL
    constrained = NULL
    X_weights = NULL
  } else {
    X_dummy = .make_dummy(X_dummy, nr)
    # X_dummy = .make_dummy(X, nr)
    X_act = .prepare_matrix(X_dummy, rw)
    out = .regress(pY, X_act)
    
    rY = out$resids
    fitted = out$fitted
    
    # ordinate
    constrained = .remove_zeros(svd(fitted))
    
    # Calculate group weights
    zcol = ifelse(is.null(Z), 0, ncol(Z_dummy))
    X_weights = 1/sqrt(colSums(X_act^2))[-c(1:zcol)]
  }
  
  # Ordinate residuals
  unconstrained = .remove_zeros(svd(rY))
  
  # To return:
  out = list(
    method = method,
    Y = Y_orig,
    Y_scaled = Y,
    Y_partial = pY,
    Y_residual = rY,
    X = X,
    X_weights = X_weights,
    Z = Z,
    row_weights = rw,
    row_names = rownames(Y_orig),
    col_weights = cw,
    col_names = colnames(Y_orig),
    group_names = colnames(X),
    partial = partial,
    fitted = fitted,
    constrained = constrained,
    unconstrained = unconstrained
  )
  attr(out, 'class') = 'pcwOrd'
  attr(out, 'method') = ifelse(as_CA, 'CA', 'PCA')
  
  return(out)
}
