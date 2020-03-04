permute_ord = function(R, permute_on=c('partial', 'initial', 'fitted'), times=999, return_permutations=FALSE) {
  # Permutational ANOVA of the pcwOrd results.
  # Arguments
  #   R: pcwOrd object
  #   permute_on: level at which to permute. 
  #   times: N permutations
  #   return_permutations: whether to return the permuted F values.
  # 
  # Notes:
  #   - initial = initial data. Equal to vegan::permutest(model='direct').
  #   - partial = residuals of partialing (default). Equal to vegan::permutest(model='reduced')
  #   - fitted = residuals of constraints. Equal to vegan::permutest(model='full')
  # 
  # Permutes the given matrix, runs through the partialing and conditioning steps, and calculates F.
  # Compares observed F to the permuted F. Returns p_val = p(Fobs < Fpermute).
  
  require(permute)
  
  stopifnot('pcwOrd' %in% class(R))
  
  permute_on = match.arg(permute_on)
  
  if (is.null(R$constrained)) stop("Can't test a model without a constraining variable!")
  
  iY = R$Y_scaled    # initial (scaled/weighted) data
  pY = R$Y_partial   # data after partialing Z
  fY = R$fitted      # data fitted to x
  rY = R$Y_residual  # residuals after fitting X
  X = R$X
  Z = R$Z
  rw_orig = R$row_weights
  nr = length(rw_orig)
  
  zcol = ifelse(is.null(Z),
                0,
                ncol(.make_dummy(Z)))
  xcol = ncol(.make_dummy(X))
  
  num_df = xcol
  denom_df = nr - num_df - zcol - 1
  residuals = sum(rY^2)
  fitted = sum(fY^2)
  partial_resids = sum(pY^2)
  total_var = sum(iY^2)
  
  F_stat = (fitted/num_df) / (residuals/denom_df)
  
  # Permute
  
  P = switch(match.arg(permute_on), 
             'partial' = pY,
             'initial' = iY,
             'fitted'  = rY 
  )
  
  permutes = shuffleSet(nr, times)
  
  F_perm = apply(permutes, 1, function(i) {
    
    P = P[i, ]
    rw = rw_orig[i]
    
    # Partial Z from P
    if (is.null(Z)) {
      pP = P
      X_dummy = X
    } else {
      Z_dummy = .make_dummy(Z)
      Z_act = .prepare_matrix(Z_dummy, rw)
      out = .regress(P, Z_act)
      
      pP = out$resids
      X_dummy = cbind(Z, X)
    }
    
    # Constrain pY by X
    X_dummy = .make_dummy(X_dummy)
    X_act = .prepare_matrix(X_dummy, rw)
    out = .regress(pP, X_act)
    
    rP = out$resids
    fP = out$fitted
    
    # Calculate terms
    perm_fitted = sum(fP^2)
    perm_resids = sum(rP^2)
    F_perm = (perm_fitted/num_df) / (perm_resids/denom_df)
    
    return(F_perm)
  })
  
  alpha = sum(F_stat >= F_perm)/(times+1)
  
  if (!return_permutations) {
    F_perm = times
  }
  
  out = list(
    permute_on = permute_on,
    total_variance = total_var,
    variance_after_partialing = partial_resids,
    fitted = fitted,
    residuals = residuals,
    num_df = num_df,
    denom_df = denom_df,
    F_stat = F_stat,
    p_val = 1-alpha,
    F_perm = F_perm
  )
  
  if (!return_permutations) out = data.frame(out)
  
  return(out)
}