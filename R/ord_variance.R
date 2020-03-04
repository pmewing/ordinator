.axis_variances = function(R) {
  # Calculates absolute and relative variances for each axis in pcwOrd object R.
  # called by ord_variance and plot_ord
  
  stopifnot(class(R) == 'pcwOrd')
  
  total = c(Total = sum(R[['Y_scaled']]^2))
  
  uc = R[['unconstrained']][['d']]^2
  names(uc) = paste0('uAxis', 1:length(uc))
  
  cc = R[['constrained']]
  if (!is.null(cc)) {
    cc = cc[['d']]^2
    names(cc) = paste0('cAxis', 1:length(cc))
  }
  
  pp = R[['partial']]
  if (!is.null(pp)) {
    pp = c(Partial = sum(pp^2))
  }
  
  out = list(total = total, 
             partialed = pp,
             constrained = cc,
             unconstrained = uc
  )
  
  out = lapply(out, function(vv) {
    rbind(
      value = vv,
      pct_group = 100*vv/sum(vv),
      pct_total = 100*vv/total
    )
  })
  
  return(out)
}

ord_variance = function(R) {
  # Return a list with absolute and relative variances of ordination axes
  # R is the pcwOrd object.
  
  axis_var = .axis_variances(R)
  
  aggregate_var = sapply(axis_var, rowSums)[-2, ]
  
  out = list(type = switch(attr(R, 'method'),
                           'CA' = 'inertia',
                           'PCA' = 'variance'),
             summary = aggregate_var)
             
  out = unlist(
    list(out,
         axis_var),
    recursive=FALSE)
  
  return(out)
}