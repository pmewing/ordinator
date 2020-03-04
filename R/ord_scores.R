#### Scores Functions ####
# ord_scores: Calculate scores for rows, columns, etc from a pcwOrd object. Returns a data.frame that can include environmental or other metadata for visualization.
# scale_scores: Scale score1 to match the range of score2, or as defined by scaling. For easier visualization.
# top_scores: Return the top n scores from ord_scores, as defined by distance from the origin.



ord_scores = function(R, 
                      choice=c('row', 'column', 'centeroid', 'biplot'), 
                      scaling=c('standard', 'principle', 'contribution'), 
                      constrained=TRUE, axes=c(1,2), add_X=FALSE, add_Z=FALSE, 
                      add_grouping=FALSE, 
                      add_label=FALSE,
                      strip_label=NULL,...) {
  # Calculate scores of an ordination.
  #
  # R: A pcwOrd object
  # choice: which scores you want to return. centeroid and biplot are for discrete and continuous constraints, respectively.
  # scaling: Which score to calculate? Generally, rows are principle, columns are standard, and if we want important columns, we look at contributions.
  #   contributions are per Greenacre, 2010. Biplots in practice.
  # constrained: return constrained axes (if they exist)?
  # axes: integer. Axes for which to return scores.
  # add_X, add_Z: Include values from X (constraints) or Z (partials) in the row scores data.frame?
  # add_grouping: Add values from a grouping to the data.frame? Rownames must match R$row_names.
  # add_label: Include a column, "LABEL", with the row/column/group names of each score?
  # strip_label: character to remove from add_label. For example ("OTU_"). 
  
  stopifnot(class(R) == 'pcwOrd')
  
  choice = match.arg(choice)
  scaling = match.arg(scaling)
  
  is_unconstrained = is.null(R$constrained)
  
  type = ifelse(constrained & !is_unconstrained,
                'constrained',
                'unconstrained')
  
  if (is_unconstrained & (choice %in% c('centeroid', 'biplot'))) {
    msg = paste('cannot plot', choice, 'as this ordination lacks constraints')
    stop(msg)
  }
  
  # need:
  # weights: by choice
  # singular values: by type
  # left- or right-sides: by choice
  # (maybe) Y - for inferring constrained row scores
  # If constrained, pull all unconstrained axes, as well.
  # If these are centeroids, need to calculate the centeroid on the unconstrained axis.
 
  if (choice == 'centeroid') {
    X_grp = as.data.frame(R$X)
    if (any(sapply(X_grp, is.numeric))) {
      stop('not sure how to calculate centeroids when mixed with continuous variables')
    }
    
    X_grp = apply(X_grp, 1, paste, collapse='_')
  }
  
  if (choice == 'biplot') {
    X_bp = .make_dummy(R$X)
    bp_names = colnames(X_bp)
    # bp_names = sapply(bp_names, function(x) {
    #   for (i in colnames(R$X)) {
    #     x = gsub(i, '', x)
    #   }
    #   return(x)
    # })
    X_bp = .prepare_matrix(X_bp, 
                           rw=R$row_weights)
    X_bp = as.matrix(X_bp)
  }
  
  wt = switch(choice,
              'row' = R$row_weights,
              'column' = R$col_weights,
              'centeroid' = tapply(R$row_weights, X_grp, mean),  # aggregate row_weights by group
              'biplot' = rep(1, ncol(X_bp))
  )
  
  
  
  V = R[[type]]$v
  U = R[[type]]$u
  d = R[[type]]$d
  if (length(R$X_weights)==0) R$X_weights = 1
  
  scores = switch(choice,
                  'row' = U,
                  'column' = V,
                  'centeroid' = apply(as.matrix(U), 2, tapply, X_grp, mean), # aggregate scores by group
                  'biplot' = sweep(as.matrix(crossprod(X_bp, U)), 1, R$X_weights, '*')
  )
  
  # Calculate unconstrained if necessary
  if (type == 'constrained') {
    uV = R$unconstrained$v
    uU = R$unconstrained$u
    ud = R$unconstrained$d
    
    ss = switch(choice, 
                'row' = uU,
                'column' = uV,
                'centeroid' = matrix(rep(0, length(ud)*nrow(scores)),
                                     nrow=nrow(scores)), 
                'biplot' = matrix(rep(0, length(ud)*nrow(scores)), 
                                  nrow=nrow(scores)))
    d = c(d, ud)
    U = cbind(U, uU)
    V = cbind(V, uV)
    scores = cbind(scores, ss)
  }
      
  if (choice=='row' & type=='constrained') {
    # infer point scores based on the column side of svd
    scores = R$Y_partial %*% V
    scores = sweep(scores, 2, d, '/')
  }
  
  out = switch(scaling,
               'contribution' = scores,
               'standard' = sweep(scores, 1, sqrt(wt), '/'),
               'principle' = 
                 sweep(
                   sweep(scores, 1, sqrt(wt), '/'), # scale standard score by SVs
                   2, d, '*')
  )
  
  if (is.null(axes) | is.character(axes)) {
    axes = 1:ncol(out)
  }
  
  # Column/axis names
  n_axes = ncol(out)
  axis_names = rep('Axis', n_axes)
  
  if (type == 'unconstrained') {
    if (!is_unconstrained) {
      axis_names = paste0(rep('u', n_axes), 
                          axis_names,
                          1:n_axes)
    } else {
      axis_names = paste0(axis_names, 1:n_axes)
    }
  }
  if (type == 'constrained') {
    axis_names = 
      paste0(
        c(rep('c', length(R$constrained$d)),
          rep('u', length(R$unconstrained$d))), 
        axis_names,
        c(1:length(R$constrained$d),
          1:length(R$unconstrained$d)))
  }
  if (!is.null(R$partial)) {
    axis_names = paste0('p',
                        axis_names)
  }
  colnames(out) = axis_names
  
  # Row names
  rownames(out) = switch(choice,
                         'row' = R$row_names,
                         'column' = R$col_names,
                         'centeroid' = rownames(out),
                         'biplot' = bp_names
  )
  
  # Select Axes
  out = out[, axes, drop=FALSE]
  out = as.data.frame(out)
  
  # Additional data
  if (add_X & !is.null(R$X)) {
    if (choice == 'column') stop("Cannot add X to column scores")
    if (choice == 'row') {
      grp = as.matrix(R$X)
      rownames(grp) = rownames(out)
      out = data.frame(out, grp)
    } else {
      out = data.frame(
        out,
        GROUP = rownames(out)
      )
    }
  }
  
  if (add_Z & !is.null(R$Z)) {
    if (choice != 'row') stop("Can only add Z to row scores")
    grp = as.matrix(R$Z)
    rownames(grp) = rownames(out)
    out = data.frame(out, grp)
  }
  
  add_grouping = as.data.frame(add_grouping)
  if (add_grouping[1,1] != FALSE | prod(dim(grouping)) > 1) {
    grp = add_grouping[rownames(out), , drop=FALSE]
    out = data.frame(out, grp)
  }
  
  if (add_label==TRUE) add_label = 15
  if (is.numeric(add_label)) {
    constr = type=='constrained'
    top_labels = top_scores(R, n=add_label, choice=choice,
                             scaling=scaling,
                             axes=axes, constrained=constr)
    
    ll = sapply(rownames(out), function(x) {
      ifelse(x %in% rownames(top_labels), x, '.')
    })
    if (!is.null(strip_label)) {
      ll = gsub(strip_label, '', ll)
    }
    ll_levels = ll[!duplicated(ll)] # preserve order
    out$LABEL = factor(ll, levels=ll_levels)
  }
  
  attr(out, 'class') = c('data.frame', 'pcwOrd_scores')
  attr(out, 'scaling') = scaling
  attr(out, 'type') = type
  attr(out, 'nAxes') = length(axes)
  return(out)
}

scale_scores = function(score1, score2=NULL, scaling=NULL) {
  # rescale score1 scores. Two methods:
  # 1. Auto - if score2 is given, then score1 will be multiplied by a constant so that it has the same root-mean-squared values as score2
  # 2. Defined - if scaling is given as a scaler. Overriden by score2 if given.
  
  stopifnot(any(class(score1) == 'pcwOrd_scores'))
  
  n_axes = attr(score1, 'nAxes')
  
  if (!is.null(score2)) {
    stopifnot(any(class(score2) == 'pcwOrd_scores'))
    
    n2 = attr(score2, 'nAxes')
    n_axes = min(n_axes, n2)
    
    if (attr(score1, 'scaling') == attr(score2, 'scaling')) {
      msg = paste('score1 and score2 both are', attr(score1, 'scaling'), 'coordinates.')
      warning(msg)
    }
    # calculating the scaling factor
    
    ranges = list(score1, score2)
    ranges = sapply(ranges, function(x) {
      x = x[, c(1:n_axes), drop=FALSE]
      x = rbind(rep(0, n_axes),
                x)
      apply(x, 2, function(y) diff(range(y)))
    }) 
    if (n_axes == 1) ranges = matrix(ranges, nrow=1)
    scaling = mean(ranges[, 2] / ranges[, 1])
    message(paste("Rescaling", attr(score1, 'scaling'), "scores by", signif(scaling, 4)))
  }
  

  score1[, c(1:n_axes)] = score1[, c(1:n_axes), drop=FALSE] * scaling
  return(score1)
}

top_scores = function(R, n=15, choice=c('column', 'row', 'centeroid', 'biplot'),
                      scaling=c('contribution', 'standard', 'principle'),
                      axes='all', constrained=TRUE) {
  # Return the top n scores from the ordination. Scores are calculated using ord_scores,
  # then (euclidian) distance from the origin is calculated. See ord_scores for argument details.
  # R is a pcwOrd ordination object.
  
  
  scores = ord_scores(R, 
                      choice=match.arg(choice),
                      scaling=match.arg(scaling),
                      constrained=constrained,
                      axes=axes)
                      
  if (axes[1] == 'all') axes = 1:ncol(scores)
  
  n = min(n, nrow(scores))
  
  dist = apply(scores, 1, function(x) sqrt(sum(x^2)))
  dist = sort(dist, decreasing=TRUE)
  dist = as.matrix(dist)[1:n, , drop=FALSE]
  colnames(dist) = paste0('mean_', match.arg(scaling))
  
  attr(dist, 'axes') = axes
  
  return(dist)
  
}
