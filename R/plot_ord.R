# Create a biplot of an ordination (uses base graphics)
# Automatically formats the plot into groups, which the user has some control over.

plot_ord = function(R,
                    axes=c(1,2), 
                    row_group=NULL, 
                    max_labels=15, 
                    main=NULL,
                    continuous_scale=NULL,
                    discrete_scale=NULL,
                    shape_scale=NULL,
                    col_text='tomato4',
                    col_arrows='black',
                    shape_legend_position=NULL,
                    color_legend_position=NULL,
                    strip_column_label=NULL) {
  # Generate a biplot (or triplot!) of a pcwOrd ordination. plot_ord automatically
  # detects constrains and plots centeroids or biplot arrows as appropriate. You can
  # also add groupings to rows. For finer control, use ord_scores to calculate scores.
  #
  # R: pcwOrd ordination object
  # axes: which axes do you want to plot?
  # row_group: Supplemental grouping or gradient for rows. Rownames must match R$row_names
  # max_labels: Maximum number of column labels to write (the ones furthest from the origin are chosen). 
  #   Other labels are replaced as periods.
  # main: plot title
  # continuous_scale, discrete_scale: Color scales for groupings/gradients of row points
  # shape_scale: Customize shapes for rows - for categorical groups
  # col_text, col_arrows: color for column labels and biplot arrows and labels.
  # shape_legend_position: absolute coordinates of the shape legend, if present. Defaults to the top-right.
  # color_legend_position: position of the color legend, if present. 
  #    Values are relative coordinates (i.e. of ordination bounding box) defining a shift from the position of
  #    the shape legend (or top-right).
  # strip_column_label: Remove a string from the column text. See ord_scores(strip_label)
  
  #### Generate scores ####
  # Defaults
  add_X = TRUE
  add_Z = TRUE
  
  # Identify grouping variables
  nX = ncol(data.frame(R$X))
  nZ = ncol(data.frame(R$Z))
  nG = ncol(data.frame(row_group))
  # if (is.null(row_group)) row_group = FALSE # for downstream compatibility
  
  if (sum(nX, nZ, nG) > 2) {
    add_X = FALSE
    add_Z = FALSE
    if (nX > 0) X = as.matrix(R$X)[, 1, drop=FALSE]
    if (nZ > 0) Z = as.matrix(R$Z)[, 1, drop=FALSE]
    if (nG > 0 & !is.null(row_group)) G = as.matrix(row_group)[R$row_names, , drop=FALSE]
    if (nG >= 2) {
      row_group = G  # manually specified takes precedence
    } else if (nG == 1) {
      first_group = data.frame(X, Z)[, 1, drop=FALSE]
      row_group = data.frame(first_group, G)
    } else {
      row_group = data.frame(X, Z)
    }
    rownames(row_group) = R$row_names
  }
  
  if (is.null(row_group)) row_group=FALSE
  row = ord_scores(R, 'row', 'principle', axes=axes, add_X=add_X, add_Z=add_Z, add_grouping=row_group)
  col = ord_scores(R, 'column', 'standard', axes=axes, add_label=max_labels, strip_label=strip_column_label)
  col = scale_scores(col, row)
  
  bip = NULL
  grp = NULL
  # is_triplot = FALSE
  if (is.factor(R$X[[1]])) {
    grp = ord_scores(R, 'centeroid', 'principle', axes=axes, add_label=TRUE)
    # is_triplot = TRUE
  }
  if (!is.null(R$X)){
    if (is.numeric(as.matrix(R$X[1]))) {
      bip = ord_scores(R, 'biplot', 'standard', axes=axes, add_X=TRUE, add_label=TRUE)
      attr(bip, 'nAxes') = min(length(axes), sum(sapply(R$X, is.numeric)))
      bip = scale_scores(bip, row)
      # is_triplot = TRUE
    }
  }
  
  #### Format colors and shapes ####
  # set scales
  if (is.null(continuous_scale)) {
    continuous_scale = colorRampPalette(colors=c('steelblue4', 
                                                 'darkgoldenrod2', 
                                                 'tomato4'))
  }
  
  if (is.null(discrete_scale)) {  # this is c('Dark2', 'Set3') from RColorBrewer
    discrete_scale = c('#1B9E77','#D95F02','#7570B3','#E7298A',
                       '#66A61E','#E6AB02','#A6761D','#666666',
                       '#8DD3C7','#FFFFB3','#BEBADA','#FB8072',
                       '#80B1D3','#FDB462','#B3DE69','#FCCDE5',
                       '#D9D9D9','#BC80BD','#CCEBC5','#FFED6F')
  }
  
  if (is.null(shape_scale)) {
    shape_scale = c(19, 17, 15, 18, 0, 2, 1, 5, 8:14)
  }
  
  # Set defaults to modify if data requires
  continuous_legend = FALSE
  discrete_legend = FALSE
  colors = discrete_scale[1]
  
  shape_legend = FALSE
  color_legend = FALSE
  shapes = 19
  
  # Select columns (of rowscores) to use for formatting
  # Main difference is continuous vs categorical (factors)
  row_groups = row[, -c(1:2), drop=FALSE] # only two scales currently possible
  is_factor = sapply(row_groups, is.factor)
  is_continuous = sapply(row_groups, is.numeric)
  
  # Continuous colors and categorical shapes
  if (any(is_continuous)) { # continuous takes precedence for colors
    # dinking around to make a continuous colorscale
    row_colors = row_groups[is_continuous][,1]
    breakpoints = findInterval(row_colors,   # via stack exchange
                               sort(row_colors))
    colors = continuous_scale(length(breakpoints))[breakpoints]
    continuous_legend = TRUE
    
    # Shapes reserved for the factors
    if (any(is_factor)) {
      row_shape = row_groups[is_factor][[1]]
      shapes = shape_scale[row_shape]
      shape_legend = TRUE
      col_group = FALSE
    }
    
    # Categorical colors and shapes
  } else if (any(is_factor)) {
    # Categorical colors
    row_colors = row_groups[is_factor][[1]]
    colors = discrete_scale[row_colors]
    
    col_group = TRUE
    
    if (is.null(grp)) {
      color_legend = TRUE
    } else if (any(sapply(levels(grp$LABEL), function(x) !(x %in% levels(row_colors))))) {
      color_legend = TRUE
      col_group = FALSE
    }
    
    if (sum(is_factor) > 1) {
      row_shape = row_groups[is_factor][[2]]
      shapes = shape_scale[row_shape]
      shape_legend=TRUE
    }
  }
  
  # Axis lables 
  variances = ord_variance(R)$unconstrained['pct_total', ]
  # rownames(variances)[2] = 'var'
  if (!is.null(R$constrained)) {
    vv = as.matrix(ord_variance(R)$constrained)
    if (ncol(vv) == 1) {
      vv = vv['pct_total', ]
    } else {
      vv = vv['pct_group', ]
    }
    # rownames(vv)[2] = 'var'
    variances = c(vv, variances)
  }
  variances = variances[axes]
  
  axis_labels = paste(colnames(row), ' [', signif(variances, 3), '%]', sep="")
  
  # Plot....
  all_scores = rbind(row[, 1:2],
                     col[, 1:2],
                     bip[, 1:2],
                     grp[, 1:2])
  plot(all_scores, 
       type='n', asp=1, 
       xlab=axis_labels[1], ylab=axis_labels[2], main=main)
  grid(col='gray80', lty=1)
  axis(1)
  axis(2)
  # axes
  abline(h=0, col='gray50')
  abline(v=0, col='gray50')
  
  text(col[, 1:2], 
       labels=col$LABEL, 
       col=col_text, 
       cex=0.6)
  
  points(row[, 1:2], 
         pch=shapes, 
         col=adjustcolor(colors, alpha.f=0.7), 
         cex=1)
  
  if (!is.null(grp)) {
    if (col_group) {
      col_group = discrete_scale[grp$LABEL]
    } else {
      col_group = '#333333'
    }
    text(grp[, 1:2], 
         label=grp$LABEL, 
         col=adjustcolor(col_group, offset=c(-0.2, -0.2, -0.2, 0)), # darken slightly for easier reading
         font=2)
  }
  
  if (!is.null(bip)) {
    arrows(x0=0, y0=0,
           x1=0.9*bip[, 1],
           y1=0.8*bip[, 2],
           length=0.1,
           col=col_arrows)
    text(bip[, 1:2],
         label=bip$LABEL,
         col=adjustcolor(col_arrows, offset=c(-0.2, -0.2, -0.2, 0)), # darken slightly for easier reading
         cex=0.7)
  }
  
  if (any(shape_legend, continuous_legend, discrete_legend, color_legend)) {
    bounds = c(ymin = min(all_scores[, 2]), # for positioning first legend
               ymax = max(all_scores[, 2]),
               xmin = min(all_scores[, 1]),
               xmax = max(all_scores[, 1]))
    # for positioning second legend relative to first
    max_range = max(c(y = bounds['ymax'] - bounds['ymin'],  # max because aspect=1
                      x = bounds['xmax'] - bounds['xmin']))
    
    base_position = c(bounds['xmin'], bounds['ymax'])
    sl_columns = 0 # shape legend columns. For auto-positioning colors
    sl_rows = 0
    
    if (shape_legend) {
      slp = shape_legend_position
      if (is.null(slp)) slp = base_position
      
      base_position = slp  # update for color legend if needed
      posX = slp[1]
      posY = slp[min(2, length(slp))]
      
      lvls = levels(row_shape)
      nlvls = length(lvls)
      sl_columns = ceiling(nlvls/6)
      sl_rows = ceiling(nlvls/sl_columns)
      
      legend(posX, posY,
             legend=lvls,
             col='black',
             pch=shape_scale[1:nlvls],
             cex=0.7,
             bty='n',
             ncol=sl_columns)
    }
    
    if (color_legend | continuous_legend) {
      clp = color_legend_position
      if (is.null(clp)) clp = c(0, -0.1)
      
      x_offset = max_range * (clp[1] + (clp[1]!=0)*0.1*sl_columns)
      y_offset = max_range * (clp[2] - (clp[2]!=0)*0.025*sl_rows)
      
      posX = base_position[1] + x_offset
      posY = base_position[2] + y_offset
      
      if (is.factor(row_colors)) {
        lvls = levels(row_colors)
        colors = discrete_scale
      } else {
        lvls = round(
          c(min(row_colors, na.rm=TRUE),
            mean(row_colors, na.rm=TRUE),
            max(row_colors, na.rm=TRUE)),
          3)
        colors = continuous_scale(3)
      }
      nlvls = length(lvls)
      nc = ceiling(nlvls/6)
      legend(posX, posY,
             legend=lvls,
             col=colors[1:nlvls],
             pch=19,
             cex=0.7,
             bty='n',
             ncol=nc)
    }
  }
  
}