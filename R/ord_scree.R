ord_scree = function(R, n=10, as_percent=TRUE, main=NULL, as_grob=FALSE) {
  # Plot variances of ordination axes.
  # n: Number of axes to plot (excluding partial and total, if present and requested)
  # as_percent: logical or 'group'. 
  #   - If FALSE, gives raw variance values and also the total variance.
  #   - If TRUE, gives percent of total variance.
  #   - If 'group', gives percent of variance within the group, and excludes partial and unconstrained axes.
  # as_grob: as a grid object (via ggplot2). Otherwise uses base graphics.
  
  vars = ord_variance(R)
  ylab = vars$type
  ylab = paste0(toupper(substr(ylab, 1, 1)),
                substr(ylab, 2, 9999))
  
  if (as_percent=='group') {
    ylab = paste('Percent of Within-Group', ylab)
  } else if (as_percent) {
    ylab = paste('Percent of Total', ylab)
  }
  
  p = ifelse(is.null(R$partial), 0, 1)
  
  if (as_percent=='group') {
    vars = vars$constrained
    n = min(n, ncol(vars))
    vars = vars['pct_group', 1:n]
    
  } else if (as_percent) {
    vars = with(vars, 
                cbind(partialed, 
                      constrained, 
                      unconstrained)
    )
    
    n = min(n + p, ncol(vars))
    rest = min(n+1, ncol(vars))
    vars = cbind(vars[, 1:n],
                 Remaining = rowSums(vars[, rest:ncol(vars)])
    )
    
    vars = vars['pct_total', ]
    
  } else {
    vars = with(vars, 
                cbind(total,
                      partialed, 
                      constrained, 
                      unconstrained)
    )
    
    n = min(n + p + 1, ncol(vars))
    vars = vars['value', 1:n]
    
  }
  
  colors = c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E')
  color_select = as.factor(substr(names(vars), 1, 2))
  
  if (as_grob) {
    if (!require(ggplot2)) {
      install.packages('ggplot2')
      require(ggplot2)
    }
    
    pltdf = data.frame(
      Y = vars,
      X = factor(names(vars), levels=names(vars)),
      grp = color_select
    )
    
    ggplot(pltdf, 
           aes(x=X,
               y=Y,
               fill=grp)) +
      scale_fill_brewer(palette='Dark2',
                        guide=FALSE) +
      geom_col() +
      ylab(ylab) +
      xlab(NULL) +
      ggtitle(main) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=90))
    
  } else {
    barplot(vars, ylab=ylab, horiz=FALSE, col=colors[color_select], border=NA, main=main)
  }
}