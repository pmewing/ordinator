library(vegan)
library(easyCODA)
library(here)
library(magrittr)
library(ggplot2)
library(microbenchmark)

dd = readRDS(here('Spider.RDS'))
load('pcwOrd_0.1.Rdata')

{ 
  spp = dd$abund %>% 
    as.matrix
  env = dd$x %>% 
    as.data.frame
  
  spp[spp==0] = 0.5
  spp = sweep(spp, 1, rowSums(spp), '/')
  spp = CLR(spp)  # weighted uses weighted geometric mean --> mean
  weights = spp$LR.wt
  spp = spp$LR
  
  cats = c('dry', 'damp', 'wet', 'soaking')
  cc = ceiling(env$soil.dry)
  env$soil.cat = factor(cats[cc],
                         levels=cats)
  
  nobs = nrow(spp)
  row_names = c(letters[], LETTERS[])[1:nobs]
  rownames(env) = row_names
  rownames(spp) = row_names
}

Y = spp
X = env[,'soil.cat', drop=FALSE]
Z = env['moss']

test_results = c()
zero_val = 1E-10

# Test weighted PCA vs easyCODA
test = PCA(Y, weight=weights)
pme = pcwOrd(Y, weight_columns=weights)

plot_ord(pme)


test = sum(round(pme$unconstrained$d, 10) != round(test$sv[test$sv>zero_val], 10))
test_results %<>% c(weighted_pca = test)
       
# Test unweighted PCA vs easyCODA
tt = PCA(Y, weight=FALSE)
pme = pcwOrd(Y)

test = tt$sv
test = sum(round(pme$unconstrained$d, 10) != round(test[test>zero_val], 10))
test_results %<>% c(unweighted_pca = test)

pme_scores = ord_scores(pme)
test_scores = tt$rowcoord[, 1:2]
test = sum(round(pme_scores, 10) != round(test_scores, 10))
test_results %<>% c(row_standard_uw = test)

pme_scores = ord_scores(pme, 'column')
test_scores = tt$colcoord[, 1:2]
test = sum(round(pme_scores, 10) != round(test_scores, 10))
test_results %<>% c(col_standard_uw = test)

pme_scores = ord_scores(pme, scaling='principle')
test_scores = tt$rowpcoord[, 1:2]
test = sum(round(pme_scores, 10) != round(test_scores, 10))
test_results %<>% c(row_principle_uw = test)


# Test weighted RDA vs easyCODA
test = RDA(Y, cov=DUMMY(as.matrix(X)), weight=weights)
pme = pcwOrd(Y, X, weight_columns=weights)

test = test$sv
test = sum(round(pme$constrained$d, 10) != round(test[test>zero_val], 10))
test_results %<>% c(weighted_rda = test)

# Test unweighted RDA vs easyCODA
test = RDA(Y, cov=DUMMY(as.matrix(X)), weight=FALSE)
pme = pcwOrd(Y, X)

test = test$sv
test = sum(round(pme$constrained$d, 10) != round(test[test>zero_val], 10))
test_results %<>% c(unweighted_rda = test)

# Test unweighted PCA vs vegan
test = rda(Y)
pme = pcwOrd(Y, as_vegan=TRUE)

test = sqrt(test$CA$eig)
test = sum(round(pme$unconstrained$d, 10) != round(test[test>zero_val], 10))
test_results %<>% c(unweighted_vegan_pca = test)

# Test unweighted RDA vs vegan
tord = rda(Y, X)
pme = pcwOrd(Y, X, as_vegan=TRUE)

test = sqrt(tord$CA$eig)
t1 = pme$unconstrained$d
test = sum(round(pme$unconstrained$d, 10) != round(test[test>zero_val], 10))
test_results %<>% c(unweighted_vegan_rda_unconst = test)

test = sqrt(tord$CCA$eig)
test = sum(round(pme$constrained$d, 10) != round(test[test>zero_val], 10))
test_results %<>% c(unweighted_vegan_rda_const = test)

# Test unweighted pcwOrd vs vegan
# The partialing step seems to add a bit of noise, thus the less stringent results. 
# Possibly due to not including Z when constraining Y by X?
# But the error is very small (around 0.05%)
tord = rda(Y, X, Z)
pme = pcwOrd(Y, X, Z, as_vegan=TRUE)

test = tord$pCCA$tot.chi
t1 = sum(pme$partial^2)
test = sum(round(t1, 10) != round(test, 10))
test_results %<>% c(uw_veg_pcwOrd_part = test)

test = sqrt(tord$CCA$eig)
t1 = pme$constrained$d
test = sum(abs(t1-test)>1E-10)
test_results %<>% c(uw_veg_pcwOrd_const = test)

test = sqrt(tord$CA$eig)
t1 = pme$unconstrained$d
test = sum(abs(t1 - test) > 1E-10)  # one extra axis 
test_results %<>% c(uw_veg_pcwOrd_unconst = test)

print(test_results)


#### Test additional functions ####

# Run weighted pcwOrd
pme = pcwOrd(Y, X, weight_columns=FALSE, as_vegan=TRUE)
ord_variance(pme)

# Calculate scores
ord_scores(pme, choice='centeroid', scaling='standard', type='constrained', axes=c(1,2))

rr = ord_scores(pme, 'row', 'standard', add_X=TRUE, add_Z=TRUE)#, add_grouping=env['REP'])
cc = ord_scores(pme, 'column', 'standard', add_label=20)
gg = ord_scores(pme, 'centeroid', 'standard', add_X=TRUE)
bp = ord_scores(pme, 'biplot', 'standard', add_X = TRUE)

plot(rbind(rr[, 1:2], cc[, 1:2], gg[, 1:2], bp[, 1:2]), type='n', asp=1)
text(cc, label='.', col='tomato4')
cols = hcl.colors(length(levels(X[[1]])), 'Dark3')
points(rr, col=cols[X[[1]]], pch=19)
txcol = as.factor(rownames(gg))
text(gg, label=rownames(gg), col=cols[as.numeric(txcol)])
arrows(x0=0, y0=0, x1=bp[,1], y1=bp[,2])

# Check scores
calc_centers = apply(rr[, 1:2], 2, tapply, X, mean)
calc_centers = calc_centers[rownames(gg),]
test_results = c(test_results, scores_OK = sum(round(calc_centers, 5) != round(gg[, 1:2], 5)))

#### Print test results ####
# want all to be zero
test_results

#### Benchmark ####
microbenchmark(
  test = PCA(Y, weight=weights),
  pme = pcwOrd(Y, weight_columns=weights),
  veg = rda(Y),
  times=30,
  unit='ms')

#### Scree plot ####
ord_scree(pme, as_grob=TRUE, as_percent=TRUE)

#### Permutation test ####

# fix this. Disagrees with VEGAN for permuted F-values

# vegan options
# direct - permutes community matrix (after centering, weighting, and/or standardizing)
# reduced - permutes residuals after partialing the prepared community matrix. This is the default.
# full - permutes residual values after constraining. 

# my options
#   - initial = initial data. 
#   - partial = residuals of partialing. 
#   - fitted = residuals of constraints. 

# So: 
# initial = direct
# partial = reduced
# fitted = fulll

pme = pcwOrd(Y, X, Z, weight_columns=FALSE, weight_rows=FALSE, as_vegan=TRUE)
veg = rda(Y, X, Z)
plot_ord(pme)

n = 999
set.seed(123)
tt=permute_ord(pme, 'fitted', times=n, return_permutations=F) %T>% print
set.seed(123)
vv=permutest(veg, permutations=n, model='full') %T>% print

set.seed(123)
tt=permute_ord(pme, 'partial', times=n, return_permutations=F) %T>% print
set.seed(123)
vv=permutest(veg, permutations=n, model='reduced') %T>% print

set.seed(123)
tt=permute_ord(pme, 'initial', times=n, return_permutations=F) %T>% print
set.seed(123)
vv=permutest(veg, permutations=n, model='direct') %T>% print

test_results
