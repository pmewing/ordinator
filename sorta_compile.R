#!/usr/bin/R
# Lazy load all functions into an .rds object

if (!require(here)) {
  install.packages('here')
  if (!require(here)) stop("Cannot load package `here`")
}

list_files = list.files(here('R'))

for (i in list_files) {
  if (grepl('.R', i)) source(here('R', i))
}

rm(i); rm(list_files)

save.image(here('pcwOrd_0.1.Rdata'))
