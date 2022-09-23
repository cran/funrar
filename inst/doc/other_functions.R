## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(funrar)

## ----load-data----------------------------------------------------------------
data("aravo", package = "ade4")

# Extract the traits of all species of `aravo`
traits = aravo$traits

head(traits)

## ----global-di----------------------------------------------------------------
# Compute a Euclidean distance matrix (because all traits are quantitative)
dist_matrix = compute_dist_matrix(traits, metric = "euclidean")

# Compute global distinctiveness
global_di = distinctiveness_global(dist_matrix)

head(global_di)

## ----global-di-2--------------------------------------------------------------
aravo_di = distinctiveness_global(dist_matrix, di_name = "aravo_di")

head(aravo_di)

## ----di-dimensions------------------------------------------------------------
di_dim = distinctiveness_dimensions(as.matrix(aravo$spe), traits, metric = "euclidean")

str(di_dim, max.level = 1)

## ----ui-dimensions------------------------------------------------------------
ui_dim = uniqueness_dimensions(as.matrix(aravo$spe), traits, metric = "euclidean")

str(ui_dim, max.level = 1)

## ----relative-----------------------------------------------------------------
aravo_site_sp = as.matrix(aravo$spe)

# There are clearly abundances and not only presence-absence in this table
aravo_site_sp[1:5, 1:5]

# Compute total abundance per site
site_abundance = rowSums(aravo_site_sp)

head(site_abundance)

# Compute a relative abundance matrix
relative_site_sp = make_relative(aravo_site_sp)

relative_site_sp[1:5, 1:5]

rel_site_abundance = rowSums(relative_site_sp)

head(rel_site_abundance)


