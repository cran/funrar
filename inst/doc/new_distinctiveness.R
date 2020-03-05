## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
species_df = data.frame(
  species =     c("a", "b", "c"),
  trait_value = c(-1, 0, 0.5)
)

species_distance = dist(c(a = -1, b = 0, c = 0.5))

species_distance

## ----new_distinctiveness------------------------------------------------------
alternative_distinctiveness = function(pres_mat, distance_obj, given_T) {
  dist_mat = as.matrix(distance_obj)
  
  kept_sp = funrar:::species_in_common(pres_mat, dist_mat)
  dist_mat = dist_mat[kept_sp, kept_sp, drop = FALSE]
  
  
  # Correspondence matrix (tracking which species we want to keep)
  corr_dist = dist_mat
  corr_dist[dist_mat > given_T] = 0
  corr_dist[dist_mat <= given_T] = 1
  diag(corr_dist) = 0
  
  # Across sites
  # given_pres is a vector of presences per site
  di_mat = apply(pres_mat, 1, function(given_pres) {
    index_mat = given_pres %*% (dist_mat * corr_dist)  # Sum of distances per species
    denom_mat = given_pres %*% corr_dist  # Number of distances considered per species
    index_mat = index_mat / denom_mat
    index_mat[given_pres == 0] = NA
    index_mat[is.nan(index_mat)] = 1
    
    return(index_mat)
  })
  
  di_mat = t(di_mat)
  dimnames(di_mat) = dimnames(pres_mat)
  
  dplyr::mutate(funrar::matrix_to_stack(di_mat, "Di"), given_range = given_T)
}

## -----------------------------------------------------------------------------
presence_matrix = matrix(c(rep(1, 3), 1, 0, 1, 1, 1, 0), nrow = 3, ncol = 3,
                         dimnames = list(site = c("s1", "s2", "s3"),
                                         species = c("a", "b", "c")))

# Test over a ragne of possible T values
all_T = lapply(seq(0.5, 1.5, length.out = 50),
       function(given_number) alternative_distinctiveness(presence_matrix,
                                                          species_distance,
                                                          given_number))

all_T = dplyr::bind_rows(all_T)

## -----------------------------------------------------------------------------
library(ggplot2)

ggplot(all_T, aes(given_range, Di, color = species)) +
    geom_line(size = 1, alpha = 1/2) +
    facet_grid(~site) +
  labs(x = "Fixed distance range",
       y = "Functional Distinctiveness",
       color = "Species")

## ----aravo_new_distinct-------------------------------------------------------
library("funrar")

data("aravo", package = "ade4")
# Site-species matrix
mat = as.matrix(aravo$spe)

# Convert matrix to presence-absence matrix
mat[mat > 0] = 1

# Example of trait table
tra = aravo$traits[, c("Height", "SLA", "N_mass")]
# Distance matrix
dist_mat = compute_dist_matrix(tra, metric = "gower")
dist_mat = (dist_mat - min(dist_mat))/diff(range(dist_mat))

names(dimnames(mat)) = c("site", "species")

# Compute different values of distinctiveness using various ranges
all_ranges = lapply(seq(0, 1, length.out = 50), function(given_range) {
  alternative_distinctiveness(mat, as.dist(dist_mat), given_range)
})

all_ranges = dplyr::bind_rows(all_ranges)

## -----------------------------------------------------------------------------
ggplot(subset(all_ranges, site %in% c("AR07", "AR51", "AR02")),
       aes(given_range, Di, group = species)) +
  geom_line(alpha = 1/3) +
  facet_wrap(~site) +
  labs(x = "Maximum Distance Range Considered\n(Trait Range)",
       y = "Functional Distinctiveness")

## ----scaled_range-------------------------------------------------------------
all_ranges = dplyr::mutate(all_ranges, scaled_Di =
                             ifelse(Di != 1, Di / given_range, Di))

ggplot(subset(all_ranges, site %in% c("AR07", "AR51", "AR02")),
       aes(given_range, scaled_Di, group = species)) +
  geom_line(alpha = 1/4) +
  facet_wrap(~site) +
  labs(x = "Considered Trait Range\n(Functional Distance)",
       y = "Scaled Functional Distinctiveness\n(over trait range)")

## ----new_distinctiveness_ab---------------------------------------------------
ab_mat = matrix(c(rep(1/3, 3), 1/6, 1/6, 4/6, 4/6, 1/6, 1/6), nrow = 3, ncol = 3,
                         dimnames = list(site = c("s1", "s2", "s3"),
                                         species = c("a", "b", "c")),
                byrow = TRUE)

alternative_distinctiveness_abundance = function(abund_mat, dist_matrix,
                                                 given_range) {
  
  dist_mat = dist_matrix
  
  kept_sp = funrar:::species_in_common(abund_mat, dist_mat)
  dist_mat = dist_mat[kept_sp, kept_sp, drop = FALSE]
  
  
  # Correspondence matrix (tracking which species we want to keep)
  corr_dist = dist_mat
  corr_dist[dist_mat > given_range] = 0
  corr_dist[dist_mat <= given_range] = 1
  diag(corr_dist) = 0
  
  # Across sites
  # given_pres is a vector of presences per site
  di_mat = apply(abund_mat, 1, function(given_ab) {
    index_mat = given_ab %*% (dist_mat * corr_dist)  # Sum of distances per species
    denom_mat = given_ab %*% corr_dist  # Number of distances considered per species
    index_mat = (index_mat / denom_mat) * (1 - denom_mat)
    index_mat[given_ab == 0 | is.na(given_ab)] = NA
    index_mat[is.nan(index_mat)] = 1
    
    return(index_mat)
  })
  
  di_mat = t(di_mat)
  dimnames(di_mat) = dimnames(abund_mat)
  
  dplyr::mutate(funrar::matrix_to_stack(di_mat, "Di"),
                given_range = given_range)
}

## -----------------------------------------------------------------------------
ab_di_all_ranges = lapply(seq(0, 1.5, length.out = 50),
       function(given_number) alternative_distinctiveness_abundance(ab_mat,
                                                          as.matrix(species_distance),
                                                          given_number))

ab_di_all_ranges = dplyr::bind_rows(ab_di_all_ranges)

ggplot(ab_di_all_ranges, aes(given_range, Di, color = species)) +
  geom_line(size = 1) +
  facet_wrap(~site, labeller = as_labeller(c(s1 = "s1 (1/3 rel. abund each)",
                                             s2 = "s2 (a=1/6, b=1/6, c=4/6)",
                                             s3 = "s3 (a=4/6, b=1/6, c=1/6)"))) +
  labs(x = "Considered Range",
       y = "Functional Distinctiveness")

