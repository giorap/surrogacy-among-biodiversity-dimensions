##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### DATA ANALYSIS (PRESENT) ##################
#########################################################################
##### -- Phylogenetic/functional metrics as a function of area -- #####
#########################################################################
##### -- Calculate phylogenetic/functional uniqueness across priority areas -- #####
#### Tetrapods
### Calculate evolutionary distinctiveness for all tetrapod species in the Americas
tetrapods_present_evol_distinct <- evol.distinct(tree_tetrapods_complete, type="fair.proportion")
# tetrapods_present_evol_distinct <- readRDS("rds/metrics/tetrapods_present_evol_distinct.rds")
### Calculate trait distinctiveness for all tetrapod species in the Americas
tetrapods_present_trait_distinct <- evol.distinct(traits_tetrapods_tree, type="fair.proportion")
# tetrapods_present_trait_distinct <- readRDS("rds/metrics/tetrapods_present_trait_distinct.rds")

### Distinctiveness of species missing from priority regions
### Evolutionary
## Based on greedy algorithm
tetrapods_present_evol_distinct_greedy <- distinct_by_area(distrib_tetrapods_1, tetrapods_present_evol_distinct, tetrapods_species_present_greedy, algorithm = "greedy")
## Based on CAZ algorithm
tetrapods_present_evol_distinct_CAZ <- distinct_by_area(distrib_tetrapods_1, tetrapods_present_evol_distinct, tetrapods_species_present_CAZ, algorithm = "CAZ")
## Based on ABF algorithm
tetrapods_present_evol_distinct_ABF <- distinct_by_area(distrib_tetrapods_1, tetrapods_present_evol_distinct, tetrapods_species_present_ABF, algorithm = "ABF")
### Trait
## Based on greedy algorithm
tetrapods_present_trait_distinct_greedy <- distinct_by_area(distrib_tetrapods_1, tetrapods_present_trait_distinct, tetrapods_species_present_greedy, algorithm = "greedy")
## Based on CAZ algorithm
tetrapods_present_trait_distinct_CAZ <- distinct_by_area(distrib_tetrapods_1, tetrapods_present_trait_distinct, tetrapods_species_present_CAZ, algorithm = "CAZ")
## Based on ABF algorithm
tetrapods_present_trait_distinct_ABF <- distinct_by_area(distrib_tetrapods_1, tetrapods_present_trait_distinct, tetrapods_species_present_ABF, algorithm = "ABF")

#### Amphibians
### Calculate evolutionary distinctiveness for all tetrapod species in the Americas
amph_present_evol_distinct <- evol.distinct(tree_amph_complete, type="fair.proportion")
### Calculate trait distinctiveness for all tetrapod species in the Americas
amph_present_trait_distinct <- evol.distinct(traits_amph_tree, type="fair.proportion")
### Distinctiveness of species missing from priority regions
### Evolutionary
## Based on greedy algorithm
amph_present_evol_distinct_greedy <- distinct_by_area(distrib_amph_1, amph_present_evol_distinct, amph_species_present_greedy)
## Based on CAZ algorithm
amph_present_evol_distinct_CAZ <- distinct_by_area(distrib_amph_1, amph_present_evol_distinct, amph_species_present_CAZ)
## Based on ABF algorithm
amph_present_evol_distinct_ABF <- distinct_by_area(distrib_amph_1, amph_present_evol_distinct, amph_species_present_ABF)
### Trait
## Based on greedy algorithm
amph_present_trait_distinct_greedy <- distinct_by_area(distrib_amph_1, amph_present_trait_distinct, amph_species_present_greedy)
## Based on CAZ algorithm
amph_present_trait_distinct_CAZ <- distinct_by_area(distrib_amph_1, amph_present_trait_distinct, amph_species_present_CAZ)
## Based on ABF algorithm
amph_present_trait_distinct_ABF <- distinct_by_area(distrib_amph_1, amph_present_trait_distinct, amph_species_present_ABF)

#### Birds
### Calculate evolutionary distinctiveness for all tetrapod species in the Americas
bird_present_evol_distinct <- evol.distinct(tree_bird_complete, type="fair.proportion")
### Calculate trait distinctiveness for all tetrapod species in the Americas
bird_present_trait_distinct <- evol.distinct(traits_bird_tree, type="fair.proportion")
### Distinctiveness of species missing from priority regions
### Evolutionary
## Based on greedy algorithm
bird_present_evol_distinct_greedy <- distinct_by_area(distrib_bird_1, bird_present_evol_distinct, bird_species_present_greedy)
## Based on CAZ algorithm
bird_present_evol_distinct_CAZ <- distinct_by_area(distrib_bird_1, bird_present_evol_distinct, bird_species_present_CAZ)
## Based on ABF algorithm
bird_present_evol_distinct_ABF <- distinct_by_area(distrib_bird_1, bird_present_evol_distinct, bird_species_present_ABF)
### Trait
## Based on greedy algorithm
bird_present_trait_distinct_greedy <- distinct_by_area(distrib_bird_1, bird_present_trait_distinct, bird_species_present_greedy)
## Based on CAZ algorithm
bird_present_trait_distinct_CAZ <- distinct_by_area(distrib_bird_1, bird_present_trait_distinct, bird_species_present_CAZ)
## Based on ABF algorithm
bird_present_trait_distinct_ABF <- distinct_by_area(distrib_bird_1, bird_present_trait_distinct, bird_species_present_ABF)

#### Mammals
### Calculate evolutionary distinctiveness for all tetrapod species in the Americas
mamm_present_evol_distinct <- evol.distinct(tree_mamm_complete, type="fair.proportion")
### Calculate trait distinctiveness for all tetrapod species in the Americas
mamm_present_trait_distinct <- evol.distinct(traits_mamm_tree, type="fair.proportion")
### Distinctiveness of species missing from priority regions
### Evolutionary
## Based on greedy algorithm
mamm_present_evol_distinct_greedy <- distinct_by_area(distrib_mamm_1, mamm_present_evol_distinct, mamm_species_present_greedy)
## Based on CAZ algorithm
mamm_present_evol_distinct_CAZ <- distinct_by_area(distrib_mamm_1, mamm_present_evol_distinct, mamm_species_present_CAZ)
## Based on ABF algorithm
mamm_present_evol_distinct_ABF <- distinct_by_area(distrib_mamm_1, mamm_present_evol_distinct, mamm_species_present_ABF)
### Trait
## Based on greedy algorithm
mamm_present_trait_distinct_greedy <- distinct_by_area(distrib_mamm_1, mamm_present_trait_distinct, mamm_species_present_greedy)
## Based on CAZ algorithm
mamm_present_trait_distinct_CAZ <- distinct_by_area(distrib_mamm_1, mamm_present_trait_distinct, mamm_species_present_CAZ)
## Based on ABF algorithm
mamm_present_trait_distinct_ABF <- distinct_by_area(distrib_mamm_1, mamm_present_trait_distinct, mamm_species_present_ABF)

#### Reptiles
### Calculate evolutionary distinctiveness for all tetrapod species in the Americas
rept_present_evol_distinct <- evol.distinct(tree_rept_complete, type="fair.proportion")
### Calculate trait distinctiveness for all tetrapod species in the Americas
rept_present_trait_distinct <- evol.distinct(traits_rept_tree, type="fair.proportion")
### Distinctiveness of species missing from priority regions
### Evolutionary
## Based on greedy algorithm
rept_present_evol_distinct_greedy <- distinct_by_area(distrib_rept_1, rept_present_evol_distinct, rept_species_present_greedy)
## Based on CAZ algorithm
rept_present_evol_distinct_CAZ <- distinct_by_area(distrib_rept_1, rept_present_evol_distinct, rept_species_present_CAZ)
## Based on ABF algorithm
rept_present_evol_distinct_ABF <- distinct_by_area(distrib_rept_1, rept_present_evol_distinct, rept_species_present_ABF)
### Trait
## Based on greedy algorithm
rept_present_trait_distinct_greedy <- distinct_by_area(distrib_rept_1, rept_present_trait_distinct, rept_species_present_greedy)
## Based on CAZ algorithm
rept_present_trait_distinct_CAZ <- distinct_by_area(distrib_rept_1, rept_present_trait_distinct, rept_species_present_CAZ)
## Based on ABF algorithm
rept_present_trait_distinct_ABF <- distinct_by_area(distrib_rept_1, rept_present_trait_distinct, rept_species_present_ABF)

##### -- Calculate phylogenetic/functional diversity metrics across the study area -- #####
#### Tetrapods
### Phylogenetic diversity
comp_data_tetrapods_phylo <- comparative.comm(tree_tetrapods_complete, as.matrix(distrib_tetrapods_1[, -c(1:2)]))
# Mean Phylogenetic Distance
tetrapods_phylo_mpd <- .mpd(comp_data_tetrapods_phylo)
# Faith's Phylogenetic Diversity
tetrapods_phylo_pd <- .pd(comp_data_tetrapods_phylo)
# Helmus' PSV
tetrapods_phylo_psv <- .psv(comp_data_tetrapods_phylo)
# Combine measures
tetrapods_PD <- data.frame(distrib_tetrapods_1[dimnames(comp_data_tetrapods_phylo$comm)[[1]], c("x", "y")], MPD = tetrapods_phylo_mpd, PD = tetrapods_phylo_pd[, 1], PD.IVS = tetrapods_phylo_pd[, 2], PSV = tetrapods_phylo_psv)
# Save serialized .rds object
saveRDS(tetrapods_PD, "rds/metrics/tetrapods_PD.rds")
### Functional diversity
comp_data_tetrapods_trait <- comparative.comm(traits_tetrapods_tree, as.matrix(distrib_tetrapods_1[, -c(1:2)]))
tetrapods_trait_psv <- .psv(comp_data_tetrapods_trait)
tetrapods_trait_mpd <- .mpd(comp_data_tetrapods_trait)
tetrapods_trait_pd <- .pd(comp_data_tetrapods_trait)
tetrapods_FD <- data.frame(distrib_tetrapods_1[dimnames(comp_data_tetrapods_trait$comm)[[1]], c("x", "y")], MPD = tetrapods_trait_mpd, PD = tetrapods_trait_pd[, 1], PD.IVS = tetrapods_trait_pd[, 2], PSV = tetrapods_trait_psv)
saveRDS(tetrapods_FD, "rds/metrics/tetrapods_FD.rds")

#### Amphibians
### Phylogenetic diversity
comp_data_amph_phylo <- comparative.comm(tree_amph_complete, as.matrix(distrib_amph_1[, -c(1:2)]))
amph_phylo_mpd <- .mpd(comp_data_amph_phylo)
amph_phylo_psv <- .psv(comp_data_amph_phylo)
amph_phylo_pd <- .pd(comp_data_amph_phylo)
amph_PD <- data.frame(distrib_amph_1[dimnames(comp_data_amph_phylo$comm)[[1]], c("x", "y")], MPD = amph_phylo_mpd, PD = amph_phylo_pd[, 1], PD.IVS = amph_phylo_pd[, 2], PSV = amph_phylo_psv)
saveRDS(amph_PD, "rds/metrics/amph_PD.rds")
### Functional diversity
traits_amph_tree$tip.label <- traits_amph_complete$species
comp_data_amph_trait <- comparative.comm(traits_amph_tree, as.matrix(distrib_amph_1[, -c(1:2)]))
amph_trait_psv <- .psv(comp_data_amph_trait)
amph_trait_mpd <- .mpd(comp_data_amph_trait)
amph_trait_pd <- .pd(comp_data_amph_trait)
amph_FD <- data.frame(distrib_amph_1[dimnames(comp_data_amph_trait$comm)[[1]], c("x", "y")], MPD = amph_trait_mpd, PD = amph_trait_pd[, 1], PD.IVS = amph_trait_pd[, 2], PSV = amph_trait_psv)
saveRDS(amph_FD, "rds/metrics/amph_FD.rds")

### Birds
### Phylogenetic diversity
comp_data_bird_phylo <- comparative.comm(tree_bird_complete, as.matrix(distrib_bird_1[, -c(1:2)]))
bird_phylo_mpd <- .mpd(comp_data_bird_phylo)
bird_phylo_psv <- .psv(comp_data_bird_phylo)
bird_phylo_pd <- .pd(comp_data_bird_phylo)
bird_PD <- data.frame(distrib_bird_1[dimnames(comp_data_bird_phylo$comm)[[1]], c("x", "y")], MPD = bird_phylo_mpd, PD = bird_phylo_pd[, 1], PD.IVS = bird_phylo_pd[, 2], PSV = bird_phylo_psv)
saveRDS(bird_PD, "rds/metrics/bird_PD.rds")
### Functional diversity
traits_bird_tree$tip.label <- traits_bird_complete$species
comp_data_bird_trait <- comparative.comm(traits_bird_tree, as.matrix(distrib_bird_1[, -c(1:2)]))
bird_trait_psv <- .psv(comp_data_bird_trait)
bird_trait_mpd <- .mpd(comp_data_bird_trait)
bird_trait_pd <- .pd(comp_data_bird_trait)
bird_FD <- data.frame(distrib_bird_1[dimnames(comp_data_bird_trait$comm)[[1]], c("x", "y")], MPD = bird_trait_mpd, PD = bird_trait_pd[, 1], PD.IVS = bird_trait_pd[, 2], PSV = bird_trait_psv)
saveRDS(bird_FD, "rds/metrics/bird_FD.rds")

### Mammals
### Phylogenetic diversity
comp_data_mamm_phylo <- comparative.comm(tree_mamm_complete, as.matrix(distrib_mamm_1[, -c(1:2)]))
mamm_phylo_mpd <- .mpd(comp_data_mamm_phylo)
mamm_phylo_psv <- .psv(comp_data_mamm_phylo)
mamm_phylo_pd <- .pd(comp_data_mamm_phylo)
mamm_PD <- data.frame(distrib_mamm_1[dimnames(comp_data_mamm_phylo$comm)[[1]], c("x", "y")], MPD = mamm_phylo_mpd, PD = mamm_phylo_pd[, 1], PD.IVS = mamm_phylo_pd[, 2], PSV = mamm_phylo_psv)
saveRDS(mamm_PD, "rds/metrics/mamm_PD.rds")
### Functional diversity
traits_mamm_tree$tip.label <- traits_mamm_complete$species
comp_data_mamm_trait <- comparative.comm(traits_mamm_tree, as.matrix(distrib_mamm_1[, -c(1:2)]))
mamm_trait_psv <- .psv(comp_data_mamm_trait)
mamm_trait_mpd <- .mpd(comp_data_mamm_trait)
mamm_trait_pd <- .pd(comp_data_mamm_trait)
mamm_FD <- data.frame(distrib_mamm_1[dimnames(comp_data_mamm_trait$comm)[[1]], c("x", "y")], MPD = mamm_trait_mpd, PD = mamm_trait_pd[, 1], PD.IVS = mamm_trait_pd[, 2], PSV = mamm_trait_psv)
saveRDS(mamm_FD, "rds/metrics/mamm_FD.rds")

### Reptiles
### Phylogenetic diversity
comp_data_rept_phylo <- comparative.comm(tree_rept_complete, as.matrix(distrib_rept_1[, -c(1:2)]))
rept_phylo_mpd <- .mpd(comp_data_rept_phylo)
rept_phylo_psv <- .psv(comp_data_rept_phylo)
rept_phylo_pd <- .pd(comp_data_rept_phylo)
rept_PD <- data.frame(distrib_rept_1[dimnames(comp_data_rept_phylo$comm)[[1]], c("x", "y")], MPD = rept_phylo_mpd, PD = rept_phylo_pd[, 1], PD.IVS = rept_phylo_pd[, 2], PSV = rept_phylo_psv)
saveRDS(rept_PD, "rds/metrics/rept_PD.rds")
### Functional diversity
traits_rept_tree$tip.label <- traits_rept_complete$species
comp_data_rept_trait <- comparative.comm(traits_rept_tree, as.matrix(distrib_rept_1[, -c(1:2)]))
rept_trait_psv <- .psv(comp_data_rept_trait)
rept_trait_mpd <- .mpd(comp_data_rept_trait)
rept_trait_pd <- .pd(comp_data_rept_trait)
rept_FD <- data.frame(distrib_rept_1[dimnames(comp_data_rept_trait$comm)[[1]], c("x", "y")], MPD = rept_trait_mpd, PD = rept_trait_pd[, 1], PD.IVS = rept_trait_pd[, 2], PSV = rept_trait_psv)
saveRDS(rept_FD, "rds/metrics/rept_FD.rds")