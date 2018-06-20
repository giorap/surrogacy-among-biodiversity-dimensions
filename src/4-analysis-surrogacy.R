##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### DATA ANALYSIS (PRESENT) ##################
####################################
##### -- Surrogacy analyses -- #####
####################################
##### -- Randomized cell sequences -- #####
##### Generate 1000 random cell sequences, this is done for each taxon x algorithm combination
#### Tetrapods
### Greedy
tetrapods_present_greedy_random <- random_ranks(tetrapods_species_present_greedy, algorithm = "greedy", runs = 1000)
### CAZ
tetrapods_present_CAZ_random <- random_ranks(tetrapods_species_present_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
tetrapods_present_ABF_random <- random_ranks(tetrapods_species_present_ABF, algorithm = "ABF", runs = 1000)
#### Amphibians
### Greedy
amph_present_greedy_random <- random_ranks(amph_species_present_greedy, algorithm = "greedy", runs = 1000)
### CAZ
amph_present_CAZ_random <- random_ranks(amph_species_present_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
amph_present_ABF_random <- random_ranks(amph_species_present_ABF, algorithm = "ABF", runs = 1000)
#### Birds
### Greedy
bird_present_greedy_random <- random_ranks(bird_species_present_greedy, algorithm = "greedy", runs = 1000)
### CAZ
bird_present_CAZ_random <- random_ranks(bird_species_present_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
bird_present_ABF_random <- random_ranks(bird_species_present_ABF, algorithm = "ABF", runs = 1000)
#### Mammals
### Greedy
mamm_present_greedy_random <- random_ranks(mamm_species_present_greedy, algorithm = "greedy", runs = 1000)
### CAZ
mamm_present_CAZ_random <- random_ranks(mamm_species_present_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
mamm_present_ABF_random <- random_ranks(mamm_species_present_ABF, algorithm = "ABF", runs = 1000)
#### Reptiles
### Greedy
rept_present_greedy_random <- random_ranks(rept_species_present_greedy, algorithm = "greedy", runs = 1000)
### CAZ/ABF
rept_present_CAZ_random <- random_ranks(rept_species_present_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
rept_present_ABF_random <- random_ranks(rept_species_present_ABF, algorithm = "ABF", runs = 1000)

##### -- Target accumulation curves -- #####
#### Tetrapods
#### Species as a surrogate for phylo target
### Greedy
tetrapods_phylo_present_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_tetrapods_1, target_ranks = tetrapods_phylo_present_greedy, surrogate_ranks = tetrapods_species_present_greedy, random_ranks = tetrapods_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(tetrapods_phylo_present_greedy_curves, "rds/curves/tetrapods_phylo_present_greedy_curves.rds")
# tetrapods_phylo_present_greedy_curves <- readRDS("rds/curves/tetrapods_phylo_present_greedy_curves.rds")
### CAZ
tetrapods_phylo_present_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_tetrapods_1, target_ranks = tetrapods_phylo_present_CAZ, surrogate_ranks = tetrapods_species_present_CAZ, random_ranks = tetrapods_present_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_tetrapods_complete)
# saveRDS(tetrapods_phylo_present_CAZ_curves, "rds/curves/tetrapods_phylo_present_CAZ_curves.rds")
# tetrapods_phylo_present_CAZ_curves <- readRDS("rds/curves/tetrapods_phylo_present_CAZ_curves.rds")
### ABF
tetrapods_phylo_present_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_tetrapods_1, target_ranks = tetrapods_phylo_present_ABF, surrogate_ranks = tetrapods_species_present_ABF, random_ranks = tetrapods_present_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_tetrapods_complete)
# saveRDS(tetrapods_phylo_present_ABF_curves, "rds/curves/tetrapods_phylo_present_ABF_curves.rds")
# tetrapods_phylo_present_ABF_curves <- readRDS("rds/curves/tetrapods_phylo_present_ABF_curves.rds")
#### Species as a surrogate for trait target
### Greedy
tetrapods_trait_nodes_present_greedy_curves <- get_sai_curves(target_matrix = traits_tetrapods_cell_by_node, target_ranks = tetrapods_trait_nodes_present_greedy, surrogate_ranks = tetrapods_species_present_greedy, random_ranks = tetrapods_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(tetrapods_trait_nodes_present_greedy_curves, "rds/curves/tetrapods_trait_nodes_present_greedy_curves.rds")
# tetrapods_trait_nodes_present_greedy_curves <- readRDS("rds/curves/tetrapods_trait_nodes_present_greedy_curves.rds")
### CAZ
tetrapods_trait_nodes_present_CAZ_curves <- get_sai_curves(target_matrix = traits_tetrapods_cell_by_node, target_ranks = tetrapods_trait_nodes_present_CAZ, surrogate_ranks = tetrapods_species_present_CAZ, random_ranks = tetrapods_present_CAZ_random, target = "phylo", algorithm = "CAZ", branch_lengths = trait_branch_lengths_tetrapods_present)
# saveRDS(tetrapods_trait_nodes_present_CAZ_curves, "rds/curves/tetrapods_trait_nodes_present_CAZ_curves.rds")
# tetrapods_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/tetrapods_trait_nodes_present_CAZ_curves.rds")
### ABF
tetrapods_trait_nodes_present_ABF_curves <- get_sai_curves(target_matrix = traits_tetrapods_cell_by_node, target_ranks = tetrapods_trait_nodes_present_ABF, surrogate_ranks = tetrapods_species_present_ABF, random_ranks = tetrapods_present_ABF_random, target = "phylo", algorithm = "ABF", branch_lengths = trait_branch_lengths_tetrapods_present)
# saveRDS(tetrapods_trait_nodes_present_ABF_curves, "rds/curves/tetrapods_trait_nodes_present_ABF_curves.rds")
# tetrapods_trait_nodes_present_ABF_curves <- readRDS("rds/curves/tetrapods_trait_nodes_present_ABF_curves.rds")

#### Amphibians
#### Species as a surrogate for phylo target
### Greedy
amph_phylo_present_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_amph_1, target_ranks = amph_phylo_present_greedy, surrogate_ranks = amph_species_present_greedy, random_ranks = amph_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(amph_phylo_present_greedy_curves, "rds/curves/amph_phylo_present_greedy_curves.rds")
# amph_phylo_present_greedy_curves <- readRDS("rds/curves/amph_phylo_present_greedy_curves.rds")
### CAZ
amph_phylo_present_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_amph_1, target_ranks = amph_phylo_present_CAZ, surrogate_ranks = amph_species_present_CAZ, random_ranks = amph_present_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_amph_complete)
# saveRDS(amph_phylo_present_CAZ_curves, "rds/curves/amph_phylo_present_CAZ_curves.rds")
# amph_phylo_present_CAZ_curves <- readRDS("rds/curves/amph_phylo_present_CAZ_curves.rds")
### ABF
amph_phylo_present_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_amph_1, target_ranks = amph_phylo_present_ABF, surrogate_ranks = amph_species_present_ABF, random_ranks = amph_present_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_amph_complete)
# saveRDS(amph_phylo_present_ABF_curves, "rds/curves/amph_phylo_present_ABF_curves.rds")
# amph_phylo_present_ABF_curves <- readRDS("rds/curves/amph_phylo_present_ABF_curves.rds")
#### Species as a surrogate for trait target
#### Trait nodes
### Greedy
amph_trait_nodes_present_greedy_curves <- get_sai_curves(target_matrix = traits_amph_cell_by_node, target_ranks = amph_trait_nodes_present_greedy, surrogate_ranks = amph_species_present_greedy, random_ranks = amph_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(amph_trait_nodes_present_greedy_curves, "rds/curves/amph_trait_nodes_present_greedy_curves.rds")
# amph_trait_nodes_present_greedy_curves <- readRDS("rds/curves/amph_trait_nodes_present_greedy_curves.rds")
### CAZ
amph_trait_nodes_present_CAZ_curves <- get_sai_curves(target_matrix = traits_amph_cell_by_node, target_ranks = amph_trait_nodes_present_CAZ, surrogate_ranks = amph_species_present_CAZ, random_ranks = amph_present_CAZ_random, target = "phylo", algorithm = "CAZ", branch_lengths = trait_branch_lengths_amph_present)
# saveRDS(amph_trait_nodes_present_CAZ_curves, "rds/curves/amph_trait_nodes_present_CAZ_curves.rds")
# amph_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/amph_trait_nodes_present_CAZ_curves.rds")
### ABF
amph_trait_nodes_present_ABF_curves <- get_sai_curves(target_matrix = traits_amph_cell_by_node, target_ranks = amph_trait_nodes_present_ABF, surrogate_ranks = amph_species_present_ABF, random_ranks = amph_present_ABF_random, target = "phylo", algorithm = "ABF", branch_lengths = trait_branch_lengths_amph_present)
# saveRDS(amph_trait_nodes_present_ABF_curves, "rds/curves/amph_trait_nodes_present_ABF_curves.rds")
# amph_trait_nodes_present_ABF_curves <- readRDS("rds/curves/amph_trait_nodes_present_ABF_curves.rds")

#### Birds
#### Species as a surrogate for phylo target
### Greedy
bird_phylo_present_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_bird_1, target_ranks = bird_phylo_present_greedy, surrogate_ranks = bird_species_present_greedy, random_ranks = bird_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(bird_phylo_present_greedy_curves, "rds/curves/bird_phylo_present_greedy_curves.rds")
# bird_phylo_present_greedy_curves <- readRDS("rds/curves/bird_phylo_present_greedy_curves.rds")
### CAZ
bird_phylo_present_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_bird_1, target_ranks = bird_phylo_present_CAZ, surrogate_ranks = bird_species_present_CAZ, random_ranks = bird_present_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_bird_complete)
# saveRDS(bird_phylo_present_CAZ_curves, "rds/curves/bird_phylo_present_CAZ_curves.rds")
# bird_phylo_present_CAZ_curves <- readRDS("rds/curves/bird_phylo_present_CAZ_curves.rds")
### ABF
bird_phylo_present_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_bird_1, target_ranks = bird_phylo_present_ABF, surrogate_ranks = bird_species_present_ABF, random_ranks = bird_present_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_bird_complete)
# saveRDS(bird_phylo_present_ABF_curves, "rds/curves/bird_phylo_present_ABF_curves.rds")
# bird_phylo_present_ABF_curves <- readRDS("rds/curves/bird_phylo_present_ABF_curves.rds")
#### Species as a surrogate for trait target
### Greedy
bird_trait_nodes_present_greedy_curves <- get_sai_curves(target_matrix = traits_bird_cell_by_node, target_ranks = bird_trait_nodes_present_greedy, surrogate_ranks = bird_species_present_greedy, random_ranks = bird_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(bird_trait_nodes_present_greedy_curves, "rds/curves/bird_trait_nodes_present_greedy_curves.rds")
# bird_trait_nodes_present_greedy_curves <- readRDS("rds/curves/bird_trait_nodes_present_greedy_curves.rds")
### CAZ
bird_trait_nodes_present_CAZ_curves <- get_sai_curves(target_matrix = traits_bird_cell_by_node, target_ranks = bird_trait_nodes_present_CAZ, surrogate_ranks = bird_species_present_CAZ, random_ranks = bird_present_CAZ_random, target = "phylo", algorithm = "CAZ", branch_lengths = trait_branch_lengths_bird_present)
# saveRDS(bird_trait_nodes_present_CAZ_curves, "rds/curves/bird_trait_nodes_present_CAZ_curves.rds")
# bird_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/bird_trait_nodes_present_CAZ_curves.rds")
### ABF
bird_trait_nodes_present_ABF_curves <- get_sai_curves(target_matrix = traits_bird_cell_by_node, target_ranks = bird_trait_nodes_present_ABF, surrogate_ranks = bird_species_present_ABF, random_ranks = bird_present_ABF_random, target = "phylo", algorithm = "ABF", branch_lengths = trait_branch_lengths_bird_present)
# saveRDS(bird_trait_nodes_present_ABF_curves, "rds/curves/bird_trait_nodes_present_ABF_curves.rds")
# bird_trait_nodes_present_ABF_curves <- readRDS("rds/curves/bird_trait_nodes_present_ABF_curves.rds")

#### Mammals
#### Species as a surrogate for phylo target
### Greedy
mamm_phylo_present_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_mamm_1, target_ranks = mamm_phylo_present_greedy, surrogate_ranks = mamm_species_present_greedy, random_ranks = mamm_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(mamm_phylo_present_greedy_curves, "rds/curves/mamm_phylo_present_greedy_curves.rds")
# mamm_phylo_present_greedy_curves <- readRDS("rds/curves/mamm_phylo_present_greedy_curves.rds")
### CAZ
mamm_phylo_present_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_mamm_1, target_ranks = mamm_phylo_present_CAZ, surrogate_ranks = mamm_species_present_CAZ, random_ranks = mamm_present_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_mamm_complete)
# saveRDS(mamm_phylo_present_CAZ_curves, "rds/curves/mamm_phylo_present_CAZ_curves.rds")
# mamm_phylo_present_CAZ_curves <- readRDS("rds/curves/mamm_phylo_present_CAZ_curves.rds")
### ABF
mamm_phylo_present_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_mamm_1, target_ranks = mamm_phylo_present_ABF, surrogate_ranks = mamm_species_present_ABF, random_ranks = mamm_present_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_mamm_complete)
# saveRDS(mamm_phylo_present_ABF_curves, "rds/curves/mamm_phylo_present_ABF_curves.rds")
# mamm_phylo_present_ABF_curves <- readRDS("rds/curves/mamm_phylo_present_ABF_curves.rds")
#### Species as a surrogate for trait target
### Greedy
mamm_trait_nodes_present_greedy_curves <- get_sai_curves(target_matrix = traits_mamm_cell_by_node, target_ranks = mamm_trait_nodes_present_greedy, surrogate_ranks = mamm_species_present_greedy, random_ranks = mamm_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(mamm_trait_nodes_present_greedy_curves, "rds/curves/mamm_trait_nodes_present_greedy_curves.rds")
# mamm_trait_nodes_present_greedy_curves <- readRDS("rds/curves/mamm_trait_nodes_present_greedy_curves.rds")
### CAZ
mamm_trait_nodes_present_CAZ_curves <- get_sai_curves(target_matrix = traits_mamm_cell_by_node, target_ranks = mamm_trait_nodes_present_CAZ, surrogate_ranks = mamm_species_present_CAZ, random_ranks = mamm_present_CAZ_random, target = "phylo", algorithm = "CAZ", branch_lengths = trait_branch_lengths_mamm_present)
# saveRDS(mamm_trait_nodes_present_CAZ_curves, "rds/curves/mamm_trait_nodes_present_CAZ_curves.rds")
# mamm_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/mamm_trait_nodes_present_CAZ_curves.rds")
### ABF
mamm_trait_nodes_present_ABF_curves <- get_sai_curves(target_matrix = traits_mamm_cell_by_node, target_ranks = mamm_trait_nodes_present_ABF, surrogate_ranks = mamm_species_present_ABF, random_ranks = mamm_present_ABF_random, target = "phylo", algorithm = "ABF", branch_lengths = trait_branch_lengths_mamm_present)
# saveRDS(mamm_trait_nodes_present_ABF_curves, "rds/curves/mamm_trait_nodes_present_ABF_curves.rds")
# mamm_trait_nodes_present_ABF_curves <- readRDS("rds/curves/mamm_trait_nodes_present_ABF_curves.rds")

#### Reptiles
#### Species as a surrogate for phylo target
### Greedy
rept_phylo_present_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_rept_1, target_ranks = rept_phylo_present_greedy, surrogate_ranks = rept_species_present_greedy, random_ranks = rept_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(rept_phylo_present_greedy_curves, "rds/curves/rept_phylo_present_greedy_curves.rds")
# rept_phylo_present_greedy_curves <- readRDS("rds/curves/rept_phylo_present_greedy_curves.rds")
### CAZ
rept_phylo_present_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_rept_1, target_ranks = rept_phylo_present_CAZ, surrogate_ranks = rept_species_present_CAZ, random_ranks = rept_present_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_rept_complete)
# saveRDS(rept_phylo_present_CAZ_curves, "rds/curves/rept_phylo_present_CAZ_curves.rds")
# rept_phylo_present_CAZ_curves <- readRDS("rds/curves/rept_phylo_present_CAZ_curves.rds")
### ABF
rept_phylo_present_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_rept_1, target_ranks = rept_phylo_present_ABF, surrogate_ranks = rept_species_present_ABF, random_ranks = rept_present_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_rept_complete)
# saveRDS(rept_phylo_present_ABF_curves, "rds/curves/rept_phylo_present_ABF_curves.rds")
# rept_phylo_present_ABF_curves <- readRDS("rds/curves/rept_phylo_present_ABF_curves.rds")
#### Species as a surrogate for trait target
### Greedy
rept_trait_nodes_present_greedy_curves <- get_sai_curves(target_matrix = traits_rept_cell_by_node, target_ranks = rept_trait_nodes_present_greedy, surrogate_ranks = rept_species_present_greedy, random_ranks = rept_present_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(rept_trait_nodes_present_greedy_curves, "rds/curves/rept_trait_nodes_present_greedy_curves.rds")
# rept_trait_nodes_present_greedy_curves <- readRDS("rds/curves/rept_trait_nodes_present_greedy_curves.rds")
### CAZ
rept_trait_nodes_present_CAZ_curves <- get_sai_curves(target_matrix = traits_rept_cell_by_node, target_ranks = rept_trait_nodes_present_CAZ, surrogate_ranks = rept_species_present_CAZ, random_ranks = rept_present_CAZ_random, target = "phylo", algorithm = "CAZ", branch_lengths = trait_branch_lengths_rept_present)
# saveRDS(rept_trait_nodes_present_CAZ_curves, "rds/curves/rept_trait_nodes_present_CAZ_curves.rds")
# rept_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/rept_trait_nodes_present_CAZ_curves.rds")
### ABF
rept_trait_nodes_present_ABF_curves <- get_sai_curves(target_matrix = traits_rept_cell_by_node, target_ranks = rept_trait_nodes_present_ABF, surrogate_ranks = rept_species_present_ABF, random_ranks = rept_present_ABF_random, target = "phylo", algorithm = "ABF", branch_lengths = trait_branch_lengths_rept_present)
# saveRDS(rept_trait_nodes_present_ABF_curves, "rds/curves/rept_trait_nodes_present_ABF_curves.rds")
# rept_trait_nodes_present_ABF_curves <- readRDS("rds/curves/rept_trait_nodes_present_ABF_curves.rds")

##### -- Calculate Species Accumulation Index (SAI) -- #####
#### Tetrapods
### Phylo
tetrapods_phylo_present_greedy_sai <- calculate_sai(tetrapods_phylo_present_greedy_curves)
tetrapods_phylo_present_CAZ_sai <- calculate_sai(tetrapods_phylo_present_CAZ_curves)
tetrapods_phylo_present_ABF_sai <- calculate_sai(tetrapods_phylo_present_ABF_curves)
### Trait nodes
tetrapods_trait_nodes_present_greedy_sai <- calculate_sai(tetrapods_trait_nodes_present_greedy_curves)
tetrapods_trait_nodes_present_CAZ_sai <- calculate_sai(tetrapods_trait_nodes_present_CAZ_curves)
tetrapods_trait_nodes_present_ABF_sai <- calculate_sai(tetrapods_trait_nodes_present_ABF_curves)

#### Amphibians
### Phylo
amph_phylo_present_greedy_sai <- calculate_sai(amph_phylo_present_greedy_curves)
amph_phylo_present_CAZ_sai <- calculate_sai(amph_phylo_present_CAZ_curves)
amph_phylo_present_ABF_sai <- calculate_sai(amph_phylo_present_ABF_curves)
### Trait nodes
amph_trait_nodes_present_greedy_sai <- calculate_sai(amph_trait_nodes_present_greedy_curves)
amph_trait_nodes_present_CAZ_sai <- calculate_sai(amph_trait_nodes_present_CAZ_curves)
amph_trait_nodes_present_ABF_sai <- calculate_sai(amph_trait_nodes_present_ABF_curves)

#### Birds
### Phylo
bird_phylo_present_greedy_sai <- calculate_sai(bird_phylo_present_greedy_curves)
bird_phylo_present_CAZ_sai <- calculate_sai(bird_phylo_present_CAZ_curves)
bird_phylo_present_ABF_sai <- calculate_sai(bird_phylo_present_ABF_curves)
### Trait nodes
bird_trait_nodes_present_greedy_sai <- calculate_sai(bird_trait_nodes_present_greedy_curves)
bird_trait_nodes_present_CAZ_sai <- calculate_sai(bird_trait_nodes_present_CAZ_curves)
bird_trait_nodes_present_ABF_sai <- calculate_sai(bird_trait_nodes_present_ABF_curves)

#### Mammals
### Phylo
mamm_phylo_present_greedy_sai <- calculate_sai(mamm_phylo_present_greedy_curves)
mamm_phylo_present_CAZ_sai <- calculate_sai(mamm_phylo_present_CAZ_curves)
mamm_phylo_present_ABF_sai <- calculate_sai(mamm_phylo_present_ABF_curves)
### Trait nodes
mamm_trait_nodes_present_greedy_sai <- calculate_sai(mamm_trait_nodes_present_greedy_curves)
mamm_trait_nodes_present_CAZ_sai <- calculate_sai(mamm_trait_nodes_present_CAZ_curves)
mamm_trait_nodes_present_ABF_sai <- calculate_sai(mamm_trait_nodes_present_ABF_curves)

#### Reptiles
### Phylo
rept_phylo_present_greedy_sai <- calculate_sai(rept_phylo_present_greedy_curves)
rept_phylo_present_CAZ_sai <- calculate_sai(rept_phylo_present_CAZ_curves)
rept_phylo_present_ABF_sai <- calculate_sai(rept_phylo_present_ABF_curves)
### Trait nodes
rept_trait_nodes_present_greedy_sai <- calculate_sai(rept_trait_nodes_present_greedy_curves)
rept_trait_nodes_present_CAZ_sai <- calculate_sai(rept_trait_nodes_present_CAZ_curves)
rept_trait_nodes_present_ABF_sai <- calculate_sai(rept_trait_nodes_present_ABF_curves)
