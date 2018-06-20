##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
############# DATA ANALYSIS (FUTURE; DDt Scenario) ########### 
##### -- Spatial prioritizations -- #####
##### -- Greedy algorithm -- #####
##### Runs
#### Tetrapods
### Species
tetrapods_species_DDt_greedy <- greedy_iterations(distrib_tetrapods_DDt_1, runs = 10, nsites = 1000)
# Save/load serialized .rds object
# saveRDS(tetrapods_species_DDt_greedy, "output/tetrapods_species_DDt_greedy.rds")
# tetrapods_species_DDt_greedy <- readRDS("rds/greedy-runs/tetrapods_species_DDt_greedy.rds")
### Trait nodes
tetrapods_trait_nodes_DDt_greedy <- greedy_iterations(traits_tetrapods_DDt_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(tetrapods_phylo_DDt_greedy, "output/tetrapods_phylo_DDt_greedy.rds")
# tetrapods_trait_nodes_DDt_greedy <- readRDS("rds/greedy-runs/tetrapods_trait_nodes_DDt_greedy.rds")
### Phylogenetic nodes
tetrapods_phylo_DDt_greedy <- greedy_iterations(cell_by_node_tetrapods_DDt_1, runs = 10, nsites = 1000)
# saveRDS(tetrapods_phylo_DDt_greedy, "output/tetrapods_phylo_DDt_greedy.rds")
# tetrapods_phylo_DDt_greedy <- readRDS("rds/greedy-runs/tetrapods_phylo_DDt_greedy.rds")
#### Amphibians
### Species
amph_species_DDt_greedy <- greedy_iterations(distrib_amph_DDt_1, runs = 10, nsites = 1000)
# saveRDS(amph_species_DDt_greedy, "output/amph_species_DDt_greedy.rds")
# amph_species_DDt_greedy <- readRDS("rds/greedy-runs/amph_species_DDt_greedy.rds")
### Trait nodes
amph_trait_nodes_DDt_greedy <- greedy_iterations(traits_amph_DDt_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(amph_trait_nodes_DDt_greedy, "output/amph_trait_nodes_DDt_greedyy.rds")
# amph_trait_nodes_DDt_greedy <- readRDS("rds/greedy-runs/amph_trait_nodes_DDt_greedy.rds")
### Phylogenetic nodes
amph_phylo_DDt_greedy <- greedy_iterations(cell_by_node_amph_DDt_1, runs = 10, nsites = 1000)
# saveRDS(amph_phylo_DDt_greedy, "output/amph_phylo_DDt_greedy.rds")
# amph_phylo_DDt_greedy <- readRDS("rds/greedy-runs/amph_phylo_DDt_greedy.rds")
#### Birds
### Species
bird_species_DDt_greedy <- greedy_iterations(distrib_bird_DDt_1, runs = 10, nsites = 1000)
# saveRDS(bird_species_DDt_greedy, "output/bird_species_DDt_greedy.rds")
# bird_species_DDt_greedy <- readRDS("rds/greedy-runs/bird_species_DDt_greedy.rds")
### Trait nodes
bird_trait_nodes_DDt_greedy <- greedy_iterations(traits_bird_DDt_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(bird_trait_nodes_DDt_greedy, "output/bird_trait_nodes_DDt_greedy.rds")
# bird_trait_nodes_DDt_greedy <- readRDS("rds/greedy-runs/bird_trait_nodes_DDt_greedy.rds")
### Phylogenetic nodes
bird_phylo_DDt_greedy <- greedy_iterations(cell_by_node_bird_DDt_1, runs = 10, nsites = 1000)
# saveRDS(bird_phylo_DDt_greedy, "output/bird_phylo_DDt_greedy.rds")
# bird_phylo_DDt_greedy <- readRDS("rds/greedy-runs/bird_phylo_DDt_greedy.rds")
#### Mammals
### Species
mamm_species_DDt_greedy <- greedy_iterations(distrib_mamm_DDt_1, runs = 10, nsites = 1000)
# saveRDS(mamm_species_DDt_greedy, "output/mamm_species_DDt_greedy.rds")
# mamm_species_DDt_greedy <- readRDS("rds/greedy-runs/mamm_species_DDt_greedy.rds")
### Trait nodes
mamm_trait_nodes_DDt_greedy <- greedy_iterations(traits_mamm_DDt_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(mamm_trait_nodes_DDt_greedy, "output/mamm_trait_nodes_DDt_greedy.rds")
# mamm_trait_nodes_DDt_greedy <- readRDS("rds/greedy-runs/mamm_trait_nodes_DDt_greedy.rds")
### Phylogenetic nodes
mamm_phylo_DDt_greedy <- greedy_iterations(cell_by_node_mamm_DDt_1, runs = 10, nsites = 1000)
# saveRDS(mamm_phylo_DDt_greedy, "output/mamm_phylo_DDt_greedy.rds")
# mamm_phylo_DDt_greedy <- readRDS("rds/greedy-runs/mamm_phylo_DDt_greedy.rds")
#### Reptiles
### Species
rept_species_DDt_greedy <- greedy_iterations(distrib_rept_DDt_1, runs = 10, nsites = 1000)
# saveRDS(rept_species_DDt_greedy, "output/rept_species_DDt_greedy.rds")
# rept_species_DDt_greedy <- readRDS("rds/greedy-runs/rept_species_DDt_greedy.rds")
### Trait nodes
rept_trait_nodes_DDt_greedy <- greedy_iterations(traits_rept_DDt_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(rept_trait_nodes_DDt_greedy, "output/rept_trait_nodes_DDt_greedy.rds")
# rept_trait_nodes_DDt_greedy <- readRDS("rds/greedy-runs/rept_trait_nodes_DDt_greedy.rds")
### Phylogenetic nodes
rept_phylo_DDt_greedy <- greedy_iterations(cell_by_node_rept_DDt_1, runs = 10, nsites = 1000)
# saveRDS(rept_phylo_DDt_greedy, "output/rept_phylo_DDt_greedy.rds")
# rept_phylo_DDt_greedy <- readRDS("rds/greedy-runs/rept_phylo_DDt_greedy.rds")
##### Cache data
save.image("cache/analysis-DDt-greedy-run-20180331.RData")

##### -- Zonation analyses -- #####
##### Set up 
#### Set up Zonation software 
#### -- IMPORTANT -- ####
#### For Zonation analyses to run, make sure that you have installed Zonation 4, 
#### a full download can be found at https://www.helsinki.fi/en/researchgroups/metapopulation-research-centre/software,
#### then place the Zonation software's executable file ("zig4.exe") somewhere within your working directory 
#########################
### Install R zonator package, if not already installed
if (!("zonator" %in% installed.packages()[,"Package"])) install.packages("zonator")
### Create .dat template files for running core area Zonation (CAZ) or additive benefit function (ABF)
## Create template files by copying "template.dat" file included with zonator installation
zonator_path <- find.package("zonator")
setwd(zonator_path)
file.copy("extdata/template.dat", "extdata/template_CAZ.dat")
file.copy("extdata/template.dat", "extdata/template_ABF.dat")
## Update removal rule (rule 2 is the additive benefit function) within "template_ABF.dat" file
template_ABF <- readLines("extdata/template_ABF.dat")
template_ABF[2] <- "removal rule = 2"
writeLines(template_ABF, "extdata/template_ABF.dat")

#### Create directory to store outputs
dir.create("output/zonation_runs")
#### Calculate branch lengths to weight phylogenetic node-based Zonation runs
branch_lengths_tetrapods_DDt <- branching.times(tree_tetrapods_DDt)
branch_lengths_amph_DDt <- branching.times(tree_amph_DDt)
branch_lengths_bird_DDt <- branching.times(tree_bird_DDt)
branch_lengths_mamm_DDt <- branching.times(tree_mamm_DDt)
branch_lengths_rept_DDt <- branching.times(tree_rept_DDt)
#### Calculate branch lengths to weight functional node-based Zonation runs
### Tetrapods
trait_branch_lengths_tetrapods_DDt <- branching.times(traits_tetrapods_DDt_tree)
## Re-scale branch lengths by adding the absolute minimum value + 1, to avoid negative or zero branch lengths
## NOTE: this does not affect the relative length of branches among taxa
trait_branch_lengths_tetrapods_DDt <- trait_branch_lengths_tetrapods_DDt + abs(trait_branch_lengths_tetrapods_DDt[which.min(trait_branch_lengths_tetrapods_DDt)]) + 1
### Amphibians
trait_branch_lengths_amph_DDt <- branching.times(traits_amph_DDt_tree)
## Re-scale branch lengths by adding the absolute minimum value + 1, to avoid negative or zero branch lengths
trait_branch_lengths_amph_DDt <- trait_branch_lengths_amph_DDt + abs(trait_branch_lengths_amph_DDt[which.min(trait_branch_lengths_amph_DDt)]) + 1
### Birds
trait_branch_lengths_bird_DDt <- branching.times(traits_bird_DDt_tree)
## Re-scale branch lengths by adding the absolute minimum value + 1, to avoid negative or zero branch lengths
trait_branch_lengths_bird_DDt <- trait_branch_lengths_bird_DDt + abs(trait_branch_lengths_bird_DDt[which.min(trait_branch_lengths_bird_DDt)]) + 1
### Mammals
trait_branch_lengths_mamm_DDt <- branching.times(traits_mamm_DDt_tree)
## Re-scale branch lengths by adding the absolute minimum value + 1, to avoid negative or zero branch lengths
trait_branch_lengths_mamm_DDt <- trait_branch_lengths_mamm_DDt + abs(trait_branch_lengths_mamm_DDt[which.min(trait_branch_lengths_mamm_DDt)]) + 1
### Reptiles
trait_branch_lengths_rept_DDt <- branching.times(traits_rept_DDt_tree)
## Re-scale branch lengths by adding the absolute minimum value + 1, to avoid negative or zero branch lengths
trait_branch_lengths_rept_DDt <- trait_branch_lengths_rept_DDt + abs(trait_branch_lengths_rept_DDt[which.min(trait_branch_lengths_rept_DDt)]) + 1

##### Run zonation
#### Tetrapods
### Species
## CAZ
run_zonation(input = "species", group = "tetrapods", algorithm = "CAZ", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "species", group = "tetrapods", algorithm = "ABF", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "tetrapods", algorithm = "CAZ", period = "DDt", runs = 10, weight = branch_lengths_tetrapods_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "phylo", group = "tetrapods", algorithm = "ABF", period = "DDt", runs = 10, weight = branch_lengths_tetrapods_DDt, zonation_dir = "output/zonation_runs/DDt")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "tetrapods", algorithm = "CAZ", period = "DDt", runs = 10, weight = trait_branch_lengths_tetrapods_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "trait_nodes", group = "tetrapods", algorithm = "ABF", period = "DDt", runs = 10, weight = trait_branch_lengths_tetrapods_DDt, zonation_dir = "output/zonation_runs/DDt")

#### Amphibians
### Species
## CAZ
run_zonation(input = "species", group = "amph", algorithm = "CAZ", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "species", group = "amph", algorithm = "ABF", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "amph", algorithm = "CAZ", period = "DDt", runs = 10, weight = branch_lengths_amph_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "phylo", group = "amph", algorithm = "ABF", period = "DDt", runs = 10, weight = branch_lengths_amph_DDt, zonation_dir = "output/zonation_runs/DDt")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "amph", algorithm = "CAZ", period = "DDt", runs = 10, weight = trait_branch_lengths_amph_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "trait_nodes", group = "amph", algorithm = "ABF", period = "DDt", runs = 10, weight = trait_branch_lengths_amph_DDt, zonation_dir = "output/zonation_runs/DDt")

#### Birds
### Species
## CAZ
run_zonation(input = "species", group = "bird", algorithm = "CAZ", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "species", group = "bird", algorithm = "ABF", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "bird", algorithm = "CAZ", period = "DDt", runs = 10, weight = branch_lengths_bird_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "phylo", group = "bird", algorithm = "ABF", period = "DDt", runs = 10, weight = branch_lengths_bird_DDt, zonation_dir = "output/zonation_runs/DDt")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "bird", algorithm = "CAZ", period = "DDt", runs = 10, weight = trait_branch_lengths_bird_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "trait_nodes", group = "bird", algorithm = "ABF", period = "DDt", runs = 10, weight = trait_branch_lengths_bird_DDt, zonation_dir = "output/zonation_runs/DDt")

#### Mammals
### Species
## CAZ
run_zonation(input = "species", group = "mamm", algorithm = "CAZ", period = "DDt", runs = 100, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "species", group = "mamm", algorithm = "ABF", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "mamm", algorithm = "CAZ", period = "DDt", runs = 10, weight = branch_lengths_mamm_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "phylo", group = "mamm", algorithm = "ABF", period = "DDt", runs = 10, weight = branch_lengths_mamm_DDt, zonation_dir = "output/zonation_runs/DDt")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "mamm", algorithm = "CAZ", period = "DDt", runs = 10, weight = trait_branch_lengths_mamm_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "trait_nodes", group = "mamm", algorithm = "ABF", period = "DDt", runs = 10, weight = trait_branch_lengths_mamm_DDt, zonation_dir = "output/zonation_runs/DDt")

#### Reptiles
### Species
## CAZ
run_zonation(input = "species", group = "rept", algorithm = "CAZ", period = "DDt", runs = 100, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "species", group = "rept", algorithm = "ABF", period = "DDt", runs = 10, zonation_dir = "output/zonation_runs/DDt")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "rept", algorithm = "CAZ", period = "DDt", runs = 10, weight = branch_lengths_rept_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "phylo", group = "rept", algorithm = "ABF", period = "DDt", runs = 10, weight = branch_lengths_rept_DDt, zonation_dir = "output/zonation_runs/DDt")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "rept", algorithm = "CAZ", period = "DDt", runs = 10, weight = trait_branch_lengths_rept_DDt, zonation_dir = "output/zonation_runs/DDt")
## ABF
run_zonation(input = "trait_nodes", group = "rept", algorithm = "ABF", period = "DDt", runs = 10, weight = trait_branch_lengths_rept_DDt, zonation_dir = "output/zonation_runs/DDt")

##### Extract Zonation spatial prioritizations
#### Remember to re-set the current working directory to the main project folder
setwd("C:/Users/Gio/Dropbox/back-up/projects/surrogacy-among-biodiversity-dimensions")
#### Tetrapods
tetrapods_species_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_tetrapods_CAZ", sep = ""), "species_tetrapods_CAZ", runs = 10)
tetrapods_species_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_tetrapods_ABF", sep = ""), "species_tetrapods_ABF", runs = 10)
tetrapods_trait_nodes_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_tetrapods_CAZ", sep = ""), "trait_nodes_tetrapods_CAZ", runs = 10)
tetrapods_trait_nodes_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_tetrapods_ABF", sep = ""), "trait_nodes_tetrapods_ABF", runs = 10)
tetrapods_phylo_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_tetrapods_CAZ", sep = ""), "phylo_tetrapods_CAZ", runs = 10)
tetrapods_phylo_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_tetrapods_ABF", sep = ""), "phylo_tetrapods_ABF", runs = 10)
#### Amphibians
amph_species_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_amph_CAZ", sep = ""), "species_amph_CAZ", runs = 10)
amph_species_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_amph_ABF", sep = ""), "species_amph_ABF", runs = 10)
amph_trait_nodes_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_amph_CAZ", sep = ""), "trait_nodes_amph_CAZ", runs = 10)
amph_trait_nodes_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_amph_ABF", sep = ""), "trait_nodes_amph_ABF", runs = 10)
amph_phylo_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_amph_CAZ", sep = ""), "phylo_amph_CAZ", runs = 10)
amph_phylo_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_amph_ABF", sep = ""), "phylo_amph_ABF", runs = 10)
#### Birds
bird_species_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_bird_CAZ", sep = ""), "species_bird_CAZ", runs = 10)
bird_species_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_bird_ABF", sep = ""), "species_bird_ABF", runs = 10)
bird_trait_nodes_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_bird_CAZ", sep = ""), "trait_nodes_bird_CAZ", runs = 10)
bird_trait_nodes_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_bird_ABF", sep = ""), "trait_nodes_bird_ABF", runs = 10)
bird_phylo_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_bird_CAZ", sep = ""), "phylo_bird_CAZ", runs = 10)
bird_phylo_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_bird_ABF", sep = ""), "phylo_bird_ABF", runs = 10)
#### Mammals
mamm_species_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_mamm_CAZ", sep = ""), "species_mamm_CAZ", runs = 10)
mamm_species_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_mamm_ABF", sep = ""), "species_mamm_ABF", runs = 10)
mamm_trait_nodes_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_mamm_CAZ", sep = ""), "trait_nodes_mamm_CAZ", runs = 10)
mamm_trait_nodes_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_mamm_ABF", sep = ""), "trait_nodes_mamm_ABF", runs = 10)
mamm_phylo_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_mamm_CAZ", sep = ""), "phylo_mamm_CAZ", runs = 10)
mamm_phylo_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_mamm_ABF", sep = ""), "phylo_mamm_ABF", runs = 10)
#### Reptiles
rept_species_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_rept_CAZ", sep = ""), "species_rept_CAZ", runs = 10)
rept_species_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/species_rept_ABF", sep = ""), "species_rept_ABF", runs = 10)
rept_trait_nodes_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_rept_CAZ", sep = ""), "trait_nodes_rept_CAZ", runs = 10)
rept_trait_nodes_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/trait_nodes_rept_ABF", sep = ""), "trait_nodes_rept_ABF", runs = 10)
rept_phylo_DDt_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_rept_CAZ", sep = ""), "phylo_rept_CAZ", runs = 10)
rept_phylo_DDt_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/DDt/phylo_rept_ABF", sep = ""), "phylo_rept_ABF", runs = 10)

##### -- Randomized cell sequences -- #####
#### Tetrapods
### Greedy
tetrapods_DDt_greedy_random <- random_ranks(tetrapods_species_DDt_greedy, algorithm = "greedy", runs = 1000)
### CAZ
tetrapods_DDt_CAZ_random <- random_ranks(tetrapods_species_DDt_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
tetrapods_DDt_ABF_random <- random_ranks(tetrapods_species_DDt_ABF, algorithm = "ABF", runs = 1000)
#### Amphibians
### Greedy
amph_DDt_greedy_random <- random_ranks(amph_species_DDt_greedy, algorithm = "greedy", runs = 1000)
### CAZ
amph_DDt_CAZ_random <- random_ranks(amph_species_DDt_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
amph_DDt_ABF_random <- random_ranks(amph_species_DDt_ABF, algorithm = "ABF", runs = 1000)
#### Birds
### Greedy
bird_DDt_greedy_random <- random_ranks(bird_species_DDt_greedy, algorithm = "greedy", runs = 1000)
### CAZ
bird_DDt_CAZ_random <- random_ranks(bird_species_DDt_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
bird_DDt_ABF_random <- random_ranks(bird_species_DDt_ABF, algorithm = "ABF", runs = 1000)
#### Mammals
### Greedy
mamm_DDt_greedy_random <- random_ranks(mamm_species_DDt_greedy, algorithm = "greedy", runs = 1000)
### CAZ
mamm_DDt_CAZ_random <- random_ranks(mamm_species_DDt_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
mamm_DDt_ABF_random <- random_ranks(mamm_species_DDt_ABF, algorithm = "ABF", runs = 1000)
#### Reptiles
### Greedy
rept_DDt_greedy_random <- random_ranks(rept_species_DDt_greedy, algorithm = "greedy", runs = 1000)
### CAZ
rept_DDt_CAZ_random <- random_ranks(rept_species_DDt_CAZ, algorithm = "CAZ", runs = 1000)
### ABF
rept_DDt_ABF_random <- random_ranks(rept_species_DDt_ABF, algorithm = "ABF", runs = 1000)

##### -- Surrogacy analyses -- #####
##### -- Generate target accumulation curves -- #####
#### Tetrapods
#### Species as a surrogate for trait target
#### Trait nodes
### Greedy
tetrapods_trait_nodes_DDt_greedy_curves <- get_sai_curves(target_matrix = traits_tetrapods_DDt_cell_by_node, target_ranks = tetrapods_trait_nodes_DDt_greedy, surrogate_ranks = tetrapods_species_DDt_greedy, random_ranks = tetrapods_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(tetrapods_trait_nodes_DDt_greedy_curves, "rds/curves/tetrapods_trait_nodes_DDt_greedy_curves.rds")
# tetrapods_trait_nodes_DDt_greedy_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDt_greedy_curves.rds")
### CAZ
tetrapods_trait_nodes_DDt_CAZ_curves <- get_sai_curves(target_matrix = traits_tetrapods_DDt_cell_by_node, target_ranks = tetrapods_trait_nodes_DDt_CAZ, surrogate_ranks = tetrapods_species_DDt_CAZ, random_ranks = tetrapods_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", branch_lengths = trait_branch_lengths_tetrapods_DDt)
# saveRDS(tetrapods_trait_nodes_DDt_CAZ_curves, "rds/curves/tetrapods_trait_nodes_DDt_CAZ_curves.rds")
# tetrapods_trait_nodes_DDt_CAZ_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDt_CAZ_curves.rds")
### ABF
tetrapods_trait_nodes_DDt_ABF_curves <- get_sai_curves(target_matrix = traits_tetrapods_DDt_cell_by_node, target_ranks = tetrapods_trait_nodes_DDt_ABF, surrogate_ranks = tetrapods_species_DDt_ABF, random_ranks = tetrapods_DDt_ABF_random, target = "phylo", algorithm = "ABF", branch_lengths = trait_branch_lengths_tetrapods_DDt)
# saveRDS(tetrapods_trait_nodes_DDt_ABF_curves, "rds/curves/tetrapods_trait_nodes_DDt_ABF_curves.rds")
# tetrapods_trait_nodes_DDt_ABF_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDt_ABF_curves.rds")

#### Species as a surrogate for phylo target
### Greedy
tetrapods_phylo_DDt_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_tetrapods_1, target_ranks = tetrapods_phylo_DDt_greedy, surrogate_ranks = tetrapods_species_DDt_greedy, random_ranks = tetrapods_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(tetrapods_phylo_DDt_greedy_curves, "rds/curves/tetrapods_phylo_DDt_greedy_curves.rds")
# tetrapods_phylo_DDt_greedy_curves <- readRDS("rds/curves/tetrapods_phylo_DDt_greedy_curves.rds")
### CAZ
tetrapods_phylo_DDt_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_tetrapods_1, target_ranks = tetrapods_phylo_DDt_CAZ, surrogate_ranks = tetrapods_species_DDt_CAZ, random_ranks = tetrapods_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_tetrapods_complete)
# saveRDS(tetrapods_phylo_DDt_CAZ_curves, "rds/curves/tetrapods_phylo_DDt_CAZ_curves.rds")
# tetrapods_phylo_DDt_CAZ_curves <- readRDS("rds/curves/tetrapods_phylo_DDt_CAZ_curves.rds")
### ABF
tetrapods_phylo_DDt_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_tetrapods_1, target_ranks = tetrapods_phylo_DDt_ABF, surrogate_ranks = tetrapods_species_DDt_ABF, random_ranks = tetrapods_DDt_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_tetrapods_complete)
# saveRDS(tetrapods_phylo_DDt_ABF_curves, "rds/curves/tetrapods_phylo_DDt_ABF_curves.rds")
# tetrapods_phylo_DDt_ABF_curves <- readRDS("rds/curves/tetrapods_phylo_DDt_ABF_curves.rds")

#### Amphibians
#### Species as a surrogate for trait target
#### Trait nodes
### Greedy
amph_trait_nodes_DDt_greedy_curves <- get_sai_curves(target_matrix = traits_amph_cell_by_node, target_ranks = amph_trait_nodes_DDt_greedy, surrogate_ranks = amph_species_DDt_greedy, random_ranks = amph_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(amph_trait_nodes_DDt_greedy_curves, "rds/curves/amph_trait_nodes_DDt_greedy_curves.rds")
# amph_trait_nodes_DDt_greedy_curves <- readRDS("rds/curves/amph_trait_nodes_DDt_greedy_curves.rds")
### CAZ
amph_trait_nodes_DDt_CAZ_curves <- get_sai_curves(target_matrix = traits_amph_cell_by_node, target_ranks = amph_trait_nodes_DDt_CAZ, surrogate_ranks = amph_species_DDt_CAZ, random_ranks = amph_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = traits_amph_tree)
# saveRDS(amph_trait_nodes_DDt_CAZ_curves, "rds/curves/amph_trait_nodes_DDt_CAZ_curves.rds")
# amph_trait_nodes_DDt_CAZ_curves <- readRDS("rds/curves/amph_trait_nodes_DDt_CAZ_curves.rds")
### ABF
amph_trait_nodes_DDt_ABF_curves <- get_sai_curves(target_matrix = traits_amph_cell_by_node, target_ranks = amph_trait_nodes_DDt_ABF, surrogate_ranks = amph_species_DDt_ABF, random_ranks = amph_DDt_ABF_random, target = "phylo", algorithm = "ABF", tree = traits_amph_tree)
# saveRDS(amph_trait_nodes_DDt_ABF_curves, "rds/curves/amph_trait_nodes_DDt_ABF_curves.rds")
# amph_trait_nodes_DDt_ABF_curves <- readRDS("rds/curves/amph_trait_nodes_DDt_ABF_curves.rds")

#### Species as a surrogate for phylo target
### Greedy
amph_phylo_DDt_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_amph_1, target_ranks = amph_phylo_DDt_greedy, surrogate_ranks = amph_species_DDt_greedy, random_ranks = amph_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(amph_phylo_DDt_greedy_curves, "rds/curves/amph_phylo_DDt_greedy_curves.rds")
# amph_phylo_DDt_greedy_curves <- readRDS("rds/curves/amph_phylo_DDt_greedy_curves.rds")
### CAZ
amph_phylo_DDt_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_amph_1, target_ranks = amph_phylo_DDt_CAZ, surrogate_ranks = amph_species_DDt_CAZ, random_ranks = amph_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_amph_complete)
# saveRDS(amph_phylo_DDt_CAZ_curves, "rds/curves/amph_phylo_DDt_CAZ_curves.rds")
# amph_phylo_DDt_CAZ_curves <- readRDS("rds/curves/amph_phylo_DDt_CAZ_curves.rds")
### ABF
amph_phylo_DDt_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_amph_1, target_ranks = amph_phylo_DDt_ABF, surrogate_ranks = amph_species_DDt_ABF, random_ranks = amph_DDt_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_amph_complete)
# saveRDS(amph_phylo_DDt_ABF_curves, "rds/curves/amph_phylo_DDt_ABF_curves.rds")
# amph_phylo_DDt_ABF_curves <- readRDS("rds/curves/amph_phylo_DDt_ABF_curves.rds")

#### Birds
#### Species as a surrogate for trait target
#### Trait nodes
### Greedy
bird_trait_nodes_DDt_greedy_curves <- get_sai_curves(target_matrix = traits_bird_cell_by_node, target_ranks = bird_trait_nodes_DDt_greedy, surrogate_ranks = bird_species_DDt_greedy, random_ranks = bird_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(bird_trait_nodes_DDt_greedy_curves, "rds/curves/bird_trait_nodes_DDt_greedy_curves.rds")
# bird_trait_nodes_DDt_greedy_curves <- readRDS("rds/curves/bird_trait_nodes_DDt_greedy_curves.rds")
### CAZ
bird_trait_nodes_DDt_CAZ_curves <- get_sai_curves(target_matrix = traits_bird_cell_by_node, target_ranks = bird_trait_nodes_DDt_CAZ, surrogate_ranks = bird_species_DDt_CAZ, random_ranks = bird_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = traits_bird_tree)
# saveRDS(bird_trait_nodes_DDt_CAZ_curves, "rds/curves/bird_trait_nodes_DDt_CAZ_curves.rds")
# bird_trait_nodes_DDt_CAZ_curves <- readRDS("rds/curves/bird_trait_nodes_DDt_CAZ_curves.rds")
### ABF
bird_trait_nodes_DDt_ABF_curves <- get_sai_curves(target_matrix = traits_bird_cell_by_node, target_ranks = bird_trait_nodes_DDt_ABF, surrogate_ranks = bird_species_DDt_ABF, random_ranks = bird_DDt_ABF_random, target = "phylo", algorithm = "ABF", tree = traits_bird_tree)
# saveRDS(bird_trait_nodes_DDt_ABF_curves, "rds/curves/bird_trait_nodes_DDt_ABF_curves.rds")
# bird_trait_nodes_DDt_ABF_curves <- readRDS("rds/curves/bird_trait_nodes_DDt_ABF_curves.rds")

#### Species as a surrogate for phylo target
### Greedy
bird_phylo_DDt_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_bird_1, target_ranks = bird_phylo_DDt_greedy, surrogate_ranks = bird_species_DDt_greedy, random_ranks = bird_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(bird_phylo_DDt_greedy_curves, "rds/curves/bird_phylo_DDt_greedy_curves.rds")
# bird_phylo_DDt_greedy_curves <- readRDS("rds/curves/bird_phylo_DDt_greedy_curves.rds")
### CAZ
bird_phylo_DDt_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_bird_1, target_ranks = bird_phylo_DDt_CAZ, surrogate_ranks = bird_species_DDt_CAZ, random_ranks = bird_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_bird_complete)
# saveRDS(bird_phylo_DDt_CAZ_curves, "rds/curves/bird_phylo_DDt_CAZ_curves.rds")
# bird_phylo_DDt_CAZ_curves <- readRDS("rds/curves/bird_phylo_DDt_CAZ_curves.rds")
### ABF
bird_phylo_DDt_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_bird_1, target_ranks = bird_phylo_DDt_ABF, surrogate_ranks = bird_species_DDt_ABF, random_ranks = bird_DDt_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_bird_complete)
# saveRDS(bird_phylo_DDt_ABF_curves, "rds/curves/bird_phylo_DDt_ABF_curves.rds")
# bird_phylo_DDt_ABF_curves <- readRDS("rds/curves/bird_phylo_DDt_ABF_curves.rds")


#### Mammals
#### Species as a surrogate for trait target
#### Trait nodes
### Greedy
mamm_trait_nodes_DDt_greedy_curves <- get_sai_curves(target_matrix = traits_mamm_cell_by_node, target_ranks = mamm_trait_nodes_DDt_greedy, surrogate_ranks = mamm_species_DDt_greedy, random_ranks = mamm_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(mamm_trait_nodes_DDt_greedy_curves, "rds/curves/mamm_trait_nodes_DDt_greedy_curves.rds")
# mamm_trait_nodes_DDt_greedy_curves <- readRDS("rds/curves/mamm_trait_nodes_DDt_greedy_curves.rds")
### CAZ
mamm_trait_nodes_DDt_CAZ_curves <- get_sai_curves(target_matrix = traits_mamm_cell_by_node, target_ranks = mamm_trait_nodes_DDt_CAZ, surrogate_ranks = mamm_species_DDt_CAZ, random_ranks = mamm_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = traits_mamm_tree)
# saveRDS(mamm_trait_nodes_DDt_CAZ_curves, "rds/curves/mamm_trait_nodes_DDt_CAZ_curves.rds")
# mamm_trait_nodes_DDt_CAZ_curves <- readRDS("rds/curves/mamm_trait_nodes_DDt_CAZ_curves.rds")
### ABF
mamm_trait_nodes_DDt_ABF_curves <- get_sai_curves(target_matrix = traits_mamm_cell_by_node, target_ranks = mamm_trait_nodes_DDt_ABF, surrogate_ranks = mamm_species_DDt_ABF, random_ranks = mamm_DDt_ABF_random, target = "phylo", algorithm = "ABF", tree = traits_mamm_tree)
# saveRDS(mamm_trait_nodes_DDt_ABF_curves, "rds/curves/mamm_trait_nodes_DDt_ABF_curves.rds")
# mamm_trait_nodes_DDt_ABF_curves <- readRDS("rds/curves/mamm_trait_nodes_DDt_ABF_curves.rds")

#### Species as a surrogate for phylo target
### Greedy
mamm_phylo_DDt_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_mamm_1, target_ranks = mamm_phylo_DDt_greedy, surrogate_ranks = mamm_species_DDt_greedy, random_ranks = mamm_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(mamm_phylo_DDt_greedy_curves, "rds/curves/mamm_phylo_DDt_greedy_curves.rds")
# mamm_phylo_DDt_greedy_curves <- readRDS("rds/curves/mamm_phylo_DDt_greedy_curves.rds")
### CAZ
mamm_phylo_DDt_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_mamm_1, target_ranks = mamm_phylo_DDt_CAZ, surrogate_ranks = mamm_species_DDt_CAZ, random_ranks = mamm_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_mamm_complete)
# saveRDS(mamm_phylo_DDt_CAZ_curves, "rds/curves/mamm_phylo_DDt_CAZ_curves.rds")
# mamm_phylo_DDt_CAZ_curves <- readRDS("rds/curves/mamm_phylo_DDt_CAZ_curves.rds")
### ABF
mamm_phylo_DDt_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_mamm_1, target_ranks = mamm_phylo_DDt_ABF, surrogate_ranks = mamm_species_DDt_ABF, random_ranks = mamm_DDt_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_mamm_complete)
# saveRDS(mamm_phylo_DDt_ABF_curves, "rds/curves/mamm_phylo_DDt_ABF_curves.rds")
# mamm_phylo_DDt_ABF_curves <- readRDS("rds/curves/mamm_phylo_DDt_ABF_curves.rds")

#### Reptiles
#### Species as a surrogate for trait target
#### Trait nodes
### Greedy
rept_trait_nodes_DDt_greedy_curves <- get_sai_curves(target_matrix = traits_rept_cell_by_node, target_ranks = rept_trait_nodes_DDt_greedy, surrogate_ranks = rept_species_DDt_greedy, random_ranks = rept_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(rept_trait_nodes_DDt_greedy_curves, "rds/curves/rept_trait_nodes_DDt_greedy_curves.rds")
# rept_trait_nodes_DDt_greedy_curves <- readRDS("rds/curves/rept_trait_nodes_DDt_greedy_curves.rds")
### CAZ
rept_trait_nodes_DDt_CAZ_curves <- get_sai_curves(target_matrix = traits_rept_cell_by_node, target_ranks = rept_trait_nodes_DDt_CAZ, surrogate_ranks = rept_species_DDt_CAZ, random_ranks = rept_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = traits_rept_tree)
# saveRDS(rept_trait_nodes_DDt_CAZ_curves, "rds/curves/rept_trait_nodes_DDt_CAZ_curves.rds")
# rept_trait_nodes_DDt_CAZ_curves <- readRDS("rds/curves/rept_trait_nodes_DDt_CAZ_curves.rds")
### ABF
rept_trait_nodes_DDt_ABF_curves <- get_sai_curves(target_matrix = traits_rept_cell_by_node, target_ranks = rept_trait_nodes_DDt_ABF, surrogate_ranks = rept_species_DDt_ABF, random_ranks = rept_DDt_ABF_random, target = "phylo", algorithm = "ABF", tree = traits_rept_tree)
# saveRDS(rept_trait_nodes_DDt_ABF_curves, "rds/curves/rept_trait_nodes_DDt_ABF_curves.rds")
# rept_trait_nodes_DDt_ABF_curves <- readRDS("rds/curves/rept_trait_nodes_DDt_ABF_curves.rds")

#### Species as a surrogate for phylo target
### Greedy
rept_phylo_DDt_greedy_curves <- get_sai_curves(target_matrix = cell_by_node_rept_1, target_ranks = rept_phylo_DDt_greedy, surrogate_ranks = rept_species_DDt_greedy, random_ranks = rept_DDt_greedy_random, target = "phylo", algorithm = "greedy")
# saveRDS(rept_phylo_DDt_greedy_curves, "rds/curves/rept_phylo_DDt_greedy_curves.rds")
# rept_phylo_DDt_greedy_curves <- readRDS("rds/curves/rept_phylo_DDt_greedy_curves.rds")
### CAZ
rept_phylo_DDt_CAZ_curves <- get_sai_curves(target_matrix = cell_by_node_rept_1, target_ranks = rept_phylo_DDt_CAZ, surrogate_ranks = rept_species_DDt_CAZ, random_ranks = rept_DDt_CAZ_random, target = "phylo", algorithm = "CAZ", tree = tree_rept_complete)
# saveRDS(rept_phylo_DDt_CAZ_curves, "rds/curves/rept_phylo_DDt_CAZ_curves.rds")
# rept_phylo_DDt_CAZ_curves <- readRDS("rds/curves/rept_phylo_DDt_CAZ_curves.rds")
### ABF
rept_phylo_DDt_ABF_curves <- get_sai_curves(target_matrix = cell_by_node_rept_1, target_ranks = rept_phylo_DDt_ABF, surrogate_ranks = rept_species_DDt_ABF, random_ranks = rept_DDt_CAZ_random, target = "phylo", algorithm = "ABF", tree = tree_rept_complete)
# saveRDS(rept_phylo_DDt_ABF_curves, "rds/curves/rept_phylo_DDt_ABF_curves.rds")
# rept_phylo_DDt_ABF_curves <- readRDS("rds/curves/rept_phylo_DDt_ABF_curves.rds")

##### -- Calculate Species Accumulation Index (SAI) -- #####
#### Tetrapods
### Trait nodes
tetrapods_trait_nodes_DDt_greedy_sai <- calculate_sai(tetrapods_trait_nodes_DDt_greedy_curves)
tetrapods_trait_nodes_DDt_CAZ_sai <- calculate_sai(tetrapods_trait_nodes_DDt_CAZ_curves)
tetrapods_trait_nodes_DDt_ABF_sai <- calculate_sai(tetrapods_trait_nodes_DDt_ABF_curves)
### Phylo
tetrapods_phylo_DDt_greedy_sai <- calculate_sai(tetrapods_phylo_DDt_greedy_curves)
tetrapods_phylo_DDt_CAZ_sai <- calculate_sai(tetrapods_phylo_DDt_CAZ_curves)
tetrapods_phylo_DDt_ABF_sai <- calculate_sai(tetrapods_phylo_DDt_ABF_curves)

#### Amphibians
### Trait nodes
amph_trait_nodes_DDt_greedy_sai <- calculate_sai(amph_trait_nodes_DDt_greedy_curves)
amph_trait_nodes_DDt_CAZ_sai <- calculate_sai(amph_trait_nodes_DDt_CAZ_curves)
amph_trait_nodes_DDt_ABF_sai <- calculate_sai(amph_trait_nodes_DDt_ABF_curves)
### Phylo
amph_phylo_DDt_greedy_sai <- calculate_sai(amph_phylo_DDt_greedy_curves)
amph_phylo_DDt_CAZ_sai <- calculate_sai(amph_phylo_DDt_CAZ_curves)
amph_phylo_DDt_ABF_sai <- calculate_sai(amph_phylo_DDt_ABF_curves)

#### Birds
### Trait nodes
bird_trait_nodes_DDt_greedy_sai <- calculate_sai(bird_trait_nodes_DDt_greedy_curves)
bird_trait_nodes_DDt_CAZ_sai <- calculate_sai(bird_trait_nodes_DDt_CAZ_curves)
bird_trait_nodes_DDt_ABF_sai <- calculate_sai(bird_trait_nodes_DDt_ABF_curves)
### Phylo
bird_phylo_DDt_greedy_sai <- calculate_sai(bird_phylo_DDt_greedy_curves)
bird_phylo_DDt_CAZ_sai <- calculate_sai(bird_phylo_DDt_CAZ_curves)
bird_phylo_DDt_ABF_sai <- calculate_sai(bird_phylo_DDt_ABF_curves)

#### Mammals
### Trait nodes
mamm_trait_nodes_DDt_greedy_sai <- calculate_sai(mamm_trait_nodes_DDt_greedy_curves)
mamm_trait_nodes_DDt_CAZ_sai <- calculate_sai(mamm_trait_nodes_DDt_CAZ_curves)
mamm_trait_nodes_DDt_ABF_sai <- calculate_sai(mamm_trait_nodes_DDt_ABF_curves)
### Phylo
mamm_phylo_DDt_greedy_sai <- calculate_sai(mamm_phylo_DDt_greedy_curves)
mamm_phylo_DDt_CAZ_sai <- calculate_sai(mamm_phylo_DDt_CAZ_curves)
mamm_phylo_DDt_ABF_sai <- calculate_sai(mamm_phylo_DDt_ABF_curves)

#### Reptiles
### Trait nodes
rept_trait_nodes_DDt_greedy_sai <- calculate_sai(rept_trait_nodes_DDt_greedy_curves)
rept_trait_nodes_DDt_CAZ_sai <- calculate_sai(rept_trait_nodes_DDt_CAZ_curves)
rept_trait_nodes_DDt_ABF_sai <- calculate_sai(rept_trait_nodes_DDt_ABF_curves)
### Phylo
rept_phylo_DDt_greedy_sai <- calculate_sai(rept_phylo_DDt_greedy_curves)
rept_phylo_DDt_CAZ_sai <- calculate_sai(rept_phylo_DDt_CAZ_curves)
rept_phylo_DDt_ABF_sai <- calculate_sai(rept_phylo_DDt_ABF_curves)
