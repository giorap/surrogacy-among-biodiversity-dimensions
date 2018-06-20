##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### DATA ANALYSIS (PRESENT) ##################
#########################################
##### -- Spatial prioritizations -- #####
#########################################
##### -- Greedy algorithm -- #####
##### Runs
##### Output from greedy runs saved as serialized .rds objects in folder "rds/greedy-runs/"
#### Create folder
dir.create("rds/greedy-runs")
#### Tetrapods
### Species
tetrapods_species_present_greedy <- greedy_iterations(distrib_tetrapods_1, runs = 10, nsites = 1000)
# Save/load serialized .rds object
# saveRDS(tetrapods_species_present_greedy, "rds/greedy-runs/tetrapods_species_present_greedy.rds")
# tetrapods_species_present_greedy <- readRDS("rds/greedy-runs/tetrapods_species_present_greedy.rds")
### Phylogenetic branches
tetrapods_phylo_present_greedy <- greedy_iterations(cell_by_node_tetrapods_1, runs = 10, nsites = 1000)
# saveRDS(tetrapods_phylo_present_greedy, "rds/greedy-runs/tetrapods_phylo_present_greedy.rds")
# tetrapods_phylo_present_greedy <- readRDS("rds/greedy-runs/tetrapods_traits_present_greedy.rds")
### Trait branches
tetrapods_trait_nodes_present_greedy <- greedy_iterations(traits_tetrapods_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(tetrapods_trait_nodes_present_greedy, "rds/greedy-runs/tetrapods_trait_nodes_present_greedy.rds")
# tetrapods_trait_nodes_present_greedy <- readRDS("rds/greedy-runs/tetrapods_trait_nodes_present_greedy.rds")

#### Amphibians
### Species
amph_species_present_greedy <- greedy_iterations(distrib_amph_1, runs = 10, nsites = 1000)
# saveRDS(amph_species_present_greedy, "rds/greedy-runs/amph_species_present_greedy.rds")
# amph_species_present_greedy <- readRDS("rds/greedy-runs/amph_species_present_greedy.rds")
### Phylogenetic nodes
amph_phylo_present_greedy <- greedy_iterations(cell_by_node_amph_1, runs = 10, nsites = 1000)
# saveRDS(amph_phylo_present_greedy, "rds/greedy-runs/amph_phylo_present_greedy.rds")
# amph_phylo_present_greedy <- readRDS("rds/greedy-runs/amph_phylo_present_greedy.rds")
### Trait nodes
amph_trait_nodes_present_greedy <- greedy_iterations(traits_amph_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(amph_trait_nodes_present_greedy, "rds/greedy-runs/amph_trait_nodes_present_greedy.rds")
# amph_trait_nodes_present_greedy <- readRDS("rds/greedy-runs/amph_trait_nodes_present_greedy.rds")

#### Birds
### Species
bird_species_present_greedy <- greedy_iterations(distrib_bird_1, runs = 10, nsites = 1000)
# saveRDS(bird_species_present_greedy, "rds/greedy-runs/bird_species_present_greedy.rds")
# bird_species_present_greedy <- readRDS("rds/greedy-runs/bird_species_present_greedy.rds")
### Phylogenetic nodes
bird_phylo_present_greedy <- greedy_iterations(cell_by_node_bird_1, runs = 10, nsites = 1000)
# saveRDS(bird_phylo_present_greedy, "rds/greedy-runs/bird_phylo_present_greedy.rds")
# bird_phylo_present_greedy <- readRDS("rds/greedy-runs/bird_phylo_present_greedy.rds")
## Trait nodes
bird_trait_nodes_present_greedy <- greedy_iterations(traits_bird_cell_by_node, runs = 10, nsites = 1000)
# saveRDS(bird_trait_nodes_present_greedy, "rds/greedy-runs/bird_trait_nodes_present_greedy.rds")
# bird_trait_nodes_present_greedy <- readRDS("rds/greedy-runs/bird_trait_nodes_present_greedy.rds")

#### Mammals
### Species
mamm_species_present_greedy <- greedy_iterations(distrib_mamm_1, runs = 10, nsites = 1000)
# saveRDS(mamm_species_present_greedy, "rds/greedy-runs/mamm_species_present_greedy.rds")
# mamm_species_present_greedy <- readRDS("rds/greedy-runs/mamm_species_present_greedy.rds")
### Phylogenetic nodes
mamm_phylo_present_greedy <- greedy_iterations(cell_by_node_mamm_1, runs = 10, nsites = 1000)
# saveRDS(mamm_phylo_present_greedy, "rds/greedy-runs/mamm_phylo_present_greedy.rds")
# mamm_phylo_present_greedy <- readRDS("rds/greedy-runs/mamm_phylo_present_greedy.rds")
## Trait nodes
mamm_trait_nodes_present_greedy <- greedy_iterations(traits_mamm_cell_by_node, runs = 10, nsites = 1000)
saveRDS(mamm_trait_nodes_present_greedy, "mamm_trait_nodes_present_greedy.rds")
# saveRDS(mamm_trait_nodes_present_greedy, "rds/greedy-runs/mamm_trait_nodes_present_greedy.rds")
# mamm_trait_nodes_present_greedy <- readRDS("rds/greedy-runs/mamm_trait_nodes_present_greedy.rds")

#### Reptiles
### Species
rept_species_present_greedy <- greedy_iterations(distrib_rept_1, runs = 10, nsites = 1000)
# saveRDS(rept_species_present_greedy, "rds/greedy-runs/rept_species_present_greedy.rds")
# rept_species_present_greedy <- readRDS("rds/greedy-runs/rept_species_present_greedy.rds")
### Phylogenetic nodes
rept_phylo_present_greedy <- greedy_iterations(cell_by_node_rept_1, runs = 10, nsites = 1000)
# saveRDS(rept_phylo_present_greedy, "rds/greedy-runs/rept_phylo_present_greedy.rds")
# rept_phylo_present_greedy <- readRDS("rds/greedy-runs/rept_phylo_present_greedy.rds")
## Trait nodes
rept_trait_nodes_present_greedy <- greedy_iterations(traits_rept_cell_by_node, runs = 10, nsites = 1000)
saveRDS(rept_trait_nodes_present_greedy, "rept_trait_nodes_present_greedy.rds")
# saveRDS(rept_trait_nodes_present_greedy, "rds/greedy-runs/rept_trait_nodes_present_greedy.rds")
# rept_trait_nodes_present_greedy <- readRDS("rds/greedy-runs/rept_trait_nodes_present_greedy.rds")

##### -- Zonation-based prioritizations -- #####
##### Set up Zonation software 
##### -- IMPORTANT -- ##################################################################################################
#### For Zonation analyses to run, make sure that you have installed Zonation 4, 
#### a full download can be found at https://www.helsinki.fi/en/researchgroups/metapopulation-research-centre/software,
#### then place the Zonation software's executable file ("zig4.exe") somewhere within your working directory 
########################################################################################################################
##### Install R zonator package, if not already installed
if (!("zonator" %in% installed.packages()[,"Package"])) install.packages("zonator")
### Create .dat template files for running core area Zonation (CAZ) or additive benefit function (ABF)
## Create template files by copying "template.dat" file included with zonator installation
zonator_path <- find.package("zonator")
setwd(zonator_path)
file.copy("extdata/template.dat", "extdata/template_CAZ.dat")
file.copy("extdata/template.dat", "extdata/template_ABF.dat")
## Update removal rule (rule 2 is the Additive Benefit Function) within "template_ABF.dat" file
template_ABF <- readLines("extdata/template_ABF.dat")
template_ABF[2] <- "removal rule = 2"
writeLines(template_ABF, "extdata/template_ABF.dat")

#### Create directory to store outputs
dir.create("output/zonation_runs")

#### Calculate branch lengths to weight phylogenetic/functional node-based Zonation runs
### Phylogenetic branch lenghts
phylo_branch_lengths_tetrapods_present <- branching.times(tree_tetrapods_complete)
phylo_branch_lengths_amph_present <- branching.times(tree_amph_complete)
phylo_branch_lengths_bird_present <- branching.times(tree_bird_complete)
phylo_branch_lengths_mamm_present <- branching.times(tree_mamm_complete)
phylo_branch_lengths_rept_present <- branching.times(tree_rept_complete)
### Functional branch lenghts
### Tetrapods
trait_branch_lengths_tetrapods_present <- branching.times(traits_tetrapods_tree)
trait_branch_lengths_amph_present <- branching.times(traits_amph_tree)
trait_branch_lengths_bird_present <- branching.times(traits_bird_tree)
trait_branch_lengths_mamm_present <- branching.times(traits_mamm_tree)
trait_branch_lengths_rept_present <- branching.times(traits_rept_tree)

## Re-scale branch lengths by adding the absolute minimum value + 1, to avoid negative or zero weights
## This is necessary otherwise negative/zero weights will lead to particular branches being ignored within spatial prioritizations
## NOTE: this does not affect the relative length of branches among taxa and, thus, their relative weight within analyses
trait_branch_lengths_tetrapods_present <- trait_branch_lengths_tetrapods_present + abs(trait_branch_lengths_tetrapods_present[which.min(trait_branch_lengths_tetrapods_present)]) + 1
trait_branch_lengths_amph_present <- trait_branch_lengths_amph_present + abs(trait_branch_lengths_amph_present[which.min(trait_branch_lengths_amph_present)]) + 1
trait_branch_lengths_bird_present <- trait_branch_lengths_bird_present + abs(trait_branch_lengths_bird_present[which.min(trait_branch_lengths_bird_present)]) + 1
trait_branch_lengths_mamm_present <- trait_branch_lengths_mamm_present + abs(trait_branch_lengths_mamm_present[which.min(trait_branch_lengths_mamm_present)]) + 1
trait_branch_lengths_rept_present <- trait_branch_lengths_rept_present + abs(trait_branch_lengths_rept_present[which.min(trait_branch_lengths_rept_present)]) + 1

##### Run zonation
#### Tetrapods
### Species
## CAZ
run_zonation(input = "species", group = "tetrapods", algorithm = "CAZ", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "species", group = "tetrapods", algorithm = "ABF", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "tetrapods", algorithm = "CAZ", period = "present", runs = 10, weight = branch_lengths_tetrapods_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "phylo", group = "tetrapods", algorithm = "ABF", period = "present", runs = 10, weight = branch_lengths_tetrapods_present, zonation_dir = "output/zonation_runs/present")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "tetrapods", algorithm = "CAZ", period = "present", runs = 10, weight = trait_branch_lengths_tetrapods_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "trait_nodes", group = "tetrapods", algorithm = "ABF", period = "present", runs = 10, weight = trait_branch_lengths_tetrapods_present, zonation_dir = "output/zonation_runs/present")

#### Amphibians
### Species
## CAZ
run_zonation(input = "species", group = "amph", algorithm = "CAZ", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "species", group = "amph", algorithm = "ABF", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "amph", algorithm = "CAZ", period = "present", runs = 10, weight = branch_lengths_amph_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "phylo", group = "amph", algorithm = "ABF", period = "present", runs = 10, weight = branch_lengths_amph_present, zonation_dir = "output/zonation_runs/present")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "amph", algorithm = "CAZ", period = "present", runs = 10, weight = trait_branch_lengths_amph_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "trait_nodes", group = "amph", algorithm = "ABF", period = "present", runs = 10, weight = trait_branch_lengths_amph_present, zonation_dir = "output/zonation_runs/present")

#### Birds
### Species
## CAZ
run_zonation(input = "species", group = "bird", algorithm = "CAZ", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "species", group = "bird", algorithm = "ABF", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "bird", algorithm = "CAZ", period = "present", runs = 10, weight = branch_lengths_bird_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "phylo", group = "bird", algorithm = "ABF", period = "present", runs = 10, weight = branch_lengths_bird_present, zonation_dir = "output/zonation_runs/present")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "bird", algorithm = "CAZ", period = "present", runs = 10, weight = trait_branch_lengths_bird_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "trait_nodes", group = "bird", algorithm = "ABF", period = "present", runs = 10, weight = trait_branch_lengths_bird_present, zonation_dir = "output/zonation_runs/present")

#### Mammals
### Species
## CAZ
run_zonation(input = "species", group = "mamm", algorithm = "CAZ", period = "present", runs = 100, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "species", group = "mamm", algorithm = "ABF", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "mamm", algorithm = "CAZ", period = "present", runs = 10, weight = branch_lengths_mamm_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "phylo", group = "mamm", algorithm = "ABF", period = "present", runs = 10, weight = branch_lengths_mamm_present, zonation_dir = "output/zonation_runs/present")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "mamm", algorithm = "CAZ", period = "present", runs = 10, weight = trait_branch_lengths_mamm_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "trait_nodes", group = "mamm", algorithm = "ABF", period = "present", runs = 10, weight = trait_branch_lengths_mamm_present, zonation_dir = "output/zonation_runs/present")

#### Reptiles
### Species
## CAZ
run_zonation(input = "species", group = "rept", algorithm = "CAZ", period = "present", runs = 100, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "species", group = "rept", algorithm = "ABF", period = "present", runs = 10, zonation_dir = "output/zonation_runs/present")
### Phylogenetic nodes
## CAZ
run_zonation(input = "phylo", group = "rept", algorithm = "CAZ", period = "present", runs = 10, weight = branch_lengths_rept_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "phylo", group = "rept", algorithm = "ABF", period = "present", runs = 10, weight = branch_lengths_rept_present, zonation_dir = "output/zonation_runs/present")
### Trait nodes
## CAZ
run_zonation(input = "trait_nodes", group = "rept", algorithm = "CAZ", period = "present", runs = 10, weight = trait_branch_lengths_rept_present, zonation_dir = "output/zonation_runs/present")
## ABF
run_zonation(input = "trait_nodes", group = "rept", algorithm = "ABF", period = "present", runs = 10, weight = trait_branch_lengths_rept_present, zonation_dir = "output/zonation_runs/present")

##### Extract Zonation spatial prioritizations
#### Remember to re-set the current working directory to the main project folder
setwd("C:/Users/Gio/Dropbox/back-up/projects/surrogacy-among-biodiversity-dimensions")
#### Tetrapods
tetrapods_species_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_tetrapods_CAZ", sep = ""), "species_tetrapods_CAZ", runs = 10)
tetrapods_species_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_tetrapods_ABF", sep = ""), "species_tetrapods_ABF", runs = 10)
tetrapods_phylo_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_tetrapods_CAZ", sep = ""), "phylo_tetrapods_CAZ", runs = 10)
tetrapods_phylo_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_tetrapods_ABF", sep = ""), "phylo_tetrapods_ABF", runs = 10)
tetrapods_trait_nodes_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_tetrapods_CAZ", sep = ""), "trait_nodes_tetrapods_CAZ", runs = 10)
tetrapods_trait_nodes_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_tetrapods_ABF", sep = ""), "trait_nodes_tetrapods_ABF", runs = 10)

#### Amphibians
amph_species_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_amph_CAZ", sep = ""), "species_amph_CAZ", runs = 10)
amph_species_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_amph_ABF", sep = ""), "species_amph_ABF", runs = 10)
amph_phylo_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_amph_CAZ", sep = ""), "phylo_amph_CAZ", runs = 10)
amph_phylo_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_amph_ABF", sep = ""), "phylo_amph_ABF", runs = 10)
amph_trait_nodes_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_amph_CAZ", sep = ""), "trait_nodes_amph_CAZ", runs = 10)
amph_trait_nodes_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_amph_ABF", sep = ""), "trait_nodes_amph_ABF", runs = 10)

#### Birds
bird_species_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_bird_CAZ", sep = ""), "species_bird_CAZ", runs = 10)
bird_species_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_bird_ABF", sep = ""), "species_bird_ABF", runs = 10)
bird_phylo_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_bird_CAZ", sep = ""), "phylo_bird_CAZ", runs = 10)
bird_phylo_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_bird_ABF", sep = ""), "phylo_bird_ABF", runs = 10)
bird_trait_nodes_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_bird_CAZ", sep = ""), "trait_nodes_bird_CAZ", runs = 10)
bird_trait_nodes_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_bird_ABF", sep = ""), "trait_nodes_bird_ABF", runs = 10)

#### Mammals
mamm_species_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_mamm_CAZ", sep = ""), "species_mamm_CAZ", runs = 10)
mamm_species_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_mamm_ABF", sep = ""), "species_mamm_ABF", runs = 10)
mamm_phylo_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_mamm_CAZ", sep = ""), "phylo_mamm_CAZ", runs = 10)
mamm_phylo_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_mamm_ABF", sep = ""), "phylo_mamm_ABF", runs = 10)
mamm_trait_nodes_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_mamm_CAZ", sep = ""), "trait_nodes_mamm_CAZ", runs = 10)
mamm_trait_nodes_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_mamm_ABF", sep = ""), "trait_nodes_mamm_ABF", runs = 10)

#### Reptiles
rept_species_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_rept_CAZ", sep = ""), "species_rept_CAZ", runs = 10)
rept_species_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/species_rept_ABF", sep = ""), "species_rept_ABF", runs = 10)
rept_phylo_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_rept_CAZ", sep = ""), "phylo_rept_CAZ", runs = 10)
rept_phylo_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/phylo_rept_ABF", sep = ""), "phylo_rept_ABF", runs = 10)
rept_trait_nodes_present_CAZ <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_rept_CAZ", sep = ""), "trait_nodes_rept_CAZ", runs = 10)
rept_trait_nodes_present_ABF <- extract_zonation_ranks(paste(getwd(), "/output/zonation_runs/present/trait_nodes_rept_ABF", sep = ""), "trait_nodes_rept_ABF", runs = 10)


