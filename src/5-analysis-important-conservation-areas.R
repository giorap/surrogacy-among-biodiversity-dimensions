##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### DATA ANALYSIS (PRESENT) ##################
##################################################
##### -- Areas of Conservation Importance -- #####
##################################################
##### -- Calculate target biodiversity represented within areas of conservation importance -- #####
##### Estimate percentage area/targets protected
##### Tetrapods
#### Species
#### Current
#### Use a range of thresholds to estimate percentage
protected_species_tetrapods_PA <- protected_targets_species(distrib_tetrapods_1, PA_df)
protected_species_tetrapods_hotspots <- protected_targets_species(distrib_tetrapods_1, hotspots_df)
protected_species_tetrapods_eba <- protected_targets_species(distrib_tetrapods_1, eba_df)
protected_species_tetrapods_g200 <- protected_targets_species(distrib_tetrapods_1, g200_df)
protected_species_tetrapods <- rbind(protected_species_tetrapods_PA, protected_species_tetrapods_hotspots, protected_species_tetrapods_eba, protected_species_tetrapods_g200)
protected_species_tetrapods$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_species_tetrapods, "output/protected_species_tetrapods.rds")
#### Phylogeny
#### Use a range of thresholds to estimate percentage
protected_phylo_tetrapods_PA <- protected_targets_phylo(cell_by_node_tetrapods_1, PA_df, tree = tree_tetrapods_complete)
protected_phylo_tetrapods_hotspots <- protected_targets_phylo(cell_by_node_tetrapods_1, hotspots_df, tree = tree_tetrapods_complete)
protected_phylo_tetrapods_eba <- protected_targets_phylo(cell_by_node_tetrapods_1, eba_df, tree = tree_tetrapods_complete)
protected_phylo_tetrapods_g200 <- protected_targets_phylo(cell_by_node_tetrapods_1, g200_df, tree = tree_tetrapods_complete)
protected_phylo_tetrapods <- rbind(protected_phylo_tetrapods_PA, protected_phylo_tetrapods_hotspots, protected_phylo_tetrapods_eba, protected_phylo_tetrapods_g200)
protected_phylo_tetrapods$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_phylo_tetrapods, "output/protected_phylo_tetrapods.rds")
#### Trait nodes
#### Use a range of thresholds to estimate percentage
protected_trait_nodes_tetrapods_PA <- protected_targets_phylo(traits_tetrapods_cell_by_node, PA_df, branch_lengths = trait_branch_lengths_tetrapods_present)
protected_trait_nodes_tetrapods_hotspots <- protected_targets_phylo(traits_tetrapods_cell_by_node, hotspots_df, branch_lengths = trait_branch_lengths_tetrapods_present)
protected_trait_nodes_tetrapods_eba <- protected_targets_phylo(traits_tetrapods_cell_by_node, eba_df, branch_lengths = trait_branch_lengths_tetrapods_present)
protected_trait_nodes_tetrapods_g200 <- protected_targets_phylo(traits_tetrapods_cell_by_node, g200_df, branch_lengths = trait_branch_lengths_tetrapods_present)
protected_trait_nodes_tetrapods <- rbind(protected_trait_nodes_tetrapods_PA, protected_trait_nodes_tetrapods_hotspots, protected_trait_nodes_tetrapods_eba, protected_trait_nodes_tetrapods_g200)
protected_trait_nodes_tetrapods$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_trait_nodes_tetrapods, "rds/protected/protected_trait_nodes_tetrapods.rds")

##### Amphibians
#### Species
#### Current
#### Use a range of thresholds to estimate percentage
protected_species_amph_PA <- protected_targets_species(distrib_amph_1, PA_df)
protected_species_amph_hotspots <- protected_targets_species(distrib_amph_1, hotspots_df)
protected_species_amph_eba <- protected_targets_species(distrib_amph_1, eba_df)
protected_species_amph_g200 <- protected_targets_species(distrib_amph_1, g200_df)
protected_species_amph <- rbind(protected_species_amph_PA, protected_species_amph_hotspots, protected_species_amph_eba, protected_species_amph_g200)
protected_species_amph$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_species_amph, "output/protected_species_amph.rds")
#### Phylogeny
#### Use a range of thresholds to estimate percentage
protected_phylo_amph_PA <- protected_targets_phylo(cell_by_node_amph_1, PA_df, tree = tree_amph_complete)
protected_phylo_amph_hotspots <- protected_targets_phylo(cell_by_node_amph_1, hotspots_df, tree = tree_amph_complete)
protected_phylo_amph_eba <- protected_targets_phylo(cell_by_node_amph_1, eba_df, tree = tree_amph_complete)
protected_phylo_amph_g200 <- protected_targets_phylo(cell_by_node_amph_1, g200_df, tree = tree_amph_complete)
protected_phylo_amph <- rbind(protected_phylo_amph_PA, protected_phylo_amph_hotspots, protected_phylo_amph_eba, protected_phylo_amph_g200)
protected_phylo_amph$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_phylo_amph, "output/protected_phylo_amph.rds")
#### Trait nodes
#### Use a range of thresholds to estimate percentage
protected_trait_nodes_amph_PA <- protected_targets_phylo(traits_amph_cell_by_node, PA_df, branch_lengths = trait_branch_lengths_amph_present)
protected_trait_nodes_amph_hotspots <- protected_targets_phylo(traits_amph_cell_by_node, hotspots_df, branch_lengths = trait_branch_lengths_amph_present)
protected_trait_nodes_amph_eba <- protected_targets_phylo(traits_amph_cell_by_node, eba_df, branch_lengths = trait_branch_lengths_amph_present)
protected_trait_nodes_amph_g200 <- protected_targets_phylo(traits_amph_cell_by_node, g200_df, branch_lengths = trait_branch_lengths_amph_present)
protected_trait_nodes_amph <- rbind(protected_trait_nodes_amph_PA, protected_trait_nodes_amph_hotspots, protected_trait_nodes_amph_eba, protected_trait_nodes_amph_g200)
protected_trait_nodes_amph$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_trait_nodes_amph, "rds/protected/protected_trait_nodes_amph.rds")

##### Birds
#### Species
#### Current
#### Use a range of thresholds to estimate percentage
protected_species_bird_PA <- protected_targets_species(distrib_bird_1, PA_df)
protected_species_bird_hotspots <- protected_targets_species(distrib_bird_1, hotspots_df)
protected_species_bird_eba <- protected_targets_species(distrib_bird_1, eba_df)
protected_species_bird_g200 <- protected_targets_species(distrib_bird_1, g200_df)
protected_species_bird <- rbind(protected_species_bird_PA, protected_species_bird_hotspots, protected_species_bird_eba, protected_species_bird_g200)
protected_species_bird$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_species_bird, "output/protected_species_bird.rds")
#### Phylogeny
#### Use a range of thresholds to estimate percentage
protected_phylo_bird_PA <- protected_targets_phylo(cell_by_node_bird_1, PA_df, tree = tree_bird_complete)
protected_phylo_bird_hotspots <- protected_targets_phylo(cell_by_node_bird_1, hotspots_df, tree = tree_bird_complete)
protected_phylo_bird_eba <- protected_targets_phylo(cell_by_node_bird_1, eba_df, tree = tree_bird_complete)
protected_phylo_bird_g200 <- protected_targets_phylo(cell_by_node_bird_1, g200_df, tree = tree_bird_complete)
protected_phylo_bird <- rbind(protected_phylo_bird_PA, protected_phylo_bird_hotspots, protected_phylo_bird_eba, protected_phylo_bird_g200)
protected_phylo_bird$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_phylo_bird, "output/protected_phylo_bird.rds")
#### Trait nodes
#### Use a range of thresholds to estimate percentage
protected_trait_nodes_bird_PA <- protected_targets_phylo(traits_bird_cell_by_node, PA_df, branch_lengths = trait_branch_lengths_bird_present)
protected_trait_nodes_bird_hotspots <- protected_targets_phylo(traits_bird_cell_by_node, hotspots_df, branch_lengths = trait_branch_lengths_bird_present)
protected_trait_nodes_bird_eba <- protected_targets_phylo(traits_bird_cell_by_node, eba_df, branch_lengths = trait_branch_lengths_bird_present)
protected_trait_nodes_bird_g200 <- protected_targets_phylo(traits_bird_cell_by_node, g200_df, branch_lengths = trait_branch_lengths_bird_present)
protected_trait_nodes_bird <- rbind(protected_trait_nodes_bird_PA, protected_trait_nodes_bird_hotspots, protected_trait_nodes_bird_eba, protected_trait_nodes_bird_g200)
protected_trait_nodes_bird$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_trait_nodes_bird, "rds/protected/protected_trait_nodes_bird.rds")

##### Mammals
#### Species
#### Current
#### Use a range of thresholds to estimate percentage
protected_species_mamm_PA <- protected_targets_species(distrib_mamm_1, PA_df)
protected_species_mamm_hotspots <- protected_targets_species(distrib_mamm_1, hotspots_df)
protected_species_mamm_eba <- protected_targets_species(distrib_mamm_1, eba_df)
protected_species_mamm_g200 <- protected_targets_species(distrib_mamm_1, g200_df)
protected_species_mamm <- rbind(protected_species_mamm_PA, protected_species_mamm_hotspots, protected_species_mamm_eba, protected_species_mamm_g200)
protected_species_mamm$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_species_mamm, "output/protected_species_mamm.rds")
#### Phylogeny
#### Use a range of thresholds to estimate percentage
protected_phylo_mamm_PA <- protected_targets_phylo(cell_by_node_mamm_1, PA_df, tree = tree_mamm_complete)
protected_phylo_mamm_hotspots <- protected_targets_phylo(cell_by_node_mamm_1, hotspots_df, tree = tree_mamm_complete)
protected_phylo_mamm_eba <- protected_targets_phylo(cell_by_node_mamm_1, eba_df, tree = tree_mamm_complete)
protected_phylo_mamm_g200 <- protected_targets_phylo(cell_by_node_mamm_1, g200_df, tree = tree_mamm_complete)
protected_phylo_mamm <- rbind(protected_phylo_mamm_PA, protected_phylo_mamm_hotspots, protected_phylo_mamm_eba, protected_phylo_mamm_g200)
protected_phylo_mamm$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_phylo_mamm, "output/protected_phylo_mamm.rds")
#### Trait nodes
#### Use a range of thresholds to estimate percentage
protected_trait_nodes_mamm_PA <- protected_targets_phylo(traits_mamm_cell_by_node, PA_df, branch_lengths = trait_branch_lengths_mamm_present)
protected_trait_nodes_mamm_hotspots <- protected_targets_phylo(traits_mamm_cell_by_node, hotspots_df, branch_lengths = trait_branch_lengths_mamm_present)
protected_trait_nodes_mamm_eba <- protected_targets_phylo(traits_mamm_cell_by_node, eba_df, branch_lengths = trait_branch_lengths_mamm_present)
protected_trait_nodes_mamm_g200 <- protected_targets_phylo(traits_mamm_cell_by_node, g200_df, branch_lengths = trait_branch_lengths_mamm_present)
protected_trait_nodes_mamm <- rbind(protected_trait_nodes_mamm_PA, protected_trait_nodes_mamm_hotspots, protected_trait_nodes_mamm_eba, protected_trait_nodes_mamm_g200)
protected_trait_nodes_mamm$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_trait_nodes_mamm, "rds/protected/protected_trait_nodes_mamm.rds")

##### Reptiles
#### Species
#### Current
#### Use a range of thresholds to estimate percentage
protected_species_rept_PA <- protected_targets_species(distrib_rept_1, PA_df)
protected_species_rept_hotspots <- protected_targets_species(distrib_rept_1, hotspots_df)
protected_species_rept_eba <- protected_targets_species(distrib_rept_1, eba_df)
protected_species_rept_g200 <- protected_targets_species(distrib_rept_1, g200_df)
protected_species_rept <- rbind(protected_species_rept_PA, protected_species_rept_hotspots, protected_species_rept_eba, protected_species_rept_g200)
protected_species_rept$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_species_rept, "output/protected_species_rept.rds")
#### Phylogeny
#### Use a range of thresholds to estimate percentage
protected_phylo_rept_PA <- protected_targets_phylo(cell_by_node_rept_1, PA_df, tree = tree_rept_complete)
protected_phylo_rept_hotspots <- protected_targets_phylo(cell_by_node_rept_1, hotspots_df, tree = tree_rept_complete)
protected_phylo_rept_eba <- protected_targets_phylo(cell_by_node_rept_1, eba_df, tree = tree_rept_complete)
protected_phylo_rept_g200 <- protected_targets_phylo(cell_by_node_rept_1, g200_df, tree = tree_rept_complete)
protected_phylo_rept <- rbind(protected_phylo_rept_PA, protected_phylo_rept_hotspots, protected_phylo_rept_eba, protected_phylo_rept_g200)
protected_phylo_rept$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_phylo_rept, "output/protected_phylo_rept.rds")
#### Trait nodes
#### Use a range of thresholds to estimate percentage
protected_trait_nodes_rept_PA <- protected_targets_phylo(traits_rept_cell_by_node, PA_df, branch_lengths = trait_branch_lengths_rept_present)
protected_trait_nodes_rept_hotspots <- protected_targets_phylo(traits_rept_cell_by_node, hotspots_df, branch_lengths = trait_branch_lengths_rept_present)
protected_trait_nodes_rept_eba <- protected_targets_phylo(traits_rept_cell_by_node, eba_df, branch_lengths = trait_branch_lengths_rept_present)
protected_trait_nodes_rept_g200 <- protected_targets_phylo(traits_rept_cell_by_node, g200_df, branch_lengths = trait_branch_lengths_rept_present)
protected_trait_nodes_rept <- rbind(protected_trait_nodes_rept_PA, protected_trait_nodes_rept_hotspots, protected_trait_nodes_rept_eba, protected_trait_nodes_rept_g200)
protected_trait_nodes_rept$priority_areas <- rep(c("PA", "hotspots", "eba", "g200"), each = 11)
saveRDS(protected_trait_nodes_rept, "rds/protected/protected_trait_nodes_rept.rds")