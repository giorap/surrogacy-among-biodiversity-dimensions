#########################################################
##### -- Surrogacy among biodiversity dimensions -- #####
#########################################################
#################### DISPLAY ITEMS ######################
###### Code to generate all display items for the manuscript
###### "Species diversity as a surrogate for the conservation of phylogenetic and functional diversity in terrestrial vertebrates across the Americas" by Rapacciuolo et al.
##### -- Set things up -- #####
##### Create output directory
dir.create("figures")
#########################
##### -- FIGURES -- #####
#########################
##### -- FIGURE 1: SURROGACY ACROSS TETRAPODS -- #####
#### Load required objects
tetrapods_trait_nodes_present_greedy_curves <- readRDS("rds/curves/tetrapods_trait_nodes_present_greedy_curves.rds")
tetrapods_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/tetrapods_trait_nodes_present_CAZ_curves.rds")
tetrapods_trait_nodes_present_ABF_curves <- readRDS("rds/curves/tetrapods_trait_nodes_present_ABF_curves.rds")
tetrapods_phylo_present_greedy_curves <- readRDS("rds/curves/tetrapods_phylo_present_greedy_curves.rds")
tetrapods_phylo_present_CAZ_curves <- readRDS("rds/curves/tetrapods_phylo_present_CAZ_curves.rds")
tetrapods_phylo_present_ABF_curves <- readRDS("rds/curves/tetrapods_phylo_present_ABF_curves.rds")
protected_trait_nodes_tetrapods <- readRDS("rds/protected/protected_trait_nodes_tetrapods.rds")
protected_phylo_tetrapods <- readRDS("rds/protected/protected_phylo_tetrapods.rds")
##### -- Plot target accumulation curves -- #####
#### Panel A
tiff("figures/fig1A-tetrapods-phylo-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(tetrapods_phylo_present_greedy_curves, xlab = "", ylab = "", xaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_tetrapods, threshold == 50)$perc_area_protected, subset(protected_phylo_tetrapods, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel B
tiff("figures/fig1B-tetrapods-phylo-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(tetrapods_phylo_present_CAZ_curves, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_tetrapods, threshold == 50)$perc_area_protected, subset(protected_phylo_tetrapods, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel C
tiff("figures/fig1C-tetrapods-phylo-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(tetrapods_phylo_present_ABF_curves, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_tetrapods, threshold == 50)$perc_area_protected, subset(protected_phylo_tetrapods, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel D
tiff("figures/fig1D-tetrapods-trait_nodes-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(tetrapods_trait_nodes_present_greedy_curves, xlab = "", ylab = "", cex.axis = 1.6)
points(subset(protected_trait_nodes_tetrapods, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_tetrapods, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel E
tiff("figures/fig1E-tetrapods-trait_nodes-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(tetrapods_trait_nodes_present_CAZ_curves, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_tetrapods, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_tetrapods, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel F
tiff("figures/fig1F-tetrapods-trait_nodes-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(tetrapods_trait_nodes_present_ABF_curves, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_tetrapods, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_tetrapods, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()

##### -- FIGURE 2: SURROGACY CHANGE WITH AREA -- #####
### Surrogacy by area
### Phylo
tiff("figures/tetrapods-phylo-SAI-by-area.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_by_area(tetrapods_phylo_present_greedy_curves, tetrapods_phylo_present_CAZ_curves, tetrapods_phylo_present_ABF_curves, plot_vline = 17, xlab = "", ylab = "", cex.axis = 1.6, leg = TRUE)
dev.off()
### trait_nodes
tiff("figures/tetrapods-trait_nodes-SAI-by-area.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_by_area(tetrapods_trait_nodes_present_greedy_curves, tetrapods_trait_nodes_present_CAZ_curves, tetrapods_trait_nodes_present_ABF_curves, plot_vline = 17, xlab = "", ylab = "", cex.axis = 1.6)
dev.off()

##### -- FIGURE 3: DISTINCTIVENESS CHANGE WITH AREA -- #####
### Distinctiveness by area
### Phylo
tiff("figures/tetrapods-phylo-distinct-by-area.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_distinct_by_area(tetrapods_present_evol_distinct, tetrapods_present_evol_distinct_greedy, tetrapods_present_evol_distinct_CAZ, tetrapods_present_evol_distinct_ABF, plot_vline = 17, xlab = "", ylab = "", cex.axis = 1.6, measure = "total_distinct", leg = FALSE)
dev.off()
### trait_nodes
tiff("figures/tetrapods-trait_nodes-distinct-by-area.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_distinct_by_area(tetrapods_present_trait_distinct, tetrapods_present_trait_distinct_greedy, tetrapods_present_trait_distinct_CAZ, tetrapods_present_trait_distinct_ABF, plot_vline = 17, xlab = "", ylab = "", cex.axis = 1.6, measure = "total_distinct")
dev.off()

#######################################
##### -- SUPPLEMENTARY FIGURES -- #####
#######################################
##### -- SUPPLEMENTARY FIGURE 2: SURROGACY IN AMPHIBIANS -- #####
#### Load required objects
amph_trait_nodes_present_greedy_curves <- readRDS("rds/curves/amph_trait_nodes_present_greedy_curves.rds")
amph_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/amph_trait_nodes_present_CAZ_curves.rds")
amph_trait_nodes_present_ABF_curves <- readRDS("rds/curves/amph_trait_nodes_present_ABF_curves.rds")
amph_phylo_present_greedy_curves <- readRDS("rds/curves/amph_phylo_present_greedy_curves.rds")
amph_phylo_present_CAZ_curves <- readRDS("rds/curves/amph_phylo_present_CAZ_curves.rds")
amph_phylo_present_ABF_curves <- readRDS("rds/curves/amph_phylo_present_ABF_curves.rds")
protected_trait_nodes_amph <- readRDS("rds/protected/protected_trait_nodes_amph.rds")
protected_phylo_amph <- readRDS("rds/protected/protected_phylo_amph.rds")
##### -- Plot target accumulation curves -- #####
#### Panel A
tiff("figures/amph-phylo-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(amph_phylo_present_greedy_curves, area = 100, xlab = "", ylab = "", xaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_amph, threshold == 50)$perc_area_protected, subset(protected_phylo_amph, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel B
tiff("figures/amph-phylo-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(amph_phylo_present_CAZ_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_amph, threshold == 50)$perc_area_protected, subset(protected_phylo_amph, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel C
tiff("figures/amph-phylo-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(amph_phylo_present_ABF_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_amph, threshold == 50)$perc_area_protected, subset(protected_phylo_amph, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel D
tiff("figures/amph-trait_nodes-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(amph_trait_nodes_present_greedy_curves, area = 100, xlab = "", ylab = "", cex.axis = 1.6)
points(subset(protected_trait_nodes_amph, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_amph, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel E
tiff("figures/amph-trait_nodes-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(amph_trait_nodes_present_CAZ_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_amph, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_amph, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel F
tiff("figures/amph-trait_nodes-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(amph_trait_nodes_present_ABF_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_amph, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_amph, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()

##### -- SUPPLEMENTARY FIGURE 3: SURROGACY IN BIRDS  -- #####
#### Load required objects
bird_trait_nodes_present_greedy_curves <- readRDS("rds/curves/bird_trait_nodes_present_greedy_curves.rds")
bird_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/bird_trait_nodes_present_CAZ_curves.rds")
bird_trait_nodes_present_ABF_curves <- readRDS("rds/curves/bird_trait_nodes_present_ABF_curves.rds")
bird_phylo_present_greedy_curves <- readRDS("rds/curves/bird_phylo_present_greedy_curves.rds")
bird_phylo_present_CAZ_curves <- readRDS("rds/curves/bird_phylo_present_CAZ_curves.rds")
bird_phylo_present_ABF_curves <- readRDS("rds/curves/bird_phylo_present_ABF_curves.rds")
protected_trait_nodes_bird <- readRDS("rds/protected/protected_trait_nodes_bird.rds")
protected_phylo_bird <- readRDS("rds/protected/protected_phylo_bird.rds")
##### -- Plot target accumulation curves -- #####
#### Panel A
tiff("figures/bird-phylo-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(bird_phylo_present_greedy_curves, area = 100, xlab = "", ylab = "", xaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_bird, threshold == 50)$perc_area_protected, subset(protected_phylo_bird, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel B
tiff("figures/bird-phylo-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(bird_phylo_present_CAZ_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_bird, threshold == 50)$perc_area_protected, subset(protected_phylo_bird, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel C
tiff("figures/bird-phylo-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(bird_phylo_present_ABF_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_bird, threshold == 50)$perc_area_protected, subset(protected_phylo_bird, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel D
tiff("figures/bird-trait_nodes-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(bird_trait_nodes_present_greedy_curves, area = 100, xlab = "", ylab = "", cex.axis = 1.6)
points(subset(protected_trait_nodes_bird, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_bird, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel E
tiff("figures/bird-trait_nodes-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(bird_trait_nodes_present_CAZ_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_bird, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_bird, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel F
tiff("figures/bird-trait_nodes-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(bird_trait_nodes_present_ABF_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_bird, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_bird, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()

##### -- SUPPLEMENTARY FIGURE 4: SURROGACY IN MAMMALS  -- #####
#### Load required objects
mamm_trait_nodes_present_greedy_curves <- readRDS("rds/curves/mamm_trait_nodes_present_greedy_curves.rds")
mamm_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/mamm_trait_nodes_present_CAZ_curves.rds")
mamm_trait_nodes_present_ABF_curves <- readRDS("rds/curves/mamm_trait_nodes_present_ABF_curves.rds")
mamm_phylo_present_greedy_curves <- readRDS("rds/curves/mamm_phylo_present_greedy_curves.rds")
mamm_phylo_present_CAZ_curves <- readRDS("rds/curves/mamm_phylo_present_CAZ_curves.rds")
mamm_phylo_present_ABF_curves <- readRDS("rds/curves/mamm_phylo_present_ABF_curves.rds")
protected_trait_nodes_mamm <- readRDS("rds/protected/protected_trait_nodes_mamm.rds")
protected_phylo_mamm <- readRDS("rds/protected/protected_phylo_mamm.rds")
##### -- Plot target accumulation curves -- #####
#### Panel A
tiff("figures/mamm-phylo-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(mamm_phylo_present_greedy_curves, area = 100, xlab = "", ylab = "", xaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_mamm, threshold == 50)$perc_area_protected, subset(protected_phylo_mamm, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel B
tiff("figures/mamm-phylo-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(mamm_phylo_present_CAZ_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_mamm, threshold == 50)$perc_area_protected, subset(protected_phylo_mamm, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel C
tiff("figures/mamm-phylo-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(mamm_phylo_present_ABF_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_mamm, threshold == 50)$perc_area_protected, subset(protected_phylo_mamm, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel D
tiff("figures/mamm-trait_nodes-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(mamm_trait_nodes_present_greedy_curves, area = 100, xlab = "", ylab = "", cex.axis = 1.6)
points(subset(protected_trait_nodes_mamm, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_mamm, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel E
tiff("figures/mamm-trait_nodes-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(mamm_trait_nodes_present_CAZ_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_mamm, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_mamm, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel F
tiff("figures/mamm-trait_nodes-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(mamm_trait_nodes_present_ABF_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_mamm, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_mamm, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()

##### -- SUPPLEMENTARY FIGURE 5: SURROGACY IN REPTILES -- #####
#### Load required objects
rept_trait_nodes_present_greedy_curves <- readRDS("rds/curves/rept_trait_nodes_present_greedy_curves.rds")
rept_trait_nodes_present_CAZ_curves <- readRDS("rds/curves/rept_trait_nodes_present_CAZ_curves.rds")
rept_trait_nodes_present_ABF_curves <- readRDS("rds/curves/rept_trait_nodes_present_ABF_curves.rds")
rept_phylo_present_greedy_curves <- readRDS("rds/curves/rept_phylo_present_greedy_curves.rds")
rept_phylo_present_CAZ_curves <- readRDS("rds/curves/rept_phylo_present_CAZ_curves.rds")
rept_phylo_present_ABF_curves <- readRDS("rds/curves/rept_phylo_present_ABF_curves.rds")
protected_trait_nodes_rept <- readRDS("rds/protected/protected_trait_nodes_rept.rds")
protected_phylo_rept <- readRDS("rds/protected/protected_phylo_rept.rds")
##### -- Plot target accumulation curves -- #####
#### Panel A
tiff("figures/rept-phylo-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(rept_phylo_present_greedy_curves, area = 100, xlab = "", ylab = "", xaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_rept, threshold == 50)$perc_area_protected, subset(protected_phylo_rept, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel B
tiff("figures/rept-phylo-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(rept_phylo_present_CAZ_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_rept, threshold == 50)$perc_area_protected, subset(protected_phylo_rept, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel C
tiff("figures/rept-phylo-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(rept_phylo_present_ABF_curves, area = 100, xlab = "", ylab = "", xaxt = "n", yaxt = "n", cex.axis = 1.6)
points(subset(protected_phylo_rept, threshold == 50)$perc_area_protected, subset(protected_phylo_rept, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel D
tiff("figures/rept-trait_nodes-present-greedy-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(rept_trait_nodes_present_greedy_curves, area = 100, xlab = "", ylab = "", cex.axis = 1.6)
points(subset(protected_trait_nodes_rept, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_rept, threshold == 50)$perc_target_protected_greedy, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel E
tiff("figures/rept-trait_nodes-present-CAZ-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(rept_trait_nodes_present_CAZ_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_rept, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_rept, threshold == 50)$perc_target_protected_CAZ, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()
#### Panel F
tiff("figures/rept-trait_nodes-present-ABF-SAI.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot_sai_curves(rept_trait_nodes_present_ABF_curves, area = 100, xlab = "", ylab = "", yaxt = "n", cex.axis = 1.6)
points(subset(protected_trait_nodes_rept, threshold == 50)$perc_area_protected, subset(protected_trait_nodes_rept, threshold == 50)$perc_target_protected_ABF, pch = c(0, 16:18), cex = c(1.8, 1.8, 1.8, 2.2))
dev.off()

########################
##### -- TABLES -- #####
########################
##### -- TABLE 1: SURROGACY ACROSS TETRAPODS, NOW AND IN THE FUTURE -- #####
#### Read required objects
tetrapods_trait_nodes_DDnt_greedy_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDnt_greedy_curves.rds")
tetrapods_trait_nodes_DDnt_CAZ_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDnt_CAZ_curves.rds")
tetrapods_trait_nodes_DDnt_ABF_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDnt_ABF_curves.rds")
tetrapods_phylo_DDnt_greedy_curves <- readRDS("rds/curves/tetrapods_phylo_DDnt_greedy_curves.rds")
tetrapods_phylo_DDnt_CAZ_curves <- readRDS("rds/curves/tetrapods_phylo_DDnt_CAZ_curves.rds")
tetrapods_phylo_DDnt_ABF_curves <- readRDS("rds/curves/tetrapods_phylo_DDnt_ABF_curves.rds")
tetrapods_trait_nodes_DDt_greedy_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDt_greedy_curves.rds")
tetrapods_trait_nodes_DDt_CAZ_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDt_CAZ_curves.rds")
tetrapods_trait_nodes_DDt_ABF_curves <- readRDS("rds/curves/tetrapods_trait_nodes_DDt_ABF_curves.rds")
tetrapods_phylo_DDt_greedy_curves <- readRDS("rds/curves/tetrapods_phylo_DDt_greedy_curves.rds")
tetrapods_phylo_DDt_CAZ_curves <- readRDS("rds/curves/tetrapods_phylo_DDt_CAZ_curves.rds")
tetrapods_phylo_DDt_ABF_curves <- readRDS("rds/curves/tetrapods_phylo_DDt_ABF_curves.rds")
#### Calculate SAI values
tetrapods_trait_nodes_DDnt_greedy_sai <- calculate_sai(tetrapods_trait_nodes_DDnt_greedy_curves)
tetrapods_trait_nodes_DDnt_CAZ_sai <- calculate_sai(tetrapods_trait_nodes_DDnt_CAZ_curves)
tetrapods_trait_nodes_DDnt_ABF_sai <- calculate_sai(tetrapods_trait_nodes_DDnt_ABF_curves)
tetrapods_phylo_DDnt_greedy_sai <- calculate_sai(tetrapods_phylo_DDnt_greedy_curves)
tetrapods_phylo_DDnt_CAZ_sai <- calculate_sai(tetrapods_phylo_DDnt_CAZ_curves)
tetrapods_phylo_DDnt_ABF_sai <- calculate_sai(tetrapods_phylo_DDnt_ABF_curves)
tetrapods_trait_nodes_DDt_greedy_sai <- calculate_sai(tetrapods_trait_nodes_DDt_greedy_curves)
tetrapods_trait_nodes_DDt_CAZ_sai <- calculate_sai(tetrapods_trait_nodes_DDt_CAZ_curves)
tetrapods_trait_nodes_DDt_ABF_sai <- calculate_sai(tetrapods_trait_nodes_DDt_ABF_curves)
tetrapods_phylo_DDt_greedy_sai <- calculate_sai(tetrapods_phylo_DDt_greedy_curves)
tetrapods_phylo_DDt_CAZ_sai <- calculate_sai(tetrapods_phylo_DDt_CAZ_curves)
tetrapods_phylo_DDt_ABF_sai <- calculate_sai(tetrapods_phylo_DDt_ABF_curves)

##### -- TABLE 2: SURROGACY FOR INDIVIDUAL TETRAPOD CLASSES -- #####
#### Calculate SAI values
amph_trait_nodes_present_greedy_sai <- calculate_sai(amph_trait_nodes_present_greedy_curves)
amph_trait_nodes_present_CAZ_sai <- calculate_sai(amph_trait_nodes_present_CAZ_curves)
amph_trait_nodes_present_ABF_sai <- calculate_sai(amph_trait_nodes_present_ABF_curves)
amph_phylo_present_greedy_sai <- calculate_sai(amph_phylo_present_greedy_curves)
amph_phylo_present_CAZ_sai <- calculate_sai(amph_phylo_present_CAZ_curves)
amph_phylo_present_ABF_sai <- calculate_sai(amph_phylo_present_ABF_curves)
bird_trait_nodes_present_greedy_sai <- calculate_sai(bird_trait_nodes_present_greedy_curves)
bird_trait_nodes_present_CAZ_sai <- calculate_sai(bird_trait_nodes_present_CAZ_curves)
bird_trait_nodes_present_ABF_sai <- calculate_sai(bird_trait_nodes_present_ABF_curves)
bird_phylo_present_greedy_sai <- calculate_sai(bird_phylo_present_greedy_curves)
bird_phylo_present_CAZ_sai <- calculate_sai(bird_phylo_present_CAZ_curves)
bird_phylo_present_ABF_sai <- calculate_sai(bird_phylo_present_ABF_curves)
mamm_trait_nodes_present_greedy_sai <- calculate_sai(mamm_trait_nodes_present_greedy_curves)
mamm_trait_nodes_present_CAZ_sai <- calculate_sai(mamm_trait_nodes_present_CAZ_curves)
mamm_trait_nodes_present_ABF_sai <- calculate_sai(mamm_trait_nodes_present_ABF_curves)
mamm_phylo_present_greedy_sai <- calculate_sai(mamm_phylo_present_greedy_curves)
mamm_phylo_present_CAZ_sai <- calculate_sai(mamm_phylo_present_CAZ_curves)
mamm_phylo_present_ABF_sai <- calculate_sai(mamm_phylo_present_ABF_curves)
rept_trait_nodes_present_greedy_sai <- calculate_sai(rept_trait_nodes_present_greedy_curves)
rept_trait_nodes_present_CAZ_sai <- calculate_sai(rept_trait_nodes_present_CAZ_curves)
rept_trait_nodes_present_ABF_sai <- calculate_sai(rept_trait_nodes_present_ABF_curves)
rept_phylo_present_greedy_sai <- calculate_sai(rept_phylo_present_greedy_curves)
rept_phylo_present_CAZ_sai <- calculate_sai(rept_phylo_present_CAZ_curves)
rept_phylo_present_ABF_sai <- calculate_sai(rept_phylo_present_ABF_curves)
#### Assemble values
table2 <- data.frame(taxon = rep(c("amph", "bird", "mamm", "rept"), each = 3), algorithm = rep(c("greedy", "CAZ", "ABF"), times = 4), phylo_target = NA, trait_nodes_target = NA)
### Phylogenetic target
table2$phylo_target <- c(
  paste(round(amph_phylo_present_greedy_sai$summary[2], 2), ", (", round(amph_phylo_present_greedy_sai$summary[1], 2), ", ", round(amph_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(amph_phylo_present_CAZ_sai$summary[2], 2), ", (", round(amph_phylo_present_CAZ_sai$summary[1], 2), ", ", round(amph_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(amph_phylo_present_ABF_sai$summary[2], 2), ", (", round(amph_phylo_present_ABF_sai$summary[1], 2), ", ", round(amph_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = ""),
  paste(round(bird_phylo_present_greedy_sai$summary[2], 2), ", (", round(bird_phylo_present_greedy_sai$summary[1], 2), ", ", round(bird_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(bird_phylo_present_CAZ_sai$summary[2], 2), ", (", round(bird_phylo_present_CAZ_sai$summary[1], 2), ", ", round(bird_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(bird_phylo_present_ABF_sai$summary[2], 2), ", (", round(bird_phylo_present_ABF_sai$summary[1], 2), ", ", round(bird_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = ""),
  paste(round(mamm_phylo_present_greedy_sai$summary[2], 2), ", (", round(mamm_phylo_present_greedy_sai$summary[1], 2), ", ", round(mamm_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(mamm_phylo_present_CAZ_sai$summary[2], 2), ", (", round(mamm_phylo_present_CAZ_sai$summary[1], 2), ", ", round(mamm_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(mamm_phylo_present_ABF_sai$summary[2], 2), ", (", round(mamm_phylo_present_ABF_sai$summary[1], 2), ", ", round(mamm_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = ""),
  paste(round(rept_phylo_present_greedy_sai$summary[2], 2), ", (", round(rept_phylo_present_greedy_sai$summary[1], 2), ", ", round(rept_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(rept_phylo_present_CAZ_sai$summary[2], 2), ", (", round(rept_phylo_present_CAZ_sai$summary[1], 2), ", ", round(rept_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(rept_phylo_present_ABF_sai$summary[2], 2), ", (", round(rept_phylo_present_ABF_sai$summary[1], 2), ", ", round(rept_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = "")
)
### Functional target
table2$trait_nodes_target <- c(
  paste(round(amph_trait_nodes_present_greedy_sai$summary[2], 2), ", (", round(amph_trait_nodes_present_greedy_sai$summary[1], 2), ", ", round(amph_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(amph_trait_nodes_present_CAZ_sai$summary[2], 2), ", (", round(amph_trait_nodes_present_CAZ_sai$summary[1], 2), ", ", round(amph_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(amph_trait_nodes_present_ABF_sai$summary[2], 2), ", (", round(amph_trait_nodes_present_ABF_sai$summary[1], 2), ", ", round(amph_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = ""),
  paste(round(bird_trait_nodes_present_greedy_sai$summary[2], 2), ", (", round(bird_trait_nodes_present_greedy_sai$summary[1], 2), ", ", round(bird_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(bird_trait_nodes_present_CAZ_sai$summary[2], 2), ", (", round(bird_trait_nodes_present_CAZ_sai$summary[1], 2), ", ", round(bird_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(bird_trait_nodes_present_ABF_sai$summary[2], 2), ", (", round(bird_trait_nodes_present_ABF_sai$summary[1], 2), ", ", round(bird_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = ""),
  paste(round(mamm_trait_nodes_present_greedy_sai$summary[2], 2), ", (", round(mamm_trait_nodes_present_greedy_sai$summary[1], 2), ", ", round(mamm_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(mamm_trait_nodes_present_CAZ_sai$summary[2], 2), ", (", round(mamm_trait_nodes_present_CAZ_sai$summary[1], 2), ", ", round(mamm_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(mamm_trait_nodes_present_ABF_sai$summary[2], 2), ", (", round(mamm_trait_nodes_present_ABF_sai$summary[1], 2), ", ", round(mamm_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = ""),
  paste(round(rept_trait_nodes_present_greedy_sai$summary[2], 2), ", (", round(rept_trait_nodes_present_greedy_sai$summary[1], 2), ", ", round(rept_trait_nodes_present_greedy_sai$summary[3], 2), ")", sep = ""),
  paste(round(rept_trait_nodes_present_CAZ_sai$summary[2], 2), ", (", round(rept_trait_nodes_present_CAZ_sai$summary[1], 2), ", ", round(rept_trait_nodes_present_CAZ_sai$summary[3], 2), ")", sep = ""),
  paste(round(rept_trait_nodes_present_ABF_sai$summary[2], 2), ", (", round(rept_trait_nodes_present_ABF_sai$summary[1], 2), ", ", round(rept_trait_nodes_present_ABF_sai$summary[3], 2), ")", sep = "")
)

##### -- TABLE 3: SURROGACY OF IMPORTANT CONSERVATION AREAS -- #####
#### Calculate SAI values
tetrapods_phylo_WCPA_greedy_sai <- calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "greedy")
tetrapods_phylo_WCPA_CAZ_sai <- calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "CAZ")
tetrapods_phylo_WCPA_ABF_sai <- calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "ABF")
tetrapods_phylo_eba_greedy_sai <- calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "greedy")
tetrapods_phylo_eba_CAZ_sai <- calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "CAZ")
tetrapods_phylo_eba_ABF_sai <- calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "ABF")
tetrapods_phylo_hotspots_greedy_sai <- calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "greedy")
tetrapods_phylo_hotspots_CAZ_sai <- calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "CAZ")
tetrapods_phylo_hotspots_ABF_sai <- calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "ABF")
tetrapods_phylo_g200_greedy_sai <- calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "greedy")
tetrapods_phylo_g200_CAZ_sai <- calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "CAZ")
tetrapods_phylo_g200_ABF_sai <- calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "ABF")
tetrapods_trait_nodes_WCPA_greedy_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "WCPA", algorithm = "greedy")
tetrapods_trait_nodes_WCPA_CAZ_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "WCPA", algorithm = "CAZ")
tetrapods_trait_nodes_WCPA_ABF_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "WCPA", algorithm = "ABF")
tetrapods_trait_nodes_eba_greedy_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "greedy")
tetrapods_trait_nodes_eba_CAZ_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "CAZ")
tetrapods_trait_nodes_eba_ABF_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "ABF")
tetrapods_trait_nodes_hotspots_greedy_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "greedy")
tetrapods_trait_nodes_hotspots_CAZ_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "CAZ")
tetrapods_trait_nodes_hotspots_ABF_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "ABF")
tetrapods_trait_nodes_g200_greedy_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "greedy")
tetrapods_trait_nodes_g200_CAZ_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "CAZ")
tetrapods_trait_nodes_g200_ABF_sai <- calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "ABF")
