#########################################################
##### -- Surrogacy among biodiversity dimensions -- #####
#########################################################
#################### DISPLAY ITEMS ######################
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

##### -- FIGURE 3: SPATIAL MISMATCHES VS SURROGACY -- #####
#### Load required objects
distrib_tetrapods_1 <- readRDS("data/distribution/distrib_tetrapods_1.rds")
tetrapods_FD <- readRDS("rds/metrics/tetrapods_FD.rds")
tetrapods_PD <- readRDS("rds/metrics/tetrapods_PD.rds")
#### Merge objects
tetrapods_dimensions <- data.frame(distrib_tetrapods_1[, 1:2], SR = rowSums(distrib_tetrapods_1[, -c(1:2)]))
tetrapods_dimensions <- left_join(tetrapods_dimensions, tetrapods_PD, by = c("x", "y"))
tetrapods_dimensions <- left_join(tetrapods_dimensions, tetrapods_FD, by = c("x", "y"))
#### Generate maps
### Diversity metrics
## Species richness
tiff("figures/tetrapods-map-SR.tif", width = 6.0, height = 3.8, units = 'in', res = 600)
tetrapods_dimensions$SR_hotspots <- ifelse(tetrapods_dimensions$SR >= quantile(tetrapods_dimensions$SR, .83, na.rm = TRUE), 1, 0)
r <- crop(rasterFromXYZ(tetrapods_dimensions[c("x", "y", "SR_hotspots")]), extent(-18592508, -1500508, -5356115, 6918885))
raster::plot(r, legend = FALSE, col = c(grey(.75), "tomato"), axes = FALSE, bigplot = c(0, 1, 0, 1), box = FALSE)
dev.off()
## PD
tiff("figures/tetrapods-map-PD.tif", width = 6.0, height = 3.8, units = 'in', res = 600)
tetrapods_dimensions$PD_hotspots <- ifelse(tetrapods_dimensions$MPD.x >= quantile(tetrapods_dimensions$MPD.x, .83, na.rm = TRUE), 1, 0)
tetrapods_dimensions$PD_hotspots[which(is.na(tetrapods_dimensions$PD_hotspots) | tetrapods_dimensions$SR < 5)] <- 0
tetrapods_dimensions$PD_hotspots[which(tetrapods_dimensions$PD_hotspots == 1 & tetrapods_dimensions$SR_hotspots == 1)] <- 2
r <- crop(rasterFromXYZ(tetrapods_dimensions[c("x", "y", "PD_hotspots")]), extent(-18592508, -1500508, -5356115, 6918885))
raster::plot(r, legend = FALSE, col = c(grey(.75), "deepskyblue3", "orchid4"), axes = FALSE, bigplot = c(0, 1, 0, 1), box = FALSE)
dev.off()
## FD
tiff("figures/tetrapods-map-FD.tif", width = 6.0, height = 3.8, units = 'in', res = 600)
tetrapods_dimensions$FD_hotspots <- ifelse(tetrapods_dimensions$MPD.y >= quantile(tetrapods_dimensions$MPD.y, .83, na.rm = TRUE), 1, 0)
tetrapods_dimensions$FD_hotspots[which(is.na(tetrapods_dimensions$FD_hotspots) | tetrapods_dimensions$SR < 5)] <- 0
tetrapods_dimensions$FD_hotspots[which(tetrapods_dimensions$FD_hotspots == 1 & tetrapods_dimensions$SR_hotspots == 1)] <- 2
r <- crop(rasterFromXYZ(tetrapods_dimensions[c("x", "y", "FD_hotspots")]), extent(-18592508, -1500508, -5356115, 6918885))
raster::plot(r, legend = FALSE, col = c(grey(.75), "deepskyblue3"), axes = FALSE, bigplot = c(0, 1, 0, 1), box = FALSE)
dev.off()
#### Spatial priorities
### Species surrogate
tiff("figures/tetrapods-map-species-priorities.tif", width = 6.0, height = 3.8, units = 'in', res = 600)
r_species <- raster("output/zonation_runs/present/species_tetrapods_CAZ/species_tetrapods_CAZ1/species_tetrapods_CAZ1_out/species_tetrapods_CAZ1.rank.compressed.tif")
r_species[] <- ifelse(r_species[] >= .83, 1, 0)
r_species <- crop(r_species, extent(-18592508, -1500508, -5356115, 6918885))
raster::plot(r_species, legend = FALSE, col = c(grey(.75), "tomato"), axes = FALSE, bigplot = c(0, 1, 0, 1), box = FALSE)
dev.off()
### Phylogenetic target
tiff("figures/tetrapods-map-phylo-priorities.tif", width = 6.0, height = 3.8, units = 'in', res = 600)
r_phylo <- raster("output/zonation_runs/present/phylo_tetrapods_CAZ/phylo_tetrapods_CAZ1/phylo_tetrapods_CAZ1_out/phylo_tetrapods_CAZ1.rank.compressed.tif")
r_phylo[] <- ifelse(r_phylo[] >= .83, 1, 0)
r_phylo[] <- ifelse(r_phylo[] == 1 & r_species[] == 1, 2, r_phylo[])
r <- crop(r_phylo, extent(-18592508, -1500508, -5356115, 6918885))
raster::plot(r_phylo, legend = FALSE, col = c(grey(.75), "deepskyblue3", "orchid4"), axes = FALSE, bigplot = c(0, 1, 0, 1), box = FALSE)
dev.off()
### Functional target
tiff("figures/tetrapods-map-trait-priorities.tif", width = 6.0, height = 3.8, units = 'in', res = 600)
r_trait_nodes <- raster("output/zonation_runs/present/trait_nodes_tetrapods_CAZ/trait_nodes_tetrapods_CAZ1/trait_nodes_tetrapods_CAZ1_out/trait_nodes_tetrapods_CAZ1.rank.compressed.tif")
r_trait_nodes[] <- ifelse(r_trait_nodes[] >= .83, 1, 0)
r_trait_nodes[] <- ifelse(r_trait_nodes[] == 1 & r_species[] == 1, 2, r_trait_nodes[])
r <- crop(r_trait_nodes, extent(-18592508, -1500508, -5356115, 6918885))
raster::plot(r_trait_nodes, legend = FALSE, col = c(grey(.75), "deepskyblue3", "orchid4"), axes = FALSE, bigplot = c(0, 1, 0, 1), box = FALSE)
dev.off()
#### Surrogacy
### Phylogenetic target
tiff("figures/tetrapods-map-phylo-surrogacy.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot((1:nrow(tetrapods_phylo_present_CAZ_curves$optimal)/nrow(tetrapods_phylo_present_CAZ_curves$optimal)) * 100, 
     tetrapods_phylo_present_CAZ_curves$optimal[, 1], 
     type = "n", xlab = "", ylab = "",
     cex.lab = 2, cex.axis = 2, xlim = c(0, 17), ylim = c(50, 100)
     )
lines((1:nrow(tetrapods_phylo_present_CAZ_curves$optimal)/nrow(tetrapods_phylo_present_CAZ_curves$optimal)) * 100, 
      tetrapods_phylo_present_CAZ_curves$optimal[, 1], col = "deepskyblue3", lty = 1, lwd = 4)
lines((1:nrow(tetrapods_phylo_present_CAZ_curves$surrogate)/nrow(tetrapods_phylo_present_CAZ_curves$surrogate)) * 100, 
      tetrapods_phylo_present_CAZ_curves$surrogate[, 1], col = "tomato", lty = 1, lwd = 4)
dev.off()
### Functional target
tiff("figures/tetrapods-map-trait_nodes-surrogacy.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot((1:nrow(tetrapods_trait_nodes_present_CAZ_curves$optimal)/nrow(tetrapods_trait_nodes_present_CAZ_curves$optimal)) * 100, 
     tetrapods_trait_nodes_present_CAZ_curves$optimal[, 1], 
     type = "n", xlab = "", ylab = "",
     cex.lab = 2, cex.axis = 2, xlim = c(0, 17), ylim = c(50, 100)
)
lines((1:nrow(tetrapods_trait_nodes_present_CAZ_curves$optimal)/nrow(tetrapods_trait_nodes_present_CAZ_curves$optimal)) * 100, 
      tetrapods_trait_nodes_present_CAZ_curves$optimal[, 1], col = "deepskyblue3", lty = 1, lwd = 4)
lines((1:nrow(tetrapods_trait_nodes_present_CAZ_curves$surrogate)/nrow(tetrapods_trait_nodes_present_CAZ_curves$surrogate)) * 100, 
      tetrapods_trait_nodes_present_CAZ_curves$surrogate[, 1], col = "tomato", lty = 1, lwd = 4)
dev.off()


##### -- FIGURE 4: DISTINCTIVENESS CHANGE WITH AREA -- #####
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
##### -- SUPPLEMENTARY FIGURE 2 -- #####
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

##### -- SUPPLEMENTARY FIGURE 3 -- #####
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

##### -- SUPPLEMENTARY FIGURE 4 -- #####
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

##### -- SUPPLEMENTARY FIGURE 5 -- #####
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
##### -- TABLE 2 -- #####
#### Read required objects
tetrapods_trait_nodes_DDnt_greedy_curves <- readRDS("output/tetrapods_trait_nodes_DDnt_greedy_curves.rds")
tetrapods_trait_nodes_DDnt_CAZ_curves <- readRDS("output/tetrapods_trait_nodes_DDnt_CAZ_curves.rds")
tetrapods_trait_nodes_DDnt_ABF_curves <- readRDS("output/tetrapods_trait_nodes_DDnt_ABF_curves.rds")
tetrapods_phylo_DDnt_greedy_curves <- readRDS("output/tetrapods_phylo_DDnt_greedy_curves.rds")
tetrapods_phylo_DDnt_CAZ_curves <- readRDS("output/tetrapods_phylo_DDnt_CAZ_curves.rds")
tetrapods_phylo_DDnt_ABF_curves <- readRDS("output/tetrapods_phylo_DDnt_ABF_curves.rds")
tetrapods_trait_nodes_DDt_greedy_curves <- readRDS("output/tetrapods_trait_nodes_DDt_greedy_curves.rds")
tetrapods_trait_nodes_DDt_CAZ_curves <- readRDS("output/tetrapods_trait_nodes_DDt_CAZ_curves.rds")
tetrapods_trait_nodes_DDt_ABF_curves <- readRDS("output/tetrapods_trait_nodes_DDt_ABF_curves.rds")
tetrapods_phylo_DDt_greedy_curves <- readRDS("output/tetrapods_phylo_DDt_greedy_curves.rds")
tetrapods_phylo_DDt_CAZ_curves <- readRDS("output/tetrapods_phylo_DDt_CAZ_curves.rds")
tetrapods_phylo_DDt_ABF_curves <- readRDS("output/tetrapods_phylo_DDt_ABF_curves.rds")
tetrapods_trait_nodes_present_greedy <- readRDS("output/tetrapods_trait_nodes_present_greedy.rds")
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

##### -- TABLE 3 -- #####
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
sai_byTaxon_table <- data.frame(taxon = rep(c("amph", "bird", "mamm", "rept"), each = 3), algorithm = rep(c("greedy", "CAZ", "ABF"), times = 4), phylo_target = NA, trait_nodes_target = NA)
### Phylogenetic target
sai_byTaxon_table$phylo_target <- c(
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
sai_byTaxon_table$trait_nodes_target <- c(
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

##### -- SUPPLEMENTARY TABLE 3 -- #####
#### Tetrapods
#### Load required objects
distrib_tetrapods_1 <- readRDS("data/distribution/distrib_tetrapods_1.rds")
tetrapods_FD <- readRDS("rds/metrics/tetrapods_FD.rds")
tetrapods_PD <- readRDS("rds/metrics/tetrapods_PD.rds")
#### Merge objects
tetrapods_dimensions <- data.frame(distrib_tetrapods_1[, 1:2], SR = rowSums(distrib_tetrapods_1[, -c(1:2)]))
tetrapods_dimensions <- left_join(tetrapods_dimensions, tetrapods_PD, by = c("x", "y"))
tetrapods_dimensions <- left_join(tetrapods_dimensions, tetrapods_FD, by = c("x", "y"))
#### Rename columns
names(tetrapods_dimensions) <- c("x", "y", "SR", "PD.MPD", "PD.PD", "PD.PD.IVS", "FD.MPD", "FD.PD", "FD.PD.IVS", "FD.PSV")
#### Correlations
round(cor.spatial(tetrapods_dimensions$SR, tetrapods_dimensions$PD.MPD, coords = tetrapods_dimensions[c("x", "y")]), 2)
round(cor.spatial(tetrapods_dimensions$SR, tetrapods_dimensions$PD.PD, coords = tetrapods_dimensions[c("x", "y")]), 2)
round(cor.spatial(tetrapods_dimensions$SR, tetrapods_dimensions$FD.MPD, coords = tetrapods_dimensions[c("x", "y")]), 2)
round(cor.spatial(tetrapods_dimensions$SR, tetrapods_dimensions$FD.PD, coords = tetrapods_dimensions[c("x", "y")]), 2)

#### Amphibians
#### Load required objects
distrib_amph_1 <- readRDS("data/distribution/distrib_amph_1.rds")
amph_FD <- readRDS("rds/metrics/amph_FD.rds")
amph_PD <- readRDS("rds/metrics/amph_PD.rds")
#### Merge objects
amph_dimensions <- data.frame(distrib_amph_1[, 1:2], SR = rowSums(distrib_amph_1[, -c(1:2)]))
amph_dimensions <- left_join(amph_dimensions, amph_PD, by = c("x", "y"))
amph_dimensions <- left_join(amph_dimensions, amph_FD, by = c("x", "y"))
#### Rename columns
names(amph_dimensions) <- c("x", "y", "SR", "PD.MPD", "PD.PD", "PD.PD.IVS", "PD.PSV", "FD.MPD", "FD.PD", "FD.PD.IVS", "FD.PSV")
#### Correlations
round(cor.spatial(amph_dimensions$SR, amph_dimensions$PD.MPD, coords = amph_dimensions[c("x", "y")]), 2)
round(cor.spatial(amph_dimensions$SR, amph_dimensions$PD.PD, coords = amph_dimensions[c("x", "y")]), 2)
round(cor.spatial(amph_dimensions$SR, amph_dimensions$FD.MPD, coords = amph_dimensions[c("x", "y")]), 2)
round(cor.spatial(amph_dimensions$SR, amph_dimensions$FD.PD, coords = amph_dimensions[c("x", "y")]), 2)

#### Birds
#### Load required objects
distrib_bird_1 <- readRDS("data/distribution/distrib_bird_1.rds")
bird_FD <- readRDS("rds/metrics/bird_FD.rds")
bird_PD <- readRDS("rds/metrics/bird_PD.rds")
#### Merge objects
bird_dimensions <- data.frame(distrib_bird_1[, 1:2], SR = rowSums(distrib_bird_1[, -c(1:2)]))
bird_dimensions <- left_join(bird_dimensions, bird_PD, by = c("x", "y"))
bird_dimensions <- left_join(bird_dimensions, bird_FD, by = c("x", "y"))
#### Rename columns
names(bird_dimensions) <- c("x", "y", "SR", "PD.MPD", "PD.PD", "PD.PD.IVS", "PD.PSV", "FD.MPD", "FD.PD", "FD.PD.IVS", "FD.PSV")
#### Correlations
round(cor.spatial(bird_dimensions$SR, bird_dimensions$PD.MPD, coords = bird_dimensions[c("x", "y")]), 2)
round(cor.spatial(bird_dimensions$SR, bird_dimensions$PD.PD, coords = bird_dimensions[c("x", "y")]), 2)
round(cor.spatial(bird_dimensions$SR, bird_dimensions$FD.MPD, coords = bird_dimensions[c("x", "y")]), 2)
round(cor.spatial(bird_dimensions$SR, bird_dimensions$FD.PD, coords = bird_dimensions[c("x", "y")]), 2)

#### Mammals
#### Load required objects
distrib_mamm_1 <- readRDS("data/distribution/distrib_mamm_1.rds")
mamm_FD <- readRDS("rds/metrics/mamm_FD.rds")
mamm_PD <- readRDS("rds/metrics/mamm_PD.rds")
#### Merge objects
mamm_dimensions <- data.frame(distrib_mamm_1[, 1:2], SR = rowSums(distrib_mamm_1[, -c(1:2)]))
mamm_dimensions <- left_join(mamm_dimensions, mamm_PD, by = c("x", "y"))
mamm_dimensions <- left_join(mamm_dimensions, mamm_FD, by = c("x", "y"))
#### Rename columns
names(mamm_dimensions) <- c("x", "y", "SR", "PD.MPD", "PD.PD", "PD.PD.IVS", "PD.PSV", "FD.MPD", "FD.PD", "FD.PD.IVS", "FD.PSV")
#### Correlations
round(cor.spatial(mamm_dimensions$SR, mamm_dimensions$PD.MPD, coords = mamm_dimensions[c("x", "y")]), 2)
round(cor.spatial(mamm_dimensions$SR, mamm_dimensions$PD.PD, coords = mamm_dimensions[c("x", "y")]), 2)
round(cor.spatial(mamm_dimensions$SR, mamm_dimensions$FD.MPD, coords = mamm_dimensions[c("x", "y")]), 2)
round(cor.spatial(mamm_dimensions$SR, mamm_dimensions$FD.PD, coords = mamm_dimensions[c("x", "y")]), 2)

#### Reptiles
#### Load required objects
distrib_rept_1 <- readRDS("data/distribution/distrib_rept_1.rds")
rept_FD <- readRDS("rds/metrics/rept_FD.rds")
rept_PD <- readRDS("rds/metrics/rept_PD.rds")
#### Merge objects
rept_dimensions <- data.frame(distrib_rept_1[, 1:2], SR = rowSums(distrib_rept_1[, -c(1:2)]))
rept_dimensions <- left_join(rept_dimensions, rept_PD, by = c("x", "y"))
rept_dimensions <- left_join(rept_dimensions, rept_FD, by = c("x", "y"))
#### Rename columns
names(rept_dimensions) <- c("x", "y", "SR", "PD.MPD", "PD.PD", "PD.PD.IVS", "PD.PSV", "FD.MPD", "FD.PD", "FD.PD.IVS", "FD.PSV")
#### Correlations
round(cor.spatial(rept_dimensions$SR, rept_dimensions$PD.MPD, coords = rept_dimensions[c("x", "y")]), 2)
round(cor.spatial(rept_dimensions$SR, rept_dimensions$PD.PD, coords = rept_dimensions[c("x", "y")]), 2)
round(cor.spatial(rept_dimensions$SR, rept_dimensions$FD.MPD, coords = rept_dimensions[c("x", "y")]), 2)
round(cor.spatial(rept_dimensions$SR, rept_dimensions$FD.PD, coords = rept_dimensions[c("x", "y")]), 2)

##### -- SUPPLEMENTARY TABLE 4 -- #####
#### Using reptiles as an example
### Calculate rank correlations among 10 replicates of each plan
out <- list(cor(rept_species_present_greedy[, 3:ncol(rept_species_present_greedy)], method = "kendall"),
            cor(rept_trait_nodes_present_greedy[, 3:ncol(rept_trait_nodes_present_greedy)], method = "kendall"),
            cor(rept_phylo_present_greedy[, 3:ncol(rept_phylo_present_greedy)], method = "kendall"),
            cor(rept_species_present_CAZ[[1]][, 4:ncol(rept_species_present_CAZ[[1]])], method = "kendall"),
            cor(rept_trait_nodes_present_CAZ[[1]][, 4:ncol(rept_trait_nodes_present_CAZ[[1]])], method = "kendall"),
            cor(rept_phylo_present_CAZ[[1]][, 4:ncol(rept_phylo_present_CAZ[[1]])], method = "kendall"),
            cor(rept_species_present_ABF[[1]][, 4:ncol(rept_species_present_ABF[[1]])], method = "kendall"),
            cor(rept_trait_nodes_present_ABF[[1]][, 4:ncol(rept_trait_nodes_present_ABF[[1]])], method = "kendall"),
            cor(rept_phylo_present_ABF[[1]][, 4:ncol(rept_phylo_present_ABF[[1]])], method = "kendall")
)
## Summarize rank correlations
lapply(out, function(x) {
  x <- x[-c(1, 12, 23, 34, 45, 56, 67, 78, 89, 100)]
  c(as.numeric(quantile(x, .025)), median(x), as.numeric(quantile(x, .975))) %>% round(2)
})    
### Calculate rank correlations between species and phylogenetic/functional target
out <- list(cor(cbind(rept_species_present_greedy[, 3:ncol(rept_species_present_greedy)], rept_trait_nodes_present_greedy[, 3:ncol(rept_trait_nodes_present_greedy)]), method = "kendall"),
            cor(cbind(rept_species_present_greedy[, 3:ncol(rept_species_present_greedy)], rept_phylo_present_greedy[, 3:ncol(rept_trait_nodes_present_greedy)]), method = "kendall"),            
            cor(cbind(rept_species_present_CAZ[[1]][, 4:ncol(rept_species_present_CAZ[[1]])], rept_trait_nodes_present_CAZ[[1]][, 4:ncol(rept_trait_nodes_present_CAZ[[1]])]), method = "kendall"),            
            cor(cbind(rept_species_present_CAZ[[1]][, 4:ncol(rept_species_present_CAZ[[1]])], rept_phylo_present_CAZ[[1]][, 4:ncol(rept_phylo_present_CAZ[[1]])]), method = "kendall"),            
            cor(cbind(rept_species_present_ABF[[1]][, 4:ncol(rept_species_present_ABF[[1]])], rept_trait_nodes_present_ABF[[1]][, 4:ncol(rept_trait_nodes_present_ABF[[1]])]), method = "kendall"),           
            cor(cbind(rept_species_present_ABF[[1]][, 4:ncol(rept_species_present_ABF[[1]])], rept_phylo_present_ABF[[1]][, 4:ncol(rept_phylo_present_ABF[[1]])]), method = "kendall")
)
## Summarize rank correlations
lapply(out, function(x) {
  x <- x[1:10, 11:20]
  c(as.numeric(quantile(x, .025)), median(x), as.numeric(quantile(x, .975))) %>% round(2)
})

##### -- SUPPLEMENTARY TABLE 5 -- #####
sai_protected_threshold_sensitivity_table <- data.frame(taxon = rep(rep(c("WCPA", "hotspots", "eba", "g200"), each = 5), times = 2), target = rep(c("phylo", "trait_nodes"), each = 20), sai_greedy = NA, sai_CAZ = NA, sai_ABF = NA)
sai_protected_threshold_sensitivity_table$sai_greedy <-
  c(calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_greedy_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "greedy", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_greedy_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "greedy", thres = 90)[2]
  )
sai_protected_threshold_sensitivity_table$sai_CAZ <-
  c(calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_CAZ_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "CAZ", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_CAZ_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "CAZ", thres = 90)[2]
  )
sai_protected_threshold_sensitivity_table$sai_ABF <-
  c(calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "WCPA", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_phylo_present_ABF_curves, protected_phylo_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "PA", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "hotspots", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "eba", algorithm = "ABF", thres = 90)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 10)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 30)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 50)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 70)[2],
    calculate_sai_protected(tetrapods_trait_nodes_present_ABF_curves, protected_trait_nodes_tetrapods, priority_area = "g200", algorithm = "ABF", thres = 90)[2]
  )
#### Generate table output
sai_protected_threshold_sensitivity_table$sai_greedy <- round(sai_protected_threshold_sensitivity_table$sai_greedy, 2)
sai_protected_threshold_sensitivity_table$sai_CAZ <- round(sai_protected_threshold_sensitivity_table$sai_CAZ, 2)
sai_protected_threshold_sensitivity_table$sai_ABF <- round(sai_protected_threshold_sensitivity_table$sai_ABF, 2)


