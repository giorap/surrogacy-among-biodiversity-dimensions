##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################## DATA PROCESSING (PRESENT) #################
##### -- Set things up -- #####
#### -- Load packages -- ####
my_packages <- c('dplyr', 'magrittr', 'RCurl', ### data-handling
                 'maptools', 'raster', 'sp', 'rgeos', 'rgdal', "PBSmapping", "zonator", ### spatial
                 'nodiv', 'ape', 'geiger', "PVR", "pez", "picante", ### phylogeny
                 'FD', 'geometry', ### traits
                 'pracma', 'Rmisc',  ### algorithms
                 'plot3D', 'fields', 'RColorBrewer' ### graphics
) 
lapply(my_packages, require, character.only = TRUE)

#### -- Source functions from Michael Krabbe Borregaard's R package "nodiv" -- ####
source_github("https://github.com/mkborregaard/nodiv/blob/master/R/Prepare_data.R")
source_github("https://github.com/mkborregaard/nodiv/blob/master/R/Node_based_analysis.R")

##### -- Create reference grid -- ######
#### Generate equal-area projection reference grid at 1 x 1 degree resolution
### Load world boundary
data(wrld_simpl)
### Determine equal-area projection string
cea_crs <- CRS("+proj=cea +lon_0=0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
### Create 1 x 1 degree raster
grid_1 <- raster(xmn = -180, xmx = -18, ymn = -90, ymx = 90, ncol = 162, nrow = 180)
grid_1[] <- 0
### Rasterize world boundary to each grid
## 1 x 1 degree land cover
grid_1 <- rasterize(x = wrld_simpl, y = grid_1, getCover = TRUE)
## Convert americas_grid to equal-area projection
grid_1_cea <- projectRaster(grid_1, crs = cea_crs, over = T)
## Select only cells with at least 20% cover of land
grid_1_cea[] <- ifelse(grid_1_cea[] >= 20, 1, NA)
## Create corresponding data.frame
grid_1_df <- as.data.frame(grid_1_cea, xy = TRUE)
##### -- Set final distribution data extent -- ######
#### Remove Siberia cells
grid_1_df$layer[which(grid_1_df$x < -18592508)] <- NA
#### Remove Antartica cells
grid_1_df$layer[which(grid_1_df$y < -5590615)] <- NA
#### Trim raster
grid_1_cea <- rasterFromXYZ(grid_1_df)

##### -- Load data -- #####
##### Distribution data
#### Tetrapods
distrib_tetrapods_1 <- readRDS("data/distribution/distrib_tetrapods_1.rds")
#### Amphibians
distrib_amph_1 <- readRDS("data/distribution/distrib_amph_1.rds")
#### Birds
distrib_bird_1 <- readRDS("data/distribution/distrib_bird_1.rds")
#### Mammals
distrib_mamm_1 <- readRDS("data/distribution/distrib_mamm_1.rds")
#### Reptiles
distrib_rept_1 <- readRDS("data/distribution/distrib_rept_1.rds")

##### Trait data
#### Tetrapods
traits_tetrapods_complete <- readRDS("data/traits/traits_tetrapods_complete.rds")
#### Amphibians
traits_amph_complete <- readRDS("data/traits/traits_amph_complete.rds")
#### Birds
traits_bird_complete <- readRDS("data/traits/traits_bird_complete.rds")
#### Mammals
traits_mamm_complete <- readRDS("data/traits/traits_mamm_complete.rds")
#### Reptiles
traits_rept_complete <- readRDS("data/traits/traits_rept_complete.rds")
#### Remove species with imputed trait values
### No amphibians imputed
### Birds
## Load data
traits_imputed_bird <- readRDS("data/traits/traits_bird_imputed.rds")
## Filter imputed species out from complete traits data
traits_bird_noImputed <- traits_bird_complete %>% dplyr::filter(!(species %in% traits_imputed_bird$species))
### Mammals
## Load data
traits_imputed_mamm <- readRDS("data/traits/traits_mamm_imputed.rds")
## Filter imputed species out from complete traits data
traits_mamm_noImputed <- traits_mamm_complete %>% dplyr::filter(!(species %in% traits_imputed_mamm$species))
### Reptiles
## Load data
traits_imputed_rept <- readRDS("data/traits/traits_rept_imputed.rds")
## Filter imputed species out from complete traits data
traits_rept_noImputed <- traits_rept_complete %>% dplyr::filter(!(species %in% traits_imputed_rept$species))
#### Remove species with imputed values from tetrapod-level trait dataset
traits_tetrapods_noImputed <- traits_tetrapods_complete %>% dplyr::filter(!(species %in% c(traits_imputed_bird$species, traits_imputed_mamm$species, traits_imputed_rept$species)))

##### Phylogenetic data
#### Tetrapods
tree_tetrapods_complete <- readRDS("data/phylogeny/tree_tetrapods_complete.rds")
#### Amphibians
tree_amph_complete <- readRDS("data/phylogeny/tree_amph_complete.rds")
#### Birds
tree_bird_complete <- readRDS("data/phylogeny/tree_bird_complete.rds")
#### Mammals
tree_mamm_complete <- readRDS("data/phylogeny/tree_mamm_complete.rds")
#### Reptiles
tree_rept_complete <- readRDS("data/phylogeny/tree_rept_complete.rds")

##### Extinction risk data
#### Tetrapods
risk <- readRDS("data/risk/risk.rds")

##### -- Generate phylogenetic cell-by-node matrices -- ######
##### Generate node-by-species matrices using Create_node_by_species_matrix() from Michael Krabbe Borregaard's library(nodiv)
##### Tetrapods
#### Generate node-by-species matrix
node_matrix_tetrapods <- Create_node_by_species_matrix(tree_tetrapods_complete)
### Rotate matrix and turn into data.frame
nodes_tetrapods <- data.frame(species = row.names(t(node_matrix_tetrapods)), t(node_matrix_tetrapods))
### Edit row and column names
row.names(nodes_tetrapods) <- 1:nrow(nodes_tetrapods)
names(nodes_tetrapods) <- gsub("X", "", names(nodes_tetrapods))
##### Amphibians
#### Generate node-by-species matrix
node_matrix_amph <- Create_node_by_species_matrix(tree_amph_complete)
### Rotate matrix and turn into data.frame
nodes_amph <- data.frame(species = row.names(t(node_matrix_amph)), t(node_matrix_amph))
### Edit row and column names
row.names(nodes_amph) <- 1:nrow(nodes_amph)
names(nodes_amph) <- gsub("X", "", names(nodes_amph))
##### Birds
#### Generate node-by-species matrix
node_matrix_bird <- Create_node_by_species_matrix(tree_bird_complete)
### Rotate matrix and turn into data.frame
nodes_bird <- data.frame(species = row.names(t(node_matrix_bird)), t(node_matrix_bird))
### Edit row and column names
row.names(nodes_bird) <- 1:nrow(nodes_bird)
names(nodes_bird) <- gsub("X", "", names(nodes_bird))
##### Mammals
#### Generate node-by-species matrix
node_matrix_mamm <- Create_node_by_species_matrix(tree_mamm_complete)
### Rotate matrix and turn into data.frame
nodes_mamm <- data.frame(species = row.names(t(node_matrix_mamm)), t(node_matrix_mamm))
### Edit row and column names
row.names(nodes_mamm) <- 1:nrow(nodes_mamm)
names(nodes_mamm) <- gsub("X", "", names(nodes_mamm))
##### Reptiles
#### Generate node-by-species matrix
node_matrix_rept <- Create_node_by_species_matrix(tree_rept_complete)
### Rotate matrix and turn into data.frame
nodes_rept <- data.frame(species = row.names(t(node_matrix_rept)), t(node_matrix_rept))
### Edit row and column names
row.names(nodes_rept) <- 1:nrow(nodes_rept)
names(nodes_rept) <- gsub("X", "", names(nodes_rept))

##### Generate phylogenetic cell-by-node matrices
cell_by_node_amph_1 <- matrix_from_nodes(nodes_amph, distrib_amph_1)
cell_by_node_bird_1 <- matrix_from_nodes(nodes_bird, distrib_bird_1)
cell_by_node_mamm_1 <- matrix_from_nodes(nodes_mamm, distrib_mamm_1)
cell_by_node_rept_1 <- matrix_from_nodes(nodes_rept, distrib_rept_1)
cell_by_node_tetrapods_1 <- matrix_from_nodes(nodes_tetrapods, distrib_tetrapods_1)

##### -- Generate functional cell-by-node matrices -- ######
##### Tetrapods
#### Calculate distance matrix
#### NOTE: using dist(method = "euclidean") because tetrapods include a single continuous variable
traits_tetrapods_matrix <- matrix(traits_tetrapods_complete[, -1], dimnames = list(traits_tetrapods_complete$species, names(traits_tetrapods_complete)[-1]))
traits_tetrapods_dist <- dist(traits_tetrapods_matrix, method = "euclidean")
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_tetrapods_tree <- hclust(traits_tetrapods_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_tetrapods_node_matrix <- Create_node_by_species_matrix(traits_tetrapods_tree)
### Rotate matrix and turn into data.frame
traits_tetrapods_nodes <- data.frame(species = row.names(t(traits_tetrapods_node_matrix)), t(traits_tetrapods_node_matrix))
### Edit row and column names
row.names(traits_tetrapods_nodes) <- 1:nrow(traits_tetrapods_nodes)
names(traits_tetrapods_nodes) <- gsub("X", "", names(traits_tetrapods_nodes))
##### Amphibians
#### Calculate distance matrix
#### NOTE: using gowdis() because amphibians include binary/factor variables
traits_amph_dist <- gowdis(traits_amph_complete[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_amph_tree <- hclust(traits_amph_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_amph_node_matrix <- Create_node_by_species_matrix(traits_amph_tree)
### Rotate matrix and turn into data.frame
traits_amph_nodes <- data.frame(species = traits_amph_complete$species, t(traits_amph_node_matrix))
### Edit row and column names
row.names(traits_amph_nodes) <- 1:nrow(traits_amph_nodes)
names(traits_amph_nodes) <- gsub("X", "", names(traits_amph_nodes))
##### Birds
#### Calculate distance matrix
#### NOTE: using gowdis() because birds include binary/factor variables
traits_bird_dist <- gowdis(traits_bird_complete[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_bird_tree <- hclust(traits_bird_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_bird_node_matrix <- Create_node_by_species_matrix(traits_bird_tree)
### Rotate matrix and turn into data.frame
traits_bird_nodes <- data.frame(species = traits_bird_complete$species, t(traits_bird_node_matrix))
### Edit row and column names
row.names(traits_bird_nodes) <- 1:nrow(traits_bird_nodes)
names(traits_bird_nodes) <- gsub("X", "", names(traits_bird_nodes))
##### Mammals
#### Calculate distance matrix
#### NOTE: using gowdis() because mammals include binary/factor variables
traits_mamm_dist <- gowdis(traits_mamm_complete[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_mamm_tree <- hclust(traits_mamm_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_mamm_node_matrix <- Create_node_by_species_matrix(traits_mamm_tree)
### Rotate matrix and turn into data.frame
traits_mamm_nodes <- data.frame(species = traits_mamm_complete$species, t(traits_mamm_node_matrix))
### Edit row and column names
row.names(traits_mamm_nodes) <- 1:nrow(traits_mamm_nodes)
names(traits_mamm_nodes) <- gsub("X", "", names(traits_mamm_nodes))
##### Reptiles
#### Calculate distance matrix
#### NOTE: using gowdis() because reptiles have binary variables
traits_rept_dist <- gowdis(traits_rept_complete[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_rept_tree <- hclust(traits_rept_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_rept_node_matrix <- Create_node_by_species_matrix(traits_rept_tree)
### Rotate matrix and turn into data.frame
traits_rept_nodes <- data.frame(species = traits_rept_complete$species, t(traits_rept_node_matrix))
### Edit row and column names
row.names(traits_rept_nodes) <- 1:nrow(traits_rept_nodes)
names(traits_rept_nodes) <- gsub("X", "", names(traits_rept_nodes))

##### Generate functional cell-by-node matrices
traits_tetrapods_cell_by_node <- matrix_from_nodes(traits_tetrapods_nodes, distrib_tetrapods_1)
traits_amph_cell_by_node <- matrix_from_nodes(traits_amph_nodes, distrib_amph_1)
traits_bird_cell_by_node <- matrix_from_nodes(traits_bird_nodes, distrib_bird_1)
traits_mamm_cell_by_node <- matrix_from_nodes(traits_mamm_nodes, distrib_mamm_1)
traits_rept_cell_by_node <- matrix_from_nodes(traits_rept_nodes, distrib_rept_1)

##### -- Generate species/trait/phylogenetic rasters for input into Zonation software -- ######
##### Tetrapods
### Species
dir.create("data/input_rasters", showWarnings = FALSE)
dir.create("data/input_rasters/present", showWarnings = FALSE)
dir.create("data/input_rasters/present/species", showWarnings = FALSE)
rasters_from_matrix(distrib_tetrapods_1, out_dir = "data/input_rasters/present/species/tetrapods", out_name = "tetrapods_")
#### Phylogenetic nodes
#### Generate rasters from matrices
dir.create("data/input_rasters/present/phylo", showWarnings = FALSE)
rasters_from_matrix(cell_by_node_tetrapods_1, out_dir = "data/input_rasters/present/phylo/tetrapods", out_name = "tetrapods_")
#### Trait nodes
#### Generate rasters from matrices
dir.create("data/input_rasters/present/trait_nodes", showWarnings = FALSE)
rasters_from_matrix(traits_tetrapods_cell_by_node, out_dir = "data/input_rasters/present/trait_nodes/tetrapods", out_name = "tetrapods_")
#### Trait bins
#### Generate rasters from matrices
dir.create("data/input_rasters/present/traits", showWarnings = FALSE)
rasters_from_matrix(cell_by_trait_tetrapods_1, out_dir = "data/input_rasters/present/traits/tetrapods", out_name = "tetrapods_")

##### Amphibians
#### Species
dir.create("data/input_rasters/present/species/amph", showWarnings = FALSE)
rasters_from_matrix(distrib_amph_1, out_dir = "data/input_rasters/present/species/amph", out_name = "amph_")
#### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_amph_1, out_dir = "data/input_rasters/present/phylo/amph", out_name = "amph_")
#### Trait nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_amph_cell_by_node, out_dir = "data/input_rasters/present/trait_nodes/amph", out_name = "amph_")
#### Trait bins
### Generate rasters from matrices
rasters_from_matrix(cell_by_trait_amph_1, out_dir = "data/input_rasters/present/traits/amph", out_name = "amph_")

##### Birds
#### Species
dir.create("data/input_rasters/present/species/bird", showWarnings = FALSE)
rasters_from_matrix(distrib_bird_1, out_dir = "data/input_rasters/present/species/bird", out_name = "bird_")
#### Phylogenetic nodes
#### Generate rasters from matrices
rasters_from_matrix(cell_by_node_bird_1, out_dir = "data/input_rasters/present/phylo/bird", out_name = "bird_")
#### Trait nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_bird_cell_by_node, out_dir = "data/input_rasters/present/trait_nodes/bird", out_name = "bird_")
#### Trait bins
#### Generate rasters from matrices
rasters_from_matrix(cell_by_trait_bird_1, out_dir = "data/input_rasters/present/traits/bird", out_name = "bird_")

##### Mammals
#### Species
dir.create("data/input_rasters/present/species/mamm", showWarnings = FALSE)
rasters_from_matrix(distrib_mamm_1, out_dir = "data/input_rasters/present/species/mamm", out_name = "mamm_")
#### Phylogenetic nodes
#### Generate rasters from matrices
rasters_from_matrix(cell_by_node_mamm_1, out_dir = "data/input_rasters/present/phylo/mamm", out_name = "mamm_")
#### Trait nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_mamm_cell_by_node, out_dir = "data/input_rasters/present/trait_nodes/mamm", out_name = "mamm_")
#### Trait bins
#### Generate rasters from matrices
rasters_from_matrix(cell_by_trait_mamm_1, out_dir = "data/input_rasters/present/traits/mamm", out_name = "mamm_")

##### Reptiles
#### Species
dir.create("data/input_rasters/present/species/rept", showWarnings = FALSE)
rasters_from_matrix(distrib_rept_1, out_dir = "data/input_rasters/present/species/rept", out_name = "rept_")
#### Phylogenetic nodes
#### Generate rasters from matrices
rasters_from_matrix(cell_by_node_rept_1, out_dir = "data/input_rasters/present/phylo/rept", out_name = "rept_")
#### Trait nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_rept_cell_by_node, out_dir = "data/input_rasters/present/trait_nodes/rept", out_name = "rept_")
#### Trait bins
#### Generate rasters from matrices
rasters_from_matrix(cell_by_trait_rept_1, out_dir = "data/input_rasters/present/traits/rept", out_name = "rept_")

##### -- Load rasters for areas of conservation importance -- #####
##### IUCN PAs from IUCN (2017) and Biodiversity Hotspots, Endemic Bird Areas and G200 areas from Brooks et al (2006)
#### IUCN PAs
PA_raster <- readRDS("data/conservation-areas/PA_raster.rds")
PA_df <- as.data.frame(PA_raster, xy = TRUE)
names(PA_df)[3] <- "PA"
#### Biodiversity hotspots
hotspots_raster <- readRDS("data/conservation-areas/hotspots_raster.rds")
hotspots_df <- as.data.frame(hotspots_raster, xy = TRUE)
names(hotspots_df)[3] <- "hotspots"
hotspots_df$hotspots[which(is.na(hotspots_df$hotspots))] <- 0
hotspots_df$hotspots <- hotspots_df$hotspots * 100
##### Endemic bird areas
eba_raster <- readRDS("data/conservation-areas/eba_raster.rds")
eba_df <- as.data.frame(eba_raster, xy = TRUE)
names(eba_df)[3] <- "eba"
eba_df$eba[which(is.na(eba_df$eba))] <- 0
eba_df$eba <- eba_df$eba * 100
#### G200 areas
g200_raster <- readRDS("data/conservation-areas/g200_raster.rds")
g200_df <- as.data.frame(g200_raster, xy = TRUE)
names(g200_df)[3] <- "g200"
g200_df$g200[which(is.na(g200_df$g200))] <- 0
g200_df$g200 <- g200_df$g200 * 100
#### Combine conservation regions data.frames
global_cons_areas <- data.frame(PA_df, hotspots = hotspots_df$hotspots, eba = eba_df$eba, g200 = g200_df$g200)

