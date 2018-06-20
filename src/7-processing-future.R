##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### DATA PROCESSING (FUTURE) ################# 
##### -- FUTURE SCENARIO 1: all DD/NE species not threatened -- #####
##### Identify threatened species
threatened_DDnt <- risk %>% dplyr::filter(risk == "VU" | risk == "EN" | risk == "CR")

##### Isolate thratened species based on each scenario
### Amphibians
species_names_amph_DDnt <- threatened_DDnt %>% dplyr::filter(class == "Amphibians")
species_names_amph_DDnt <- species_names_amph_DDnt$species
### Birds
species_names_bird_DDnt <- threatened_DDnt %>% dplyr::filter(class == "Birds")
species_names_bird_DDnt <- species_names_bird_DDnt$species
### Mammals
species_names_mamm_DDnt <- threatened_DDnt %>% dplyr::filter(class == "Mammals")
species_names_mamm_DDnt <- species_names_mamm_DDnt$species
### Reptiles
species_names_rept_DDnt <- threatened_DDnt %>% dplyr::filter(class == "Reptiles")
species_names_rept_DDnt <- species_names_rept_DDnt$species
### Tetrapods
species_names_tetrapods_DDnt <- threatened_DDnt$species

##### Remove threatened species from data objects
#### Amphibians
### Distrib
distrib_amph_DDnt_1 <- distrib_amph_1[setdiff(names(distrib_amph_1), species_names_amph_DDnt)]
### Traits
traits_amph_DDnt <- traits_amph_complete[!(traits_amph_complete$species %in% species_names_amph_DDnt), ]
### Phylo
tree_amph_DDnt <- drop.tip(tree_amph_complete, species_names_amph_DDnt)
#### Birds
### Distrib
distrib_bird_DDnt_1 <- distrib_bird_1[setdiff(names(distrib_bird_1), species_names_bird_DDnt)]
### Traits
traits_bird_DDnt <- traits_bird_complete[!(traits_bird_complete$species %in% species_names_bird_DDnt), ]
### Phylo
tree_bird_DDnt <- drop.tip(tree_bird_complete, species_names_bird_DDnt)
#### Mammals
### Distrib
distrib_mamm_DDnt_1 <- distrib_mamm_1[setdiff(names(distrib_mamm_1), species_names_mamm_DDnt)]
### Traits
traits_mamm_DDnt <- traits_mamm_complete[!(traits_mamm_complete$species %in% species_names_mamm_DDnt), ]
### Phylo
tree_mamm_DDnt <- drop.tip(tree_mamm_complete, species_names_mamm_DDnt)
#### Reptiles
### Distrib
distrib_rept_DDnt_1 <- distrib_rept_1[setdiff(names(distrib_rept_1), species_names_rept_DDnt)]
### Traits
traits_rept_DDnt <- traits_rept_complete[!(traits_rept_complete$species %in% species_names_rept_DDnt), ]
### Phylo
tree_rept_DDnt <- drop.tip(tree_rept_complete, species_names_rept_DDnt)
#### Tetrapods
### Distrib
distrib_tetrapods_DDnt_1 <- distrib_tetrapods_1[setdiff(names(distrib_tetrapods_1), species_names_tetrapods_DDnt)]
### Traits
traits_tetrapods_DDnt <- traits_tetrapods_complete[!(traits_tetrapods_complete$species %in% species_names_tetrapods_DDnt), ]
### Phylo
tree_tetrapods_DDnt <- drop.tip(tree_tetrapods_complete, species_names_tetrapods_DDnt)

##### -- Generate phylogenetic cell-by-node matrices -- ######
##### Generate node-by-species matrices using Create_node_by_species_matrix() from Michael Krabbe Borregaard's library(nodiv)
##### Tetrapods
#### Generate node-by-species matrix
node_matrix_tetrapods_DDnt <- Create_node_by_species_matrix(tree_tetrapods_DDnt)
### Rotate matrix and turn into data.frame
nodes_tetrapods_DDnt <- data.frame(species = row.names(t(node_matrix_tetrapods_DDnt)), t(node_matrix_tetrapods_DDnt))
### Edit row and column names
row.names(nodes_tetrapods_DDnt) <- 1:nrow(nodes_tetrapods_DDnt)
names(nodes_tetrapods_DDnt) <- gsub("X", "", names(nodes_tetrapods_DDnt))
##### Amphibians
#### Generate node-by-species matrix
node_matrix_amph_DDnt <- Create_node_by_species_matrix(tree_amph_DDnt)
### Rotate matrix and turn into data.frame
nodes_amph_DDnt <- data.frame(species = row.names(t(node_matrix_amph_DDnt)), t(node_matrix_amph_DDnt))
### Edit row and column names
row.names(nodes_amph_DDnt) <- 1:nrow(nodes_amph_DDnt)
names(nodes_amph_DDnt) <- gsub("X", "", names(nodes_amph_DDnt))
##### Birds
#### Generate node-by-species matrix
node_matrix_bird_DDnt <- Create_node_by_species_matrix(tree_bird_DDnt)
### Rotate matrix and turn into data.frame
nodes_bird_DDnt <- data.frame(species = row.names(t(node_matrix_bird_DDnt)), t(node_matrix_bird_DDnt))
### Edit row and column names
row.names(nodes_bird_DDnt) <- 1:nrow(nodes_bird_DDnt)
names(nodes_bird_DDnt) <- gsub("X", "", names(nodes_bird_DDnt))
##### Mammals
#### Generate node-by-species matrix
node_matrix_mamm_DDnt <- Create_node_by_species_matrix(tree_mamm_DDnt)
### Rotate matrix and turn into data.frame
nodes_mamm_DDnt <- data.frame(species = row.names(t(node_matrix_mamm_DDnt)), t(node_matrix_mamm_DDnt))
### Edit row and column names
row.names(nodes_mamm_DDnt) <- 1:nrow(nodes_mamm_DDnt)
names(nodes_mamm_DDnt) <- gsub("X", "", names(nodes_mamm_DDnt))
##### Reptiles
#### Generate node-by-species matrix
node_matrix_rept_DDnt <- Create_node_by_species_matrix(tree_rept_DDnt)
### Rotate matrix and turn into data.frame
nodes_rept_DDnt <- data.frame(species = row.names(t(node_matrix_rept_DDnt)), t(node_matrix_rept_DDnt))
### Edit row and column names
row.names(nodes_rept_DDnt) <- 1:nrow(nodes_rept_DDnt)
names(nodes_rept_DDnt) <- gsub("X", "", names(nodes_rept_DDnt))

##### Generate phylogenetic cell-by-node matrices
cell_by_node_amph_DDnt_1 <- matrix_from_nodes(nodes_amph_DDnt, distrib_amph_DDnt_1)
cell_by_node_bird_DDnt_1 <- matrix_from_nodes(nodes_bird_DDnt, distrib_bird_DDnt_1)
cell_by_node_mamm_DDnt_1 <- matrix_from_nodes(nodes_mamm_DDnt, distrib_mamm_DDnt_1)
cell_by_node_rept_DDnt_1 <- matrix_from_nodes(nodes_rept_DDnt, distrib_rept_DDnt_1)
cell_by_node_tetrapods_DDnt_1 <- matrix_from_nodes(nodes_tetrapods_DDnt, distrib_tetrapods_DDnt_1)

##### -- Generate functional cell-by-node matrices -- ######
##### Tetrapods
#### Calculate distance matrix
#### NOTE: using dist(method = "euclidean") because tetrapods_DDnt include a single continuous variable
traits_tetrapods_DDnt_matrix <- matrix(traits_tetrapods_DDnt[, -1], dimnames = list(traits_tetrapods_DDnt$species, names(traits_tetrapods_DDnt)[-1]))
traits_tetrapods_DDnt_dist <- dist(traits_tetrapods_DDnt_matrix, method = "euclidean")
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_tetrapods_DDnt_tree <- hclust(traits_tetrapods_DDnt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_tetrapods_DDnt_node_matrix <- Create_node_by_species_matrix(traits_tetrapods_DDnt_tree)
### Rotate matrix and turn into data.frame
traits_tetrapods_DDnt_nodes <- data.frame(species = row.names(t(traits_tetrapods_DDnt_node_matrix)), t(traits_tetrapods_DDnt_node_matrix))
### Edit row and column names
row.names(traits_tetrapods_DDnt_nodes) <- 1:nrow(traits_tetrapods_DDnt_nodes)
names(traits_tetrapods_DDnt_nodes) <- gsub("X", "", names(traits_tetrapods_DDnt_nodes))
##### Amphibians
#### Calculate distance matrix
#### NOTE: using gowdis() because amph_DDntibians include binary/factor variables
traits_amph_DDnt_dist <- gowdis(traits_amph_DDnt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_amph_DDnt_tree <- hclust(traits_amph_DDnt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_amph_DDnt_node_matrix <- Create_node_by_species_matrix(traits_amph_DDnt_tree)
### Rotate matrix and turn into data.frame
traits_amph_DDnt_nodes <- data.frame(species = traits_amph_DDnt$species, t(traits_amph_DDnt_node_matrix))
### Edit row and column names
row.names(traits_amph_DDnt_nodes) <- 1:nrow(traits_amph_DDnt_nodes)
names(traits_amph_DDnt_nodes) <- gsub("X", "", names(traits_amph_DDnt_nodes))
##### Birds
#### Calculate distance matrix
#### NOTE: using gowdis() because bird_DDnts include binary/factor variables
traits_bird_DDnt_dist <- gowdis(traits_bird_DDnt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_bird_DDnt_tree <- hclust(traits_bird_DDnt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_bird_DDnt_node_matrix <- Create_node_by_species_matrix(traits_bird_DDnt_tree)
### Rotate matrix and turn into data.frame
traits_bird_DDnt_nodes <- data.frame(species = traits_bird_DDnt$species, t(traits_bird_DDnt_node_matrix))
### Edit row and column names
row.names(traits_bird_DDnt_nodes) <- 1:nrow(traits_bird_DDnt_nodes)
names(traits_bird_DDnt_nodes) <- gsub("X", "", names(traits_bird_DDnt_nodes))
##### Mammals
#### Calculate distance matrix
#### NOTE: using gowdis() because mamm_DDntals include binary/factor variables
traits_mamm_DDnt_dist <- gowdis(traits_mamm_DDnt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_mamm_DDnt_tree <- hclust(traits_mamm_DDnt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_mamm_DDnt_node_matrix <- Create_node_by_species_matrix(traits_mamm_DDnt_tree)
### Rotate matrix and turn into data.frame
traits_mamm_DDnt_nodes <- data.frame(species = traits_mamm_DDnt$species, t(traits_mamm_DDnt_node_matrix))
### Edit row and column names
row.names(traits_mamm_DDnt_nodes) <- 1:nrow(traits_mamm_DDnt_nodes)
names(traits_mamm_DDnt_nodes) <- gsub("X", "", names(traits_mamm_DDnt_nodes))
##### Reptiles
#### Calculate distance matrix
#### NOTE: using gowdis() because rept_DDntiles have binary variables
traits_rept_DDnt_dist <- gowdis(traits_rept_DDnt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_rept_DDnt_tree <- hclust(traits_rept_DDnt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_rept_DDnt_node_matrix <- Create_node_by_species_matrix(traits_rept_DDnt_tree)
### Rotate matrix and turn into data.frame
traits_rept_DDnt_nodes <- data.frame(species = traits_rept_DDnt$species, t(traits_rept_DDnt_node_matrix))
### Edit row and column names
row.names(traits_rept_DDnt_nodes) <- 1:nrow(traits_rept_DDnt_nodes)
names(traits_rept_DDnt_nodes) <- gsub("X", "", names(traits_rept_DDnt_nodes))

##### Generate functional cell-by-node matrices
traits_tetrapods_DDnt_cell_by_node <- matrix_from_nodes(traits_tetrapods_DDnt_nodes, distrib_tetrapods_DDnt_1)
traits_amph_DDnt_cell_by_node <- matrix_from_nodes(traits_amph_DDnt_nodes, distrib_amph_DDnt_1)
traits_bird_DDnt_cell_by_node <- matrix_from_nodes(traits_bird_DDnt_nodes, distrib_bird_DDnt_1)
traits_mamm_DDnt_cell_by_node <- matrix_from_nodes(traits_mamm_DDnt_nodes, distrib_mamm_DDnt_1)
traits_rept_DDnt_cell_by_node <- matrix_from_nodes(traits_rept_DDnt_nodes, distrib_rept_DDnt_1)

##### -- Generate species/trait categories/phylogenetic node rasters for input into Zonation software -- ######
##### Tetrapods
### Species
dir.create("data/input_rasters/DDnt", showWarnings = FALSE)
dir.create("data/input_rasters/DDnt/species", showWarnings = FALSE)
rasters_from_matrix(distrib_tetrapods_DDnt_1, out_dir = "data/input_rasters/DDnt/species/tetrapods", out_name = "tetrapods_")
##### Phylogenetic nodes
#### Generate rasters from matrices
dir.create("data/input_rasters/DDnt/phylo", showWarnings = FALSE)
rasters_from_matrix(cell_by_node_tetrapods_DDnt_1, out_dir = "data/input_rasters/DDnt/phylo/tetrapods", out_name = "tetrapods_")
##### Functional nodes
#### Generate rasters from matrices
dir.create("data/input_rasters/DDnt/trait_nodes", showWarnings = FALSE)
rasters_from_matrix(traits_tetrapods_DDnt_cell_by_node, out_dir = "data/input_rasters/DDnt/trait_nodes/tetrapods", out_name = "tetrapods_")

##### Amphibians
### Species
dir.create("data/input_rasters/DDnt/species/amph", showWarnings = FALSE)
rasters_from_matrix(distrib_amph_DDnt_1, out_dir = "data/input_rasters/DDnt/species/amph", out_name = "amph_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_amph_DDnt_1, out_dir = "data/input_rasters/DDnt/phylo/amph", out_name = "amph_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_amph_DDnt_cell_by_node, out_dir = "data/input_rasters/DDnt/trait_nodes/amph", out_name = "amph_")

##### Birds
### Species
dir.create("data/input_rasters/DDnt/species/bird", showWarnings = FALSE)
rasters_from_matrix(distrib_bird_DDnt_1, out_dir = "data/input_rasters/DDnt/species/bird", out_name = "bird_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_bird_DDnt_1, out_dir = "data/input_rasters/DDnt/phylo/bird", out_name = "bird_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_bird_DDnt_cell_by_node, out_dir = "data/input_rasters/DDnt/trait_nodes/bird", out_name = "bird_")

##### Mammals
### Species
dir.create("data/input_rasters/DDnt/species/mamm", showWarnings = FALSE)
rasters_from_matrix(distrib_mamm_DDnt_1, out_dir = "data/input_rasters/DDnt/species/mamm", out_name = "mamm_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_mamm_DDnt_1, out_dir = "data/input_rasters/DDnt/phylo/mamm", out_name = "mamm_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_mamm_DDnt_cell_by_node, out_dir = "data/input_rasters/DDnt/trait_nodes/mamm", out_name = "mamm_")

##### Reptiles
### Species
dir.create("data/input_rasters/DDnt/species/rept", showWarnings = FALSE)
rasters_from_matrix(distrib_rept_DDnt_1, out_dir = "data/input_rasters/DDnt/species/rept", out_name = "rept_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_rept_DDnt_1, out_dir = "data/input_rasters/DDnt/phylo/rept", out_name = "rept_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_rept_DDnt_cell_by_node, out_dir = "data/input_rasters/DDnt/trait_nodes/rept", out_name = "rept_")

##### -- FUTURE SCENARIO 2: all DD/NE species threatened -- #####
##### Identify threatened species
threatened_DDt <- risk %>% dplyr::filter(risk == "VU" | risk == "EN" | risk == "CR" | risk == "DD" | risk == "NE")

##### Isolate thratened species based on each scenario
#### Amphibians
species_names_amph_DDt <- threatened_DDt %>% dplyr::filter(class == "Amphibians")
species_names_amph_DDt <- species_names_amph_DDt$species
#### Birds
species_names_bird_DDt <- threatened_DDt %>% dplyr::filter(class == "Birds")
species_names_bird_DDt <- species_names_bird_DDt$species
#### Mammals
species_names_mamm_DDt <- threatened_DDt %>% dplyr::filter(class == "Mammals")
species_names_mamm_DDt <- species_names_mamm_DDt$species
#### Reptiles
species_names_rept_DDt <- threatened_DDt %>% dplyr::filter(class == "Reptiles")
species_names_rept_DDt <- species_names_rept_DDt$species
#### Tetrapods
species_names_tetrapods_DDt <- threatened_DDt$species

##### Remove threatened species from data objects
#### Amphibians
### Distrib
distrib_amph_DDt_1 <- distrib_amph_1[setdiff(names(distrib_amph_1), species_names_amph_DDt)]
### Traits
traits_amph_DDt <- traits_amph_complete[!(traits_amph_complete$species %in% species_names_amph_DDt), ]
### Phylo
tree_amph_DDt <- drop.tip(tree_amph_complete, species_names_amph_DDt)
#### Birds
### Distrib
distrib_bird_DDt_1 <- distrib_bird_1[setdiff(names(distrib_bird_1), species_names_bird_DDt)]
### Traits
traits_bird_DDt <- traits_bird_complete[!(traits_bird_complete$species %in% species_names_bird_DDt), ]
### Phylo
tree_bird_DDt <- drop.tip(tree_bird_complete, species_names_bird_DDt)
#### Mammals
### Distrib
distrib_mamm_DDt_1 <- distrib_mamm_1[setdiff(names(distrib_mamm_1), species_names_mamm_DDt)]
### Traits
traits_mamm_DDt <- traits_mamm_complete[!(traits_mamm_complete$species %in% species_names_mamm_DDt), ]
### Phylo
tree_mamm_DDt <- drop.tip(tree_mamm_complete, species_names_mamm_DDt)
#### Reptiles
### Distrib
distrib_rept_DDt_1 <- distrib_rept_1[setdiff(names(distrib_rept_1), species_names_rept_DDt)]
### Traits
traits_rept_DDt <- traits_rept_complete[!(traits_rept_complete$species %in% species_names_rept_DDt), ]
### Phylo
tree_rept_DDt <- drop.tip(tree_rept_complete, species_names_rept_DDt)
#### Tetrapods
### Distrib
distrib_tetrapods_DDt_1 <- distrib_tetrapods_1[setdiff(names(distrib_tetrapods_1), species_names_tetrapods_DDt)]
### Traits
traits_tetrapods_DDt <- traits_tetrapods_complete[!(traits_tetrapods_complete$species %in% species_names_tetrapods_DDt), ]
### Phylo
tree_tetrapods_DDt <- drop.tip(tree_tetrapods_complete, species_names_tetrapods_DDt)

##### -- Generate phylogenetic cell-by-node matrices -- ######
##### Generate node-by-species matrices using Create_node_by_species_matrix() from Michael Krabbe Borregaard's library(nodiv)
##### Tetrapods
#### Generate node-by-species matrix
node_matrix_tetrapods_DDt <- Create_node_by_species_matrix(tree_tetrapods_DDt)
### Rotate matrix and turn into data.frame
nodes_tetrapods_DDt <- data.frame(species = row.names(t(node_matrix_tetrapods_DDt)), t(node_matrix_tetrapods_DDt))
### Edit row and column names
row.names(nodes_tetrapods_DDt) <- 1:nrow(nodes_tetrapods_DDt)
names(nodes_tetrapods_DDt) <- gsub("X", "", names(nodes_tetrapods_DDt))
##### Amphibians
#### Generate node-by-species matrix
node_matrix_amph_DDt <- Create_node_by_species_matrix(tree_amph_DDt)
### Rotate matrix and turn into data.frame
nodes_amph_DDt <- data.frame(species = row.names(t(node_matrix_amph_DDt)), t(node_matrix_amph_DDt))
### Edit row and column names
row.names(nodes_amph_DDt) <- 1:nrow(nodes_amph_DDt)
names(nodes_amph_DDt) <- gsub("X", "", names(nodes_amph_DDt))
##### Birds
#### Generate node-by-species matrix
node_matrix_bird_DDt <- Create_node_by_species_matrix(tree_bird_DDt)
### Rotate matrix and turn into data.frame
nodes_bird_DDt <- data.frame(species = row.names(t(node_matrix_bird_DDt)), t(node_matrix_bird_DDt))
### Edit row and column names
row.names(nodes_bird_DDt) <- 1:nrow(nodes_bird_DDt)
names(nodes_bird_DDt) <- gsub("X", "", names(nodes_bird_DDt))
##### Mammals
#### Generate node-by-species matrix
node_matrix_mamm_DDt <- Create_node_by_species_matrix(tree_mamm_DDt)
### Rotate matrix and turn into data.frame
nodes_mamm_DDt <- data.frame(species = row.names(t(node_matrix_mamm_DDt)), t(node_matrix_mamm_DDt))
### Edit row and column names
row.names(nodes_mamm_DDt) <- 1:nrow(nodes_mamm_DDt)
names(nodes_mamm_DDt) <- gsub("X", "", names(nodes_mamm_DDt))
##### Reptiles
#### Generate node-by-species matrix
node_matrix_rept_DDt <- Create_node_by_species_matrix(tree_rept_DDt)
### Rotate matrix and turn into data.frame
nodes_rept_DDt <- data.frame(species = row.names(t(node_matrix_rept_DDt)), t(node_matrix_rept_DDt))
### Edit row and column names
row.names(nodes_rept_DDt) <- 1:nrow(nodes_rept_DDt)
names(nodes_rept_DDt) <- gsub("X", "", names(nodes_rept_DDt))

##### Generate phylogenetic cell-by-node matrices
cell_by_node_amph_DDt_1 <- matrix_from_nodes(nodes_amph_DDt, distrib_amph_DDt_1)
cell_by_node_bird_DDt_1 <- matrix_from_nodes(nodes_bird_DDt, distrib_bird_DDt_1)
cell_by_node_mamm_DDt_1 <- matrix_from_nodes(nodes_mamm_DDt, distrib_mamm_DDt_1)
cell_by_node_rept_DDt_1 <- matrix_from_nodes(nodes_rept_DDt, distrib_rept_DDt_1)
cell_by_node_tetrapods_DDt_1 <- matrix_from_nodes(nodes_tetrapods_DDt, distrib_tetrapods_DDt_1)

##### -- Generate functional cell-by-node matrices -- ######
##### Tetrapods
#### Calculate distance matrix
#### NOTE: using dist(method = "euclidean") because tetrapods_DDnt include a single continuous variable
traits_tetrapods_DDt_matrix <- matrix(traits_tetrapods_DDt[, -1], dimnames = list(traits_tetrapods_DDt$species, names(traits_tetrapods_DDt)[-1]))
traits_tetrapods_DDt_dist <- dist(traits_tetrapods_DDt_matrix, method = "euclidean")
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_tetrapods_DDt_tree <- hclust(traits_tetrapods_DDt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_tetrapods_DDt_node_matrix <- Create_node_by_species_matrix(traits_tetrapods_DDt_tree)
### Rotate matrix and turn into data.frame
traits_tetrapods_DDt_nodes <- data.frame(species = row.names(t(traits_tetrapods_DDt_node_matrix)), t(traits_tetrapods_DDt_node_matrix))
### Edit row and column names
row.names(traits_tetrapods_DDt_nodes) <- 1:nrow(traits_tetrapods_DDt_nodes)
names(traits_tetrapods_DDt_nodes) <- gsub("X", "", names(traits_tetrapods_DDt_nodes))
##### Amphibians
#### Calculate distance matrix
#### NOTE: using gowdis() because amph_DDtibians include binary/factor variables
traits_amph_DDt_dist <- gowdis(traits_amph_DDt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_amph_DDt_tree <- hclust(traits_amph_DDt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_amph_DDt_node_matrix <- Create_node_by_species_matrix(traits_amph_DDt_tree)
### Rotate matrix and turn into data.frame
traits_amph_DDt_nodes <- data.frame(species = traits_amph_DDt$species, t(traits_amph_DDt_node_matrix))
### Edit row and column names
row.names(traits_amph_DDt_nodes) <- 1:nrow(traits_amph_DDt_nodes)
names(traits_amph_DDt_nodes) <- gsub("X", "", names(traits_amph_DDt_nodes))
##### Birds
#### Calculate distance matrix
#### NOTE: using gowdis() because bird_DDts include binary/factor variables
traits_bird_DDt_dist <- gowdis(traits_bird_DDt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_bird_DDt_tree <- hclust(traits_bird_DDt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_bird_DDt_node_matrix <- Create_node_by_species_matrix(traits_bird_DDt_tree)
### Rotate matrix and turn into data.frame
traits_bird_DDt_nodes <- data.frame(species = traits_bird_DDt$species, t(traits_bird_DDt_node_matrix))
### Edit row and column names
row.names(traits_bird_DDt_nodes) <- 1:nrow(traits_bird_DDt_nodes)
names(traits_bird_DDt_nodes) <- gsub("X", "", names(traits_bird_DDt_nodes))
##### Mammals
#### Calculate distance matrix
#### NOTE: using gowdis() because mamm_DDtals include binary/factor variables
traits_mamm_DDt_dist <- gowdis(traits_mamm_DDt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_mamm_DDt_tree <- hclust(traits_mamm_DDt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_mamm_DDt_node_matrix <- Create_node_by_species_matrix(traits_mamm_DDt_tree)
### Rotate matrix and turn into data.frame
traits_mamm_DDt_nodes <- data.frame(species = traits_mamm_DDt$species, t(traits_mamm_DDt_node_matrix))
### Edit row and column names
row.names(traits_mamm_DDt_nodes) <- 1:nrow(traits_mamm_DDt_nodes)
names(traits_mamm_DDt_nodes) <- gsub("X", "", names(traits_mamm_DDt_nodes))
##### Reptiles
#### Calculate distance matrix
#### NOTE: using gowdis() because rept_DDtiles have binary variables
traits_rept_DDt_dist <- gowdis(traits_rept_DDt[, -1])
#### Calculate dendrogram from distance matrix, then convert dendrogram to phylo object
traits_rept_DDt_tree <- hclust(traits_rept_DDt_dist, method="average") %>% as.phylo()
#### Generate node-by-species matrix from phylo object
traits_rept_DDt_node_matrix <- Create_node_by_species_matrix(traits_rept_DDt_tree)
### Rotate matrix and turn into data.frame
traits_rept_DDt_nodes <- data.frame(species = traits_rept_DDt$species, t(traits_rept_DDt_node_matrix))
### Edit row and column names
row.names(traits_rept_DDt_nodes) <- 1:nrow(traits_rept_DDt_nodes)
names(traits_rept_DDt_nodes) <- gsub("X", "", names(traits_rept_DDt_nodes))

##### Generate functional cell-by-node matrices
traits_tetrapods_DDt_cell_by_node <- matrix_from_nodes(traits_tetrapods_DDt_nodes, distrib_tetrapods_DDt_1)
traits_amph_DDt_cell_by_node <- matrix_from_nodes(traits_amph_DDt_nodes, distrib_amph_DDt_1)
traits_bird_DDt_cell_by_node <- matrix_from_nodes(traits_bird_DDt_nodes, distrib_bird_DDt_1)
traits_mamm_DDt_cell_by_node <- matrix_from_nodes(traits_mamm_DDt_nodes, distrib_mamm_DDt_1)
traits_rept_DDt_cell_by_node <- matrix_from_nodes(traits_rept_DDt_nodes, distrib_rept_DDt_1)

##### -- Generate species/trait categories/phylogenetic node rasters for input into Zonation software -- ######
##### Tetrapods
### Species
dir.create("data/input_rasters/DDt", showWarnings = FALSE)
dir.create("data/input_rasters/DDt/species", showWarnings = FALSE)
rasters_from_matrix(distrib_tetrapods_DDt_1, out_dir = "data/input_rasters/DDt/species/tetrapods", out_name = "tetrapods_")
##### Phylogenetic nodes
#### Generate rasters from matrices
dir.create("data/input_rasters/DDt/phylo", showWarnings = FALSE)
rasters_from_matrix(cell_by_node_tetrapods_DDt_1, out_dir = "data/input_rasters/DDt/phylo/tetrapods", out_name = "tetrapods_")
##### Functional nodes
#### Generate rasters from matrices
dir.create("data/input_rasters/DDt/trait_nodes", showWarnings = FALSE)
rasters_from_matrix(traits_tetrapods_DDt_cell_by_node, out_dir = "data/input_rasters/DDt/trait_nodes/tetrapods", out_name = "tetrapods_")

##### Amphibians
### Species
dir.create("data/input_rasters/DDt/species/amph", showWarnings = FALSE)
rasters_from_matrix(distrib_amph_DDt_1, out_dir = "data/input_rasters/DDt/species/amph", out_name = "amph_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_amph_DDt_1, out_dir = "data/input_rasters/DDt/phylo/amph", out_name = "amph_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_amph_DDt_cell_by_node, out_dir = "data/input_rasters/DDt/trait_nodes/amph", out_name = "amph_")

##### Birds
### Species
dir.create("data/input_rasters/DDt/species/bird", showWarnings = FALSE)
rasters_from_matrix(distrib_bird_DDt_1, out_dir = "data/input_rasters/DDt/species/bird", out_name = "bird_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_bird_DDt_1, out_dir = "data/input_rasters/DDt/phylo/bird", out_name = "bird_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_bird_DDt_cell_by_node, out_dir = "data/input_rasters/DDt/trait_nodes/bird", out_name = "bird_")

##### Mammals
### Species
dir.create("data/input_rasters/DDt/species/mamm", showWarnings = FALSE)
rasters_from_matrix(distrib_mamm_DDt_1, out_dir = "data/input_rasters/DDt/species/mamm", out_name = "mamm_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_mamm_DDt_1, out_dir = "data/input_rasters/DDt/phylo/mamm", out_name = "mamm_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_mamm_DDt_cell_by_node, out_dir = "data/input_rasters/DDt/trait_nodes/mamm", out_name = "mamm_")

##### Reptiles
### Species
dir.create("data/input_rasters/DDt/species/rept", showWarnings = FALSE)
rasters_from_matrix(distrib_rept_DDt_1, out_dir = "data/input_rasters/DDt/species/rept", out_name = "rept_")
### Phylogenetic nodes
## Generate rasters from matrices
rasters_from_matrix(cell_by_node_rept_DDt_1, out_dir = "data/input_rasters/DDt/phylo/rept", out_name = "rept_")
##### Functional nodes
#### Generate rasters from matrices
rasters_from_matrix(traits_rept_DDt_cell_by_node, out_dir = "data/input_rasters/DDt/trait_nodes/rept", out_name = "rept_")
