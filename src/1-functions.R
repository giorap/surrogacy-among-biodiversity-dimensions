##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
########################## FUNCTIONS #########################
######################################
##### -- Processing functions -- #####
######################################
##### -- source_github() -- #####
##### Functions for sourcing github functions
source_github <- function(u) {
  # load package
  require(RCurl)
  
  # read script lines from website
  script <- getURL(u, ssl.verifypeer = FALSE)
  
  # parase lines and evaluate in the global environment
  eval(parse(text = script))
}

##### -- Create_node_by_species_matrix() -- #####
##### Create_node_by_species_matrix() function was written by Michael Krabbe Borregaard within the R package "nodiv"
Create_node_by_species_matrix <- function(tree)
{
  # create a matrix with 0s and 1s indicating which species descend from each node
  nodespecies <- matrix(0, nrow = Nnode(tree), ncol = Ntip(tree))
  colnames(nodespecies) <- tree$tip.label
  rownames(nodespecies) <- nodenumbers(tree)
  
  
  for ( i in 1:Nnode(tree))
  {
    nodespecies[i, Node_spec(tree, nodenumbers(tree)[i], names = FALSE)] <- 1
  }
  
  return(nodespecies)
}

##### -- matrix_from_nodes() -- #####
##### Generate site-by-node data.frame from node-by-species variable data.frame
matrix_from_nodes <- function(nodes = NA, distrib_data = NA){
  cell_by_node <- data.frame(distrib_data[, 1:2])
  for(i in 2:ncol(nodes)){
    temp_vector <- nodes[, c(1, i)]
    node_species <- subset(temp_vector, temp_vector[, 2] == 1)$species
    if (length(node_species) > 1){
      cell_by_node$new <- as.numeric(rowSums(distrib_data[, c(which(names(distrib_data) %in% node_species))]) > 0)
      names(cell_by_node)[ncol(cell_by_node)] <- names(temp_vector)[2]
    }
    if (length(node_species) == 1){
      cell_by_node$new <- distrib_data[, node_species]
      names(cell_by_node)[ncol(cell_by_node)] <- names(temp_vector)[2]
    }
  }
  return(cell_by_node)
}

##### -- rasters_from_matrix() -- #####
##### Map categories/nodes aross the study area given a data.frame with ncol >= 3 with coordinates + categories
rasters_from_matrix <- function(out_matrix = NA,
                                crs = cea_crs,
                                out_dir = "output",
                                out_name = NULL){
  dir.create(out_dir, showWarnings = FALSE)
  nvar <- ncol(out_matrix)
  for (i in 3:nvar){
    out_raster <- rasterFromXYZ(out_matrix[, c(1:2, i)], crs = crs)
    writeRaster(out_raster, filename = paste(getwd(), "/", out_dir, "/", out_name, names(out_matrix)[i], sep = ""), format = "GTiff", overwrite = TRUE)
  }
}

##################################################
##### -- Spatial prioritization functions -- #####
##################################################
##### -- Greedy algorithm -- #####
##### Modified from functions written by: Francesco Maria Sabatini. Sapienza, University of Rome
##### published as a supplementary material (Appendix 1) to the paper:  
##### Sabatini F.M. et al. (2016) One taxon does not fit all: herb-layer diversity and stand structural complexity are weak predictors of biodiversity in Fagus sylvatica forests. Ecological Indicators
##### DOI: 10.1016/j.ecolind.2016.04.012
##### -- greedy_algorithm() --#####
#### Sub-routine to calculate the optimal curve, i.e. the curve which accumulates species at the steepest rate
greedy_algorithm <- function(input, nsites = NULL){
  x <- decostand(input, "pa")
  final.order <- which.max(rowSums(x))
  if (is.null(nsites)) nsites <- nrow(x)
  combined <- x[final.order,]
  richness.acc <- sum(combined)
  for(i in 2:nsites){
    y <- decostand(t(t(x) + combined), "pa")
    tmp <- which(rowSums(y) == max(rowSums(y)))  
    max.i <- tmp[which(!tmp %in% final.order)]
    if(length(max.i)==1) 
      final.order <- c(final.order, max.i)
    else
      final.order <- c(final.order, sample(max.i, 1))           
    combined <- y[final.order[i],]
    richness.acc <- c(richness.acc, sum(combined)) 
  }
  return(data.frame(sites=1:nsites, richness.acc, rank=final.order))
}

##### -- greedy_iterations -- #####
greedy_iterations <- function(target_matrix, runs = 10, nsites = 500){
  runs <- runs - 1
  full_greedy_run <- greedy_algorithm(as.matrix(target_matrix[, -c(1:2)]))
  greedy_ranks <- data.frame(full_greedy_run, replicate(runs, full_greedy_run$rank))
  for (run in 1:runs){
    greedy_ranks[1:nsites, run + 3] <- greedy_algorithm(as.matrix(target_matrix[, -c(1:2)]), nsites = nsites)$rank
  }
  names(greedy_ranks)[-c(1:2)] <- paste("rank", 1:(runs + 1), sep = "")
  return(greedy_ranks)
}

##### -- Zonation analyses -- #####
##### This code requires installation of the "zonator" R package by Atte Moilanen, Joona Lehtomaki et al. 
##### The latest CRAN version can be installed using: install.packages("zonator")
##### -- create_zproject_edit() -- #####
##### Modified version of zonator::create_zproject, to call alternative templates files to run either core area Zonation (CAZ) or additive benefit function (ABF)
create_zproject_edit <- function(name, dir, variants, dat_template_file = NULL, algorithm = c("CAZ", "ABF"),
                                 spp_template_file = NULL, spp_template_dir = NULL,
                                 overwrite = FALSE, debug = FALSE, ...) {
  if (!file.exists(dir)) {
    stop("Directory ", dir, " provided does not exist.")
  }
  
  # Create the new location
  project_dir <- file.path(dir, name)
  if (file.exists(project_dir)) {
    if (overwrite) {
      if (debug) message("Removing existing directory ", project_dir)
      unlink(project_dir, recursive = TRUE, force = TRUE)
    } else {
      stop("Project ", project_dir, " already exists and overwrite is off")
    }
  }
  
  if (debug) message("Creating a project directory ", project_dir)
  dir.create(project_dir)
  
  # Create an empty README file for the project
  if (debug) message("Creating an empty README file")
  file.create(file.path(project_dir, "README.md"), showWarnings = FALSE)
  
  # Create the variant subfolders with content
  for (variant in variants) {
    variant_dir <- file.path(project_dir, variant)
    if (debug) message("Creating a variant directory ", variant_dir)
    dir.create(variant_dir)
    
    # If no templates are provided, use the ones shipped with zonator. Change
    # the filenames to match the variant.
    if (is.null(dat_template_file)) {
      dat_template_file <- system.file("extdata", paste("template_", algorithm, ".dat", sep = ""),
                                       package = "zonator")
    }
    dat_to <- file.path(variant_dir, paste0(variant, ".dat"))
    
    # If no templates are provided, use the ones shipped with zonator. Change
    # the filenames to match the variant.
    if (is.null(spp_template_file) & is.null(spp_template_dir)) {
      spp_template_file <- system.file("extdata", "template.spp",
                                       package = "zonator")
    }
    
    # Define the target variant spp file path
    spp_to <- file.path(variant_dir, paste0(variant, ".spp"))
    
    # Copy the templates to the new variant folder
    if (debug) message("Copying template dat-file ", dat_template_file,
                       " to variant directory ", variant_dir)
    if (file.exists(dat_template_file)) {
      file.copy(from = dat_template_file, to = dat_to, overwrite = TRUE)
    } else {
      stop("dat-file template ", dat_template_file, " not found")
    }
    # Work out the details depending if using a template file or a
    # directory of input rasters.
    if (!is.null(spp_template_dir)) {
      # We may have multiple directories
      if (all(sapply(spp_template_dir, function(x) file.exists(x)))) {
        if (debug) {
          if (length(spp_template_dir) > 1) {
            dir_msg <- paste("Creating a spp file from rasters in directories ",
                             paste(spp_template_dir, collapse = ", "))
          } else{
            dir_msg <- paste("Creating a spp file from rasters in directory ",
                             spp_template_dir)
          }
          message(dir_msg)
        }
        create_spp(filename = spp_to, spp_file_dir = spp_template_dir, ...)
      } else {
        stop("Spp template dir ", spp_template_dir, " not found.")
      }
    } else if (!is.null(spp_template_file)) {
      if (file.exists(spp_template_file)) {
        if (debug) {
          message("Copying template spp-file  ", spp_template_file,
                  " to variant directory ", variant_dir)
        }
        file.copy(from = spp_template_file, to = spp_to, overwrite = TRUE)
      } else {
        stop("Input template spp-file ", spp_template_file, " not found!")
      }
    }
    
    # Create to output folder
    output_dir <- file.path(variant_dir, paste0(variant, "_out"))
    if (debug) message("Creating an output directory ", output_dir)
    dir.create(output_dir, recursive = TRUE)
    # Create a bat file, first read the template content
    bat_from <- system.file("extdata", "template.bat", package = "zonator")
    cmd_sequence <- scan(file = bat_from, "character", sep = " ",
                         quiet = TRUE)
    # Replace tokens with actual (relative) paths
    dat_relative <- gsub(paste0(project_dir, .Platform$file.sep), "", dat_to)
    spp_relative <- gsub(paste0(project_dir, .Platform$file.sep), "", spp_to)
    output_dir_relative <- gsub(paste0(project_dir, .Platform$file.sep), "",
                                output_dir)
    cmd_sequence <- gsub("INPUT_DAT", dat_relative, cmd_sequence)
    cmd_sequence <- gsub("INPUT_SPP", spp_relative, cmd_sequence)
    cmd_sequence <- gsub("OUTPUT", file.path(output_dir_relative,
                                             paste0(variant, ".txt")),
                         cmd_sequence)
    # Write bat-file
    bat_to <- file.path(project_dir, paste0(variant, ".bat"))
    if (debug) message("Writing bat file ", bat_to)
    cat(paste0(paste(cmd_sequence, collapse = " "), "\n"), file = bat_to)
  }
  
  return(invisible(NULL))
}

##### -- run_zonation() -- #####
##### Wrapper function for running multiple iterations (select number of runs using "runs" argument) of the Zonation software
##### using species (input = "species"), trait categories (input = "traits") or phylogenetic nodes (input = "phylo") as input rasters
##### and either the CAZ or ABF algorithms ("algorithm" argument).
##### The "period" argument specifies which period ("present", "DDnt" or "DDt") the input rasters should be selected from
##### This wrapper function uses the create_zproject() (the "_edit" version above) and run_bat() functions of the zonator R package
run_zonation <- function(input = c("species", "traits", "trait_nodes", "phylo"), 
                         runs = 100, 
                         group = "tetrapods",
                         algorithm = c("CAZ", "ABF"),
                         period = c("present", "DDnt", "DDt"),
                         working_dir = "C:/Users/Gio/Dropbox/back-up/projects/surrogacy-among-biodiversity-dimensions", 
                         zonation_dir = "output/zonation_runs",
                         ...){
  #### Match input
  input <- match.arg(input)
  algorithm <- match.arg(algorithm)
  period <- match.arg(period)
  #### Set working directory
  setwd(working_dir)
  #### Create directory to store outputs
  dir.create(zonation_dir, showWarnings = FALSE)
  #### Specify ID
  id <- paste(input, group, algorithm, sep = "_")
  #### Define project settings
  zonation_project <- create_zproject_edit(name = id, dir = zonation_dir, 
                                           variants = paste(id, 1:runs, sep = ""),
                                           algorithm = algorithm,
                                           spp_template_dir = paste(working_dir, "/data/input_rasters/", period, "/", input, "/", group, sep = ""),
                                           spp_file_pattern = group,
                                           overwrite = FALSE,
                                           ...)
  #### Run Zonation
  for (run in 1:runs){
    temp_dir <- paste(working_dir, "/output/zonation_runs/", period, "/", id, sep = "")
    setwd(temp_dir)
    run_bat(paste(id, run, ".bat", sep = ""), exe = paste(working_dir, "/data/zig4", sep = ""))
  }
}

##### -- extract_zonation_output() -- #####
##### Function to extract output rasters for spatial prioritizations generated from Zonation runs
extract_zonation_ranks <- function(zonation_dir, variant_name, runs = 100){
  zonation_rasters <- vector("list", runs)
  ### Extract Zonation rasters from n runs
  for (run in 1:runs){
    zonation_rasters[[run]] <- raster(paste(zonation_dir, "/", variant_name, run, "/", variant_name, run, "_out/", variant_name, run, ".rank.compressed.tif", sep = ""))
  }
  ### Extract ranks from
  zonation_ranks <- lapply(zonation_rasters, function(x){
    target_rank_df <- x %>% as.data.frame(xy = TRUE) %>% mutate(cellID = 1:length(x[])) %>% dplyr::filter(complete.cases(.))
    names(target_rank_df)[3] <- "rank"
    target_rank_df <- target_rank_df[order(target_rank_df$rank, decreasing = TRUE), ]
    target_rank_df
  })
  zonation_ranks <- data.frame(zonation_ranks[[1]], do.call("cbind", lapply(zonation_ranks[-1], "[[", 4)))
  names(zonation_ranks)[-c(1:2)] <- c("rank_priority", paste("rank", 1:(runs), sep = ""))
  ### Merge zonation ranks with the full reference coordinates: this will enable assigning the correct rows with the target matrix later on
  zonation_ranks <- list(zonation_ranks = zonation_ranks, reference_coordinates = as.data.frame(zonation_rasters[[1]], xy = TRUE)[, 1:2])
  return(zonation_ranks) 
}

##### -- Random selection of grid cells -- #####
##### -- random_ranks -- #####
##### Generate N (specified in "runs") random sequences of grid cell numbers (based on the total number of grid cells in a target matrix) 
##### and store them within a data.frame
random_ranks <- function(target_ranks, runs = 1000, algorithm = c("greedy", "CAZ", "ABF")){
  algorithm <- match.arg(algorithm)
  if (algorithm == "CAZ" | algorithm == "ABF") target_ranks <- target_ranks$zonation_ranks 
  random_ranks_df <- matrix(NA, nrow = length(target_ranks$rank1), ncol = runs) %>% as.data.frame()
  for(run in 1:runs){
    random_ranks_df[, run] <- c(sample(target_ranks$rank1, length(target_ranks$rank1)))
  }
  names(random_ranks_df) <- paste("rank", 1:(runs), sep = "")
  return(random_ranks_df)
}

####################################
##### -- Surrogacy analyses -- #####
####################################
##### -- Functions for estimating proportion of target diversity represented as a function area -- #####
##### Greedy algorithm
#### The function is the same regardless of whether the target is species, trait categories, or nodes.
target_accumulation_greedy <- function(target_matrix){
  cumulative_target <- apply(target_matrix, 2, cumsum)
  out_curves <- 100 * (rowSums(cumulative_target > 0)/ncol(target_matrix))  
  return(out_curves)
}
##### CAZ
##### Equations derived from Pollock et al. (2017) Nature (doi: 10.1038/nature22368)
#### Trait categories
#### Note this is the same function as that used for the greedy algorithm
traits_CAZ <- target_accumulation_greedy
#### Phylogenetic nodes
phylo_CAZ <- function(target_matrix, branches){
  cumulative_target <- apply(target_matrix, 2, cumsum)
  out_curves <- 100 * (rowSums(t(t(apply(cumulative_target > 0, 2, as.numeric)) * branches))/sum(branches))
  return(out_curves)
}
##### ABF
##### Equations derived from Pollock et al. (2017) Nature (doi: 10.1038/nature22368)
#### Trait categories
traits_ABF <- function(target_matrix, ranges){
  cumulative_target <- apply(target_matrix, 2, cumsum)
  out_curves <- (1/ncol(target_matrix)) * rowSums(t(t(cumulative_target)/ranges)) * 100
  return(out_curves)
}
#### Phylogenetic nodes
phylo_ABF <- function(target_matrix, ranges, branches){
  cumulative_target <- apply(target_matrix, 2, cumsum)
  transformed_target <- t(t(apply(cumulative_target > 0, 2, as.numeric)) * branches) * t(t(cumulative_target)/ranges)
  out_curves <- (1/sum(branches)) * rowSums(transformed_target) * 100
  return(out_curves)
}

##### -- Species Accumulation Index -- #####
##### -- get_SAI_curves -- #####
##### Derive accumulation curves based on spatial prioritization ranks and the appropriate target diversity function
get_sai_curves <- function(target_matrix,
                           target_ranks,
                           surrogate_ranks,
                           random_ranks,
                           target = c("traits", "phylo"), 
                           algorithm = c("greedy", "CAZ", "ABF"),
                           runs = 10,
                           branch_lengths = NULL){
  #### Match arguments
  target <- match.arg(target)
  algorithm <- match.arg(algorithm)
  #### Calculate range sizes
  range_sizes <- colSums(target_matrix[, -c(1:2)])
  #### For algorithms CAZ or ABF, merge the target matrix back to the extent of full raster extent
  if (algorithm == "CAZ" | algorithm == "ABF"){
    #### Update objects
    reference_coordinates <- target_ranks$reference_coordinates
    target_ranks <- target_ranks$zonation_ranks
    surrogate_ranks <- surrogate_ranks$zonation_ranks
    target_matrix <- left_join(reference_coordinates, target_matrix, by = c("x", "y"))    
  }
  #### Subset colums and rows from rank data.frames
  target_ranks <- target_ranks[paste("rank", 1:runs, sep = "")] 
  surrogate_ranks <- surrogate_ranks[paste("rank", 1:runs, sep = "")] 
  #### Estimate curves
  ### Optimal
  ## Create output object
  O_curves <- data.frame(matrix(NA, nrow = nrow(target_ranks), ncol = ncol(target_ranks)))
  names(O_curves) <- paste("curve", 1:ncol(O_curves), sep = "")
  ## Calculate accumulation of appropriate target and algorithm combination
  for (i in 1:ncol(O_curves)){
    if (target == "traits" & algorithm == "greedy") O_curves[, i] <- target_accumulation_greedy(target_matrix[target_ranks[, i], -c(1:2)])
    if (target == "traits" & algorithm == "CAZ") O_curves[, i] <- traits_CAZ(target_matrix[target_ranks[, i], -c(1:2)])
    if (target == "traits" & algorithm == "ABF") O_curves[, i] <- traits_ABF(target_matrix[target_ranks[, i], -c(1:2)], ranges = range_sizes)
    if (target == "phylo" & algorithm == "greedy") O_curves[, i] <- target_accumulation_greedy(target_matrix[target_ranks[, i], -c(1:2)])
    if (target == "phylo" & algorithm == "CAZ") O_curves[, i] <- phylo_CAZ(target_matrix[target_ranks[, i], -c(1:2)], branches = branch_lengths)
    if (target == "phylo" & algorithm == "ABF") O_curves[, i] <- phylo_ABF(target_matrix[target_ranks[, i], -c(1:2)], ranges = range_sizes, branches = branch_lengths)
  }
  # Filter curves
  O_curves <- O_curves %>% dplyr::filter(complete.cases(.))
  ### Surrogate
  ## Create output object
  S_curves <- data.frame(matrix(NA, nrow = nrow(surrogate_ranks), ncol = ncol(surrogate_ranks)))
  names(S_curves) <- paste("curve", 1:ncol(S_curves), sep = "")
  ## Calculate accumulation of appropriate target and algorithm combination
  for (i in 1:ncol(S_curves)){
    if (target == "traits" & algorithm == "greedy") S_curves[, i] <- target_accumulation_greedy(target_matrix[surrogate_ranks[, i], -c(1:2)])
    if (target == "traits" & algorithm == "CAZ") S_curves[, i] <- traits_CAZ(target_matrix[surrogate_ranks[, i], -c(1:2)])
    if (target == "traits" & algorithm == "ABF") S_curves[, i] <- traits_ABF(target_matrix[surrogate_ranks[, i], -c(1:2)], ranges = range_sizes)
    if (target == "phylo" & algorithm == "greedy") S_curves[, i] <- target_accumulation_greedy(target_matrix[surrogate_ranks[, i], -c(1:2)])
    if (target == "phylo" & algorithm == "CAZ") S_curves[, i] <- phylo_CAZ(target_matrix[surrogate_ranks[, i], -c(1:2)], branches = branch_lengths)
    if (target == "phylo" & algorithm == "ABF") S_curves[, i] <- phylo_ABF(target_matrix[surrogate_ranks[, i], -c(1:2)], ranges = range_sizes, branches = branch_lengths)
  }
  # Filter curves
  S_curves <- S_curves %>% dplyr::filter(complete.cases(.))
  ### Random
  ## Create output object
  R_curves <- random_ranks
  names(R_curves) <- paste("curve", 1:ncol(R_curves), sep = "")
  for (i in 1:ncol(R_curves)){
    if (target == "traits" & algorithm == "greedy") R_curves[, i] <- target_accumulation_greedy(target_matrix[random_ranks[, i], -c(1:2)])
    if (target == "traits" & algorithm == "CAZ") R_curves[, i] <- traits_CAZ(target_matrix[random_ranks[, i], -c(1:2)])
    if (target == "traits" & algorithm == "ABF") R_curves[, i] <- traits_ABF(target_matrix[random_ranks[, i], -c(1:2)], ranges = range_sizes)
    if (target == "phylo" & algorithm == "greedy") R_curves[, i] <- target_accumulation_greedy(target_matrix[random_ranks[, i], -c(1:2)])
    if (target == "phylo" & algorithm == "CAZ") R_curves[, i] <- phylo_CAZ(target_matrix[random_ranks[, i], -c(1:2)], branches = branch_lengths)
    if (target == "phylo" & algorithm == "ABF") R_curves[, i] <- phylo_ABF(target_matrix[random_ranks[, i], -c(1:2)], ranges = range_sizes, branches = branch_lengths)
  }
  # Filter curves
  R_curves <- R_curves %>% dplyr::filter(complete.cases(.))
  ### Merge output
  curves <- list(optimal = O_curves, surrogate = S_curves, random = R_curves)
  ### Return output
  return(curves)
}

##### -- plot_sai_curves() -- #####
##### Plot derived accumulation curves 
plot_sai_curves <- function(sai_curves, area = 50, ...){
  
  ### Summarize Optimal curves
  O_mean <- rowMeans(sai_curves$optimal)
  O_CI_lower <- as.data.frame(apply(sai_curves$optimal, 1, function(x){ CI(as.vector(x))[3] }))[, 1]
  O_CI_upper <- as.data.frame(apply(sai_curves$optimal, 1, function(x){ CI(as.vector(x))[1] }))[, 1]
  ## Add an initial 0 to start curves at 0%
  O_mean <- c(0, O_mean); O_CI_lower <- c(0, O_CI_lower); O_CI_upper <- c(0, O_CI_upper)
  
  ### Summarize Surrogate curves
  S_mean <- rowMeans(sai_curves$surrogate)
  S_CI_lower <- as.data.frame(apply(sai_curves$surrogate, 1, function(x){ CI(as.vector(x))[3] }))[, 1]
  S_CI_upper <- as.data.frame(apply(sai_curves$surrogate, 1, function(x){ CI(as.vector(x))[1] }))[, 1]
  ## Add an initial 0 to start curves at 0%
  S_mean <- c(0, S_mean); S_CI_lower <- c(0, S_CI_lower); S_CI_upper <- c(0, S_CI_upper)
  
  ### Summarize Optimal curves
  R_mean <- rowMeans(sai_curves$random)
  R_CI_lower <- as.data.frame(apply(sai_curves$random, 1, function(x){ CI(as.vector(x))[3] }))[, 1]
  R_CI_upper <- as.data.frame(apply(sai_curves$random, 1, function(x){ CI(as.vector(x))[1] }))[, 1]
  ## Add an initial 0 to start curves at 0%
  R_mean <- c(0, R_mean); R_CI_lower <- c(0, R_CI_lower); R_CI_upper <- c(0, R_CI_upper)
  
  ### Generate x axis
  x_axis <- (0:nrow(sai_curves$optimal)/nrow(sai_curves$optimal)) * 100
  
  #### Generate plot
  plot(x_axis, R_CI_lower, type = "n", xlim = c(0, area), ...)
  ### Plot Random curve
  polygon(c(x_axis, rev(x_axis)), c(R_CI_lower, rev(R_CI_upper)), col = alpha("green", 0.3), border=alpha("green", 0.5))
  lines(x_axis, R_mean, col = "green", lwd = 2)
  #lines(x_axis, R_CI_lower, col = "green", lwd = 1)
  #lines(x_axis, R_CI_upper, col = "green", lwd = 1) 
  ### Plot Optimal curve
  polygon(c(x_axis, rev(x_axis)), c(O_CI_lower, rev(O_CI_upper)), col = alpha("deepskyblue3", 0.3), border=alpha("deepskyblue3", 0.5))
  lines(x_axis, O_mean, col = "deepskyblue3", lwd = 2)
  #lines(x_axis, O_CI_lower, col = "deepskyblue3", lwd = 1)
  #lines(x_axis, O_CI_upper, col = "deepskyblue3", lwd = 1)
  ### Plot Surrogate curve
  polygon(c(x_axis, rev(x_axis)), c(S_CI_lower, rev(S_CI_upper)), col = alpha("tomato", 0.3), border=alpha("tomato", 0.5))
  lines(x_axis, S_mean, col = "tomato", lwd = 2)
  #lines(x_axis, S_CI_lower, col = "tomato", lwd = 1)
  #lines(x_axis, S_CI_upper, col = "tomato", lwd = 1)
}

##### -- calculate_sai() -- #####
##### Calculate the Species Accumulation Index from the optimal/target, surrogate and random curves 
calculate_sai <- function(sai_curves){
  #### Create a percentage area object
  area <- (0:nrow(sai_curves$optimal)/nrow(sai_curves$optimal)) * 100
  #### Calculate area under curves
  ### Create output objects
  sai_areas_optimal <- sai_areas_surrogate <- sai_areas_random <- sai_values <- NULL
  ### Area under Optimal curves
  for (i in 1:ncol(sai_curves$optimal)){
    sai_areas_optimal[i] <- trapz(area, c(0, sai_curves$optimal[, i]))
  }
  ### Area under Surrogate curves
  for (i in 1:ncol(sai_curves$surrogate)){
    sai_areas_surrogate[i] <- trapz(area, c(0, sai_curves$surrogate[, i]))
  }
  ### Area under Random curves
  for (i in 1:ncol(sai_curves$random)){
    sai_areas_random[i] <- trapz(area, c(0, sai_curves$random[, i]))
  }
  sai_values <- as.vector(outer(rep(sai_areas_surrogate, times = length(sai_areas_optimal)), sai_areas_random, "-"))/as.vector(outer(rep(sai_areas_optimal, each = length(sai_areas_surrogate)), sai_areas_random, "-"))
  sai_values <- sai_values[sai_values != "-Inf" & sai_values != "NaN"]
  sai_summary <- c(quantile(sai_values, .025, na.rm = TRUE), quantile(sai_values, .5, na.rm = TRUE), quantile(sai_values, .975, na.rm = TRUE))
  names(sai_summary) <- c("CI_lower", "median", "CI_upper")
  return(list(summary = sai_summary, all_values = sai_values))
}

##### -- calculate_sai_by_area() -- #####
#### Calculate surrogacy values at fixed 1% intervals up to a maximum
calculate_sai_by_area <- function(sai_curves, max_area = 50){
  #### Create a percentage area object
  area <- (0:nrow(sai_curves$optimal)/nrow(sai_curves$optimal)) * 100
  sai_rows <- NULL
  for (i in 1:max_area) sai_rows[i] <- which.min(abs(area - i))[1]
  sai_curves <- lapply(sai_curves, function(x) x[sai_rows, ])
  target_max <- max(sai_curves$optimal)
  change_by_area <- data.frame(area = 1:max_area)
  #### Calculate SAI values by area
  ### Create output object
  sai_values <- vector("list", nrow(sai_curves$optimal))
  for (i in 1:nrow(sai_curves$optimal)){
    sai_areas_optimal <- sai_curves$optimal[i, ] %>% as.numeric()
    sai_areas_surrogate <- sai_curves$surrogate[i, ] %>% as.numeric()
    sai_areas_random <- sai_curves$random[i, ] %>% as.numeric()
    sai_values[[i]] <- as.vector(outer(rep(sai_areas_surrogate, times = length(sai_areas_optimal)), sai_areas_random, "-"))/as.vector(outer(rep(sai_areas_optimal, each = length(sai_areas_surrogate)), sai_areas_random, "-"))
    sai_values[[i]] <- ifelse(is.na(sai_values[[i]]), 1, sai_values[[i]])
    sai_values[[i]] <- ifelse(sai_values[[i]] > 1, 1, sai_values[[i]])
  }
  sai_summary <- lapply(sai_values, function(x) {
    c(quantile(x, .025, na.rm = TRUE), quantile(x, .5, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
  }
  )
  sai_summary <- data.frame(change_by_area, do.call("rbind", sai_summary))
  names(sai_summary) <- c("area", "SAI_lower", "SAI_median", "SAI_upper")
  return(sai_summary)
}

calculate_sai_protected <- function(sai_curves, protected_target, thres = 50, priority_area = "WCPA", algorithm = "greedy"){
  #### Create a percentage area object
  area <- (0:nrow(sai_curves$optimal)/nrow(sai_curves$optimal)) * 100
  #### Subset protected_target
  protected_target <- subset(protected_target, threshold == thres & priority_areas == priority_area)
  protected_field <- paste("perc_target_protected_", algorithm, sep = "")
  target_row <- which.min(abs(area - protected_target$perc_area_protected))
  #### Calculate SAI values
  ### Create output object
  sai_values <- NULL
  sai_areas_optimal <- sai_curves$surrogate[target_row, ] %>% as.numeric()
  sai_areas_surrogate <- protected_target[, protected_field] %>% as.numeric() 
  sai_areas_random <- sai_curves$random[target_row, ] %>% as.numeric()
  sai_values <- as.vector(outer(rep(sai_areas_surrogate, times = length(sai_areas_optimal)), sai_areas_random, "-"))/as.vector(outer(rep(sai_areas_optimal, each = length(sai_areas_surrogate)), sai_areas_random, "-"))
  sai_values <- ifelse(is.na(sai_values), 1, sai_values)
  sai_values <- ifelse(sai_values > 1, 1, sai_values)
  sai_summary <- c(quantile(sai_values, .025, na.rm = TRUE), quantile(sai_values, .5, na.rm = TRUE), quantile(sai_values, .975, na.rm = TRUE))
  names(sai_summary) <- c("SAI_lower", "SAI_median", "SAI_upper")
  return(sai_summary)  
}

##### -- plot_sai_by_area() -- #####
#### Plot the output of calculate_sai_by_area() for the three conservation algorithms
plot_sai_by_area <- function(greedy_curves, CAZ_curves, ABF_curves, area = 50, plot_vline = NULL, leg = FALSE, ...){
  greedy_sai <- calculate_sai_by_area(greedy_curves, max_area = area)
  CAZ_sai <- calculate_sai_by_area(CAZ_curves, max_area = area)
  ABF_sai <- calculate_sai_by_area(ABF_curves, max_area = area)
  #### Generate plot
  plot(greedy_sai$area, greedy_sai$SAI_upper, type = "n", xlim = c(0, area), ylim = c(0, 1), ...)
  ### Plot vertical line
  if (!is.null(plot_vline)) abline(v = plot_vline, lty = 2, col = grey(.3))
  ### Plot greedy SAI values
  lines(greedy_sai$area, greedy_sai$SAI_lower, col = grey(.6), lwd = 2, lty = 1, cex = 1.1)
  ### Plot greedy SAI values
  lines(CAZ_sai$area, CAZ_sai$SAI_lower, col = grey(.35), lwd = 2, lty = 2, cex = 1.1)
  ### Plot greedy SAI values
  lines(ABF_sai$area, ABF_sai$SAI_lower, col = grey(.1), lwd = 2, lty = 3, cex = 1.1)
  ### Legend
  if (isTRUE(leg)) {
    legend("topright", c("Maximize overall diversity", "Maximize overall rare diversity", "Maximize locally rare diversity"),
           col = c(grey(.6), grey(.35), grey(.1)),
           lwd = c(2, 2, 2),
           lty = c(1, 2, 3),
           bty = "n"
           )
  }
}

#################################
##### -- Metrics by area -- #####
#################################
##### Evolutionary/trait distinctiveness by area
distinct_by_area <- function(distrib_dat, species_distinct, priority_ranks, algorithm = c("greedy", "CAZ", "ABF")){
  algorithm <- match.arg(algorithm)
  #### For algorithms CAZ or ABF, merge the target matrix back to the extent of full raster extent
  if (algorithm == "CAZ" | algorithm == "ABF"){
    #### Update objects
    reference_coordinates <- priority_ranks$reference_coordinates
    priority_ranks <- priority_ranks$zonation_ranks
    distrib_dat <- left_join(reference_coordinates, distrib_dat, by = c("x", "y"))    
  }
  priority_ranks <- priority_ranks[, grep("rank1$", names(priority_ranks)):ncol(priority_ranks)] 
  out <- vector("list", ncol(priority_ranks))
  for (i in 1:length(out)){
    distrib_reordered <- distrib_dat[priority_ranks[, i], ]
    distrib_reordered_cumulative <- apply(distrib_reordered, 2, cumsum)
    distrib_reordered_cumulative <- apply(distrib_reordered_cumulative == 0, 2, as.numeric)
    missing_species <- vector("list", nrow(distrib_reordered_cumulative))
    for (k in 1:nrow(distrib_reordered_cumulative)){
      missing_species[[k]] <- names(distrib_reordered_cumulative[k, -c(1:2)][which(distrib_reordered_cumulative[k, -c(1:2)] == 1)])
    }
    missing_distinct <- lapply(missing_species, function(x){
      c(summary(species_distinct[species_distinct$Species %in% x, ]$w),
        total_distinct = sum(species_distinct[species_distinct$Species %in% x, ]$w)
      )
    })
    out[[i]] <- do.call("rbind", missing_distinct)
  }
  return(out)
}

##### -- plot_sai_by_area() -- #####
#### Plot the output of calculate_sai_by_area() for the three conservation algorithms
plot_distinct_by_area <- function(distinct_greedy, total_distinct, distinct_CAZ, distinct_ABF, area = 50, measure = "total_distinct", summary_measure = max, plot_vline = NULL, leg = FALSE, ...){
  #### Extract maximum distinctiveness values
  distinct_greedy <- (apply(do.call("cbind", lapply(distinct_greedy, function(x) x[, measure])), 1, summary_measure)/sum(total_distinct$w)) * 100
  distinct_CAZ <- (apply(do.call("cbind", lapply(distinct_CAZ, function(x) x[, measure])), 1, summary_measure)/sum(total_distinct$w)) * 100
  distinct_ABF <- (apply(do.call("cbind", lapply(distinct_ABF, function(x) x[, measure])), 1, summary_measure)/sum(total_distinct$w)) * 100
  #### Create a percentage area object
  area_greedy <- (1:length(distinct_greedy)/length(distinct_greedy)) * 100
  area_CAZ <- (1:length(distinct_CAZ)/length(distinct_CAZ)) * 100
  #### Generate plot
  plot(area_CAZ, distinct_CAZ, type = "n", xlim = c(0, area), ...)
  ### Plot vertical line
  if (!is.null(plot_vline)) abline(v = plot_vline, lty = 2, col = grey(.3))
  ### Plot greedy SAI values
  lines(area_greedy, distinct_greedy, col = grey(.6), lwd = 2, lty = 1, cex = 1.1)
  ### Plot greedy SAI values
  lines(area_CAZ, distinct_CAZ, col = grey(.35), lwd = 2, lty = 2, cex = 1.1)
  ### Plot greedy SAI values
  lines(area_CAZ, distinct_ABF, col = grey(.1), lwd = 2, lty = 3, cex = 1.1)
  ### Legend
  if (isTRUE(leg)) {
    legend("topright", c("Maximize overall diversity", "Maximize overall rare diversity", "Maximize locally rare diversity"),
           col = c(grey(.6), grey(.35), grey(.1)),
           lwd = c(2, 2, 2),
           lty = c(1, 2, 3),
           bty = "n"
    )
  }
}

##### Evolutionary/trait metric by area
metric_by_area <- function(distrib_dat, priority_ranks, phylo_tree, trait_tree, algorithm = c("greedy", "CAZ", "ABF")){
  algorithm <- match.arg(algorithm)
  #### For algorithms CAZ or ABF, merge the target matrix back to the extent of full raster extent
  if (algorithm == "CAZ" | algorithm == "ABF"){
    #### Update objects
    reference_coordinates <- priority_ranks$reference_coordinates
    priority_ranks <- priority_ranks$zonation_ranks
    distrib_dat <- left_join(reference_coordinates, distrib_dat, by = c("x", "y"))    
  }
  priority_ranks <- priority_ranks[, grep("rank1$", names(priority_ranks)):ncol(priority_ranks)] 
  out <- vector("list", ncol(priority_ranks))
  for (i in 1:length(out)){
    distrib_reordered <- distrib_dat[priority_ranks[, i], ]
    distrib_reordered <- as.matrix(distrib_reordered[, -c(1:2)])
    phylo_comp_data <- comparative.comm(phylo_tree, distrib_reordered)
    trait_comp_data <- comparative.comm(trait_tree, distrib_reordered)
    out[[i]] <- data.frame(area = (1:nrow(priority_ranks)/nrow(priority_ranks)) * 100,
                           SR = rowSums(distrib_reordered),
                           PD = .pd(phylo_comp_data),
                           FD = .pd(trait_comp_data)
    )
  }
  return(out)
}

##############################################
##### -- Important conservation areas -- #####
##############################################
##### Functions to calculate the percentage of target biodiversity represented within areas of conservation importance 
##### (i.e. Protected Areas, Biodiversity Hotspots, Endemic Bird Areas, G200 Areas)
##### Each function requires the following arguments: a site-by-biodiversity target matrix, a site-by-conservation priority matrix and a range of thresholds for inclusion/exclusion
##### -- protected_targets_species() -- #####
#### Percentage of species represented within each conservation area
protected_targets_species <- function(target_data, cons_priorities, thresholds = c(1, seq(10, 100, 10))){
  target_data <- left_join(target_data, cons_priorities, by = c("x", "y"))
  out <- data.frame(threshold = thresholds, perc_area_protected = NA, num_target_protected = NA, perc_target_protected = NA)
  for (i in 1:length(thresholds)){
    target_data_protected <- target_data[target_data[, names(cons_priorities)[3]] >= thresholds[i], -c(1:2, ncol(target_data))]
    perc_area_protected <- (nrow(target_data_protected)/nrow(target_data)) * 100
    t_all <- names(target_data_protected)
    t_protected <- names(target_data_protected)[which(colSums(target_data_protected) > 0)]
    t_unprotected <- setdiff(t_all, t_protected)
    perc_t_protected <- (length(t_protected)/length(t_all)) * 100
    out$perc_area_protected[i] <- perc_area_protected
    out$num_target_protected[i] <- length(t_protected)
    out$perc_target_protected[i] <- perc_t_protected       
  }                
  return(out)
}

#### -- protected_targets_traits() -- #####
#### Percentage of species represented within each conservation area
protected_targets_traits <- function(target_data, cons_priorities, thresholds = c(1, seq(10, 100, 10))){
  #### Calculate range sizes
  range_sizes <- colSums(target_data[, -c(1:2)])
  #### Combine conservation priorities data
  target_data <- left_join(target_data, cons_priorities, by = c("x", "y"))
  out <- data.frame(threshold = thresholds, perc_area_protected = NA, perc_target_protected_greedy = NA, perc_target_protected_CAZ = NA, perc_target_protected_ABF = NA)
  for (i in 1:length(thresholds)){
    #### Only targets within grid cells with >= threshold percentage area of each conservation priority are represented
    ### Calculate level of representation based on greedy algorithm target
    target_data_protected <- target_data[target_data[, names(cons_priorities)[3]] >= thresholds[i], -c(1:2, ncol(target_data))]
    perc_area_protected <- (nrow(target_data_protected)/nrow(target_data)) * 100
    perc_target_protected_greedy <- perc_target_protected_CAZ <- traits_CAZ(target_data_protected)
    ### Calculate level of representation based on ABF algorithm target
    perc_target_protected_ABF <- traits_ABF(target_data_protected, ranges = range_sizes)
    ### Populate output
    out$perc_area_protected[i] <- perc_area_protected
    out$perc_target_protected_greedy[i] <- perc_target_protected_greedy[length(perc_target_protected_greedy)] 
    out$perc_target_protected_CAZ[i] <- perc_target_protected_CAZ[length(perc_target_protected_CAZ)] 
    out$perc_target_protected_ABF[i] <- perc_target_protected_ABF[length(perc_target_protected_ABF)]
  }                
  return(out)
}

#### -- protected_targets_phylo() -- #####
#### Percentage of species represented within each conservation area
protected_targets_phylo <- function(target_data, cons_priorities, thresholds = c(1, seq(10, 100, 10)), branch_lengths){
  #### Calculate range sizes
  range_sizes <- colSums(target_data[, -c(1:2)]) 
  #### Combine conservation priorities data
  target_data <- left_join(target_data, cons_priorities, by = c("x", "y"))
  out <- data.frame(threshold = thresholds, perc_area_protected = NA, perc_target_protected_greedy = NA, perc_target_protected_CAZ = NA, perc_target_protected_ABF = NA)
  for (i in 1:length(thresholds)){
    #### Only targets within grid cells with >= threshold percentage area of each conservation priority are represented
    target_data_protected <- target_data[target_data[, names(cons_priorities)[3]] >= thresholds[i], -c(1:2, ncol(target_data))]
    ### Calculate level of representation based on greedy algorithm target
    # target_data_protected <- target_data_protected[which(colSums(target_data_protected) > 0)]
    perc_area_protected <- (nrow(target_data_protected)/nrow(target_data)) * 100
    t_all <- names(target_data_protected)
    t_protected <- names(target_data_protected)[which(colSums(target_data_protected) > 0)]
    t_unprotected <- setdiff(t_all, t_protected)
    perc_t_protected <- (length(t_protected)/length(t_all)) * 100
    perc_target_protected_greedy <- perc_t_protected
    ### Calculate level of representation based on CAZ algorithm target
    perc_target_protected_CAZ <- phylo_CAZ(target_data_protected, branches = branch_lengths)
    ### Calculate level of representation based on ABF algorithm target
    perc_target_protected_ABF <- phylo_ABF(target_data_protected, ranges = range_sizes, branches = branch_lengths)
    ### Populate output
    out$perc_area_protected[i] <- perc_area_protected
    out$perc_target_protected_greedy[i] <- perc_target_protected_greedy[length(perc_target_protected_greedy)] 
    out$perc_target_protected_CAZ[i] <- perc_target_protected_CAZ[length(perc_target_protected_CAZ)] 
    out$perc_target_protected_ABF[i] <- perc_target_protected_ABF[length(perc_target_protected_ABF)]
  }                
  return(out)
}


##### -- map_totals() -- #####
##### Map matrix row totals
map_totals <- function(cell_matrix, grid = grid_1_cea){
  r <- rasterFromXYZ(data.frame(cell_matrix[, c(1:2)], rowSums(cell_matrix[, -c(1:2)], na.rm = TRUE)))
  map_raster(r, study_grid = grid)
  return(r)
}

#### Map cell locations 
map_cells <- function(focal_cells = NA, study_grid = americas_grid_cea){
  raster::plot(study_grid, col = "grey")
  cell_xy <- xyFromCell(study_grid, focal_cells)
  points(cell_xy, col = "deepskyblue", pch = 16, cex = 1.5)
}

#### Map individual species from site by species matrix 
map_range <- function(focal_sp, 
                      sp_matrix, 
                      study_grid,
                      out_name = NULL,
                      out_save = FALSE){
  require(raster)
  species_raster <- study_grid
  species_raster[] <- sp_matrix[, focal_sp]
  species_raster[][which(is.na(study_grid[]))] <- NA
  #raster::plot(study_grid, col = "grey", legend = FALSE)
  raster::plot(species_raster, col = c(grey(0.8), "deepskyblue"))
  if (isTRUE(out_save)){
    pdf(paste("output/", out_name, ".pdf"))
    plot(study_grid, col = "grey")
    plot(species_raster, col = "deepskyblue", add = TRUE)
    dev.off()
  }
}

#### Map measure(s) aross the study area given a data.frame with ncol >= 3 with coordinates and values
map_measures <- function(out_matrix = NA,
                         study_grid = grid_1_cea,
                         crs = cea_crs, ...){
  study_grid <- as.data.frame(study_grid, xy = TRUE)[, 1:2]
  measure_df <- left_join(study_grid, out_matrix, by = c("x", "y"))
  nvar <- ncol(measure_df) - 2
  if (nvar > 1){
    out_rasters <- vector("list", nvar)
    for (i in 1:nvar){
      out_rasters[[i]] <- rasterFromXYZ(measure_df[, c(1:2, i + 2)], crs = crs)
      map_raster(out_rasters[[i]], ...)
    }
    names(out_rasters) <- names(measure_df[, -c(1:2)])
    return(out_rasters)
  } else {
    out_raster <- rasterFromXYZ(measure_df, crs = crs)
    map_raster(out_raster, ...)
    return(out_raster)
  }
}

#### Identify species not contained in specified grid cells
missing_species_from_cells <- function(cells, distrib_dat = distrib_tetrapods_1){
  species_in_cells <- distrib_dat[cells, -c(1:2)] %>% colSums()
  species_in_cells <- names(species_in_cells[species_in_cells > 0])
  species_not_in_cells <- setdiff(names(distrib_dat)[-c(1:2)], species_in_cells)
  return(species_not_in_cells)
}
