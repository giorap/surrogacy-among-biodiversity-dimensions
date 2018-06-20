##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### SUPPLEMENTARY ANALYSIS ###################
#########################################################
##### -- Testing for influence of imputed values -- #####
#########################################################
##### Load library
library(caper)
##### Quantify phylogenetic signal of traits with and without imputed taxa
#### Mammals
mamm_comparative_data_full <- comparative.data(data = traits_mamm_complete, phy = tree_mamm_complete, names.col = "species", na.omit = TRUE)
tree_mamm_noImputed <- drop.tip(tree_mamm_complete, traits_imputed_mamm$species)
mamm_comparative_data_noImputed <- comparative.data(data = traits_mamm_noImputed, phy = tree_mamm_noImputed, names.col = "species", na.omit = TRUE)

mamm_full_phylosig <- pgls(body_mass_log ~ 1, mamm_comparative_data_full, lambda='ML')
saveRDS(mamm_full_phylosig, "rds/metrics/mamm_full_phylosig.rds")
mamm_noImputed_phylosig <- pgls(body_mass_log ~ 1, mamm_comparative_data_noImputed, lambda='ML')
saveRDS(mamm_noImputed_phylosig, "rds/metrics/mamm_noImputed_phylosig.rds")

#### Reptiles
##### Quantify phylogenetic signal of traits with and without imputed taxa
rept_comparative_data_full <- comparative.data(data = traits_rept_complete, phy = tree_rept_complete, names.col = "species", na.omit = TRUE)
tree_rept_noImputed <- drop.tip(tree_rept_complete, traits_imputed_rept$species)
rept_comparative_data_noImputed <- comparative.data(data = traits_rept_noImputed, phy = tree_rept_noImputed, names.col = "species", na.omit = TRUE)

rept_full_phylosig <- pgls(body_mass_log ~ 1, rept_comparative_data_full, lambda='ML')
saveRDS(rept_full_phylosig, "rds/metrics/rept_full_phylosig.rds")
rept_noImputed_phylosig <- pgls(body_mass_log ~ 1, rept_comparative_data_noImputed, lambda='ML')
saveRDS(rept_noImputed_phylosig, "rds/metrics/rept_noImputed_phylosig.rds")

#### Birds
##### Quantify phylogenetic signal of traits with and without imputed taxa
bird_comparative_data_full <- comparative.data(data = traits_bird_complete, phy = tree_bird_complete, names.col = "species", na.omit = TRUE)
tree_bird_noImputed <- drop.tip(tree_bird_complete, traits_imputed_bird$species)
bird_comparative_data_noImputed <- comparative.data(data = traits_bird_noImputed, phy = tree_bird_noImputed, names.col = "species", na.omit = TRUE)

bird_full_phylosig <- pgls(body_mass_log ~ 1, bird_comparative_data_full, lambda='ML')
saveRDS(bird_full_phylosig, "rds/metrics/bird_full_phylosig.rds")
bird_noImputed_phylosig <- pgls(body_mass_log ~ 1, bird_comparative_data_noImputed, lambda='ML')
saveRDS(bird_noImputed_phylosig, "rds/metrics/bird_noImputed_phylosig.rds")
