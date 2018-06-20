##############################################################
##### -- Good surrogacy among biodiversity dimensions -- #####
##############################################################
################### SUPPLEMENTARY ANALYSIS ###################
########################################################################
##### -- Checking the taxonomic distribution of missing species -- #####
########################################################################
##### -- Load full list of tetrapod species of the Americas -- #####
tetrapods_americas_list <- readRDS("data/tetrapods_americas_list.rds")

#### What proportion of Americas species was included in the analysis?
nrow(subset(tetrapods_americas_list, included_in_analysis == TRUE))/nrow(tetrapods_americas_list)

#### Generate plot of species/family proportions with outliers
tiff("figures/families-missing-from-analysis.tif", width = 6.0, height = 6.0, units = 'in', res = 600)
plot(
  as.numeric(table(subset(tetrapods_americas_list, included_in_analysis == TRUE | included_in_analysis == FALSE)$family)/nrow(subset(tetrapods_americas_list, included_in_analysis == TRUE | included_in_analysis == FALSE))),
  as.numeric(table(subset(tetrapods_americas_list, included_in_analysis == FALSE)$family)/nrow(subset(tetrapods_americas_list, included_in_analysis == FALSE))),
  xlab = "Proportion of species in each family in the Americas",
  ylab = "Proportion of species in each family excluded from analyses",
  pch = 16, 
  cex = 0.8
)
abline(a = 0, b = 1, lty = 2)
text(as.numeric(table(subset(tetrapods_americas_list, included_in_analysis == TRUE | included_in_analysis == FALSE)$family)/nrow(subset(tetrapods_americas_list, included_in_analysis == TRUE | included_in_analysis == FALSE)))[c(246, 209, 88, 38, 213, 219, 144)],
     as.numeric(table(subset(tetrapods_americas_list, included_in_analysis == FALSE)$family)/nrow(subset(tetrapods_americas_list, included_in_analysis == FALSE)))[c(246, 209, 88, 38, 213, 219, 144)],
     labels = names(table(tetrapods_americas_list$family))[c(246, 209, 88, 38, 213, 219, 144)],
     adj = -0.1,
     col = "tomato"
     )
dev.off()

# What is the proportion of missing families?
length(setdiff(unique(subset(tetrapods_americas_list, included_in_analysis == FALSE)$family), unique(subset(tetrapods_americas_list, included_in_analysis == TRUE)$family)))/
length(levels(tetrapods_americas_list$family)[order(levels(tetrapods_americas_list$family))])

#Which families are missing from the analysis?
setdiff(unique(subset(tetrapods_americas_list, included_in_analysis == FALSE)$family), unique(subset(tetrapods_americas_list, included_in_analysis == TRUE)$family))


