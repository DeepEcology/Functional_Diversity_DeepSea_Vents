
###############################################################################################
### Script of "High functional vulnerability across deep-sea hydrothermal vent communities" ###
### Ecologycal analyses                                                                     ###
### Joan M. Alfaro (jmalfarolucas@gmail.com)                                                ###
### Department of Biology, University of Victoria (Canada)                                  ###
### November 2023                                                                           ### 
###############################################################################################

# Set directory
setwd()

# Load R packages
library("mFD") # functional trait analyses
library("vegan") # species diversity analyses
library("betapart") # functional B-diversity
library("picante") # null models
library("abind") # arrays

# Import dataset and data preparation
#####################################

sfDVent_sp <- read.csv() # Import sfDVent_sp_JA.csv dataset

sp_comm <- as.data.frame(t(sfDVent_sp[,c(7:23)])) # species occurrences data
sp_traits <- sfDVent_sp[,c(1:6)] # species traits data
trait_info <- data.frame(trait_name = colnames(sp_traits), trait_type = c("O","O","N","O","O","N")) # Trait info data

sp_traits$Relative.Adult.Mobility <- as.ordered(sp_traits$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits$Estimated.Max.Body.Size <- as.ordered(sp_traits$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits$Habitat.Complexity <- as.factor(sp_traits$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits$Chemosynthesis.Obligate <- as.ordered(sp_traits$Chemosynthesis.Obligate) # Chemo Obligate as an ordered trait
sp_traits$Zonation.From.Vent <- ordered(sp_traits$Zonation.From.Vent, levels = c("Low", "Medium", "High")) # Zonation as an ordered trait and reorder levels
sp_traits$Feeding.Mode <- as.factor(sp_traits$Feeding.Mode) # Feeding as a factor trait

# Determining vent regions
##########################

sp_dist <- betadiver(sp_comm, method = "sim") # Sorensen dissimilarity of species 
sp_comm_hc <- hclust(sp_dist, method = "average") # Cluster analysis using the algorithm
plot(sp_comm_hc, hang = -1, main = "Vent biogeographical regions", ylab= "", xlab= "") # Plot dendogram
abline(v = 0, h = 0.75, col ="red", lty = 2) # Dissimilarity threshold set at 75
points(x = c(1:17), y = rep(0,17), pch = 19, cex = 1.5,
       col = c("#b989f0","#FFA500","#753232", "#753232", "#CD2626","#CD2626","#CD2626",
               "#c95fa4", "#3acd88", "#3A5FCD", "#3A5FCD", "#d1cb54","#242323", "#54cfd1", 
               "#03520c", "#03520c","#03520c")) # Add point and color to highlight different vent biogeoprovinces

# Regroup by region
prov <- c("IndR", "ESR", "EPR", "GoC", "JdF", "Kerm", "SWP", "SWP","Izu-Mar",
          "MohR", "EPR", "MAR", "SWP", "Okin", "EPR", "MAR", "IndR") # to add as column for re-classification of vents

sp_comm$Province <- prov # add column province
sp_comm_biog = sapply(sp_comm[-512], tapply, INDEX = sp_comm$Province, sum) # create new table with species classified in bioprovinces
sp_comm <- sp_comm[,-512] # delete column of provinces
sp_comm_biog <- decostand(sp_comm_biog, "pa") # presence/abscense transformation (1/0)

## Multidimensional functional spaces and indices
#################################################

# Functional distance between species
func_dist_sp <- funct.dist(sp_traits, 
                           tr_cat = trait_info, 
                           metric = "gower", 
                           ordinal_var = "metric",
                           weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual <- quality.fspaces(func_dist_sp, 
                             deviation_weighting = c("absolute", "squarred"),
                             maxdim_pcoa = 9,
                             fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual$quality_fspaces, 2, which.min) # best number of dimesions calculated with MAD index (4D) (Maire et al. 2015)
sp_faxes <- func_qual$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices of vent regions
func_ind_sp <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                 asb_sp_w = sp_comm_biog,
                                 scaling = TRUE,
                                 details_returned = T,
                                 verbose = T) # mutidimensinal functional diversity indices in a 4D functional space of vent regions


# Functionl B-diversity 
fbdiv <- functional.betapart.core.pairwise(sp_comm_biog, sp_faxes[,c(1:4)]) # quantities needed for computing pairwise functional B-diversity 
fb <- functional.beta.pair(fbdiv, index.family = "jaccard") # Functional B-diversity (b-diversity, tunrover, nestedness), Jaccard index
fb$funct.beta.jtu # Functional turnover, Jaccard index

## Functional entities or Unique trait combinations (UTC) and derived functional indices
########################################################################################

func_ent <- sp.to.fe(sp_traits, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent) # compute UTC derived indices
func_ind_fe$asb_fdfe # UTC Functional indexes

## Functional diversity indices summary
#######################################

func_ind <- data.frame(nb_sp= func_ind_fe$asb_fdfe[,"nb_sp"]) # species richness
func_ind$Pnb_sp <- round((func_ind$nb_sp/511)*100, 4) # species richness percentage
func_ind$fdis <- func_ind_sp$functional_diversity_indices$fdis # functional dispersion 
func_ind$fred <- func_ind_fe$asb_fdfe[,"fred"] # functional redundancy
func_ind$fvuln <- func_ind_fe$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind,2), "func_ind.csv")  

# Simulate species extinctions effects on functional richness 
sim_obs <- extinct.sim(comm = sp_comm_biog, nperm = 30) # load extinct.sim function


## Null models
##############

n_perm = 1000 # number of permutations

multif_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (FDis)
fe_perm <- as.list(rep(NA, n_perm)) # list to save permutations results of FE indices (FRed, FVul)
fb_jac_perm <- as.list(rep(NA, n_perm)) # list to save Jaccard's B-diversity turnover permutations results

for(i in seq(n_perm)){
  
  sp_comm_n = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n <- alpha.fd.multidim(sp_faxes_coord = sp_faxes[,c(1:4)],
                                asb_sp_w = sp_comm_n,
                                ind_vect= "fdis",
                                scaling = T,
                                details_returned = T,
                                verbose = F) #  Null FDis
  
  multif_perm[[i]] <- multif_n$functional_diversity_indices["fdis"] # store null FDis results 
  
  fbdiv_n <- functional.betapart.core.pairwise(sp_comm_n, sp_faxes[,c(1:4)]) # quantities needed for computing the null pairwise functional B-diversity, Jaccard index
  fb_n <- functional.beta.pair(fbdiv_n, index.family = "jaccard") # Null functional B-diversity indices, Jaccard index
  fb_jac_perm[[i]] <- fb_n$funct.beta.jtu # store null Jaccard functional turnover results
  
  func_ind_fe_n <- alpha.fd.fe(sp_comm_n, sp_to_fe = func_ent) # Null FE indices
  fe_perm[[i]] <- func_ind_fe_n$asb_fdfe[,c("fred", "fvuln")] # store null FRed and FVuln indices result
  
} # null model for functional indexes

# Wrap null model results
#########################

multif_perm_a <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for null FDis results  
fe_perm_a <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for FRed and FVul inices results 
fb_jac_perm_a <- array(NA, dim = c(11,11,0)) # create an array with 11 col & rows for functional B-diversity results

for (i in c(1:n_perm)) { 
  multif_perm_a <- abind(multif_perm_a, multif_perm[[i]]) # FRic & FDis
  fe_perm_a <- abind(fe_perm_a, fe_perm[[i]]) # FE indices
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a), colnames(multif_perm_a))) # matrix to save mean values of FDis permutations
multif_perm_sd <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a), colnames(multif_perm_a))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in "fdis") {
    multif_perm_mean[i,j] = mean(multif_perm_a[i,j,])
    multif_perm_sd[i,j] = sd(multif_perm_a[i,j,])
  }
} # mean and sd values ofFDis permutations 

fe_perm_mean <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a), colnames(fe_perm_a))) # matrix to save mean values of FRed and FVuln indices permutations
fe_perm_sd <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a), colnames(fe_perm_a))) # matrix to save sd values of FRed and FVuln indices permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean[i,j] = mean(fe_perm_a[i,j,])
    fe_perm_sd[i,j] = sd(fe_perm_a[i,j,])
  }
} # mean and sd values of FE indices permutations  

for (i in c(1:n_perm)) { 
  fb_jac_mat <- as.matrix(fb_jac_perm[[i]]) # functional B-diversity distance as matrix
  fb_jac_perm_a <- abind(fb_jac_perm_a, fb_jac_mat) # functional B-diversiy
} # merge each elemnt of the lists as dimensions of the arrays

fb_jac_perm_mean = matrix(NA, ncol = 11, nrow= 11, dimnames = list(rownames(fb_jac_mat), colnames(fb_jac_mat))) # matrix to save mean values of B-div indices permutations
fb_jac_perm_sd = matrix(NA, ncol = 11, nrow= 11, dimnames = list(rownames(fb_jac_mat), colnames(fb_jac_mat))) # # matrix to save sd values of B-div indices permutations 
for (i in 1:11) {
  for (j in 1:11) {
    fb_jac_perm_mean[i,j] = mean(fb_jac_perm_a[i,j,])
    fb_jac_perm_sd[i,j] = sd(fb_jac_perm_a[i,j,])
  }
} # mean and sd values of functional B-diversity permutations

# Standardized effect size (SES)
ses_multif <- (func_ind_sp$functional_diversity_indices[,c("fdis")] - multif_perm_mean) / multif_perm_sd # SES functional multidimensional alpha diversity indices
ses_fe <- (func_ind_fe$asb_fdfe[,c("fred", "fvuln")] - fe_perm_mean) / fe_perm_sd # SES FE indices
ses_fbjac <- ((as.matrix(fb$funct.beta.jtu)) - fb_jac_perm_mean) / fb_jac_perm_sd # SES functional Jaccard's B-diversity

###############
### The end ###
###############
