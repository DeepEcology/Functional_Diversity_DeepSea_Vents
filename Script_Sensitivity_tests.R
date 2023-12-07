
###############################################################################################
### Script of "High functional vulnerability across deep-sea hydrothermal vent communities" ###
### Sensitivity analyses - Modality reduction & trait reduction                             ###
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

####################################################
## Sensitivity Analyses: trait modality reduction ##
####################################################

# Import dataset and data preparation
sfDVent_sp <- read.csv() # Import sfDVent_sp_JA.csv dataset
sp_traits_sen <- sfDVent_sp[,c(1:6)] # species traits
sp_comm <- as.data.frame(t(sfDVent_sp[,c(7:23)])) # species occurrences data
trait_info <- data.frame(trait_name = colnames(sp_traits_sen), trait_type = c("O","O","N","O","O","N")) # Trait info data

# Regroup by region
prov <- c("IndR", "ESR", "EPR", "GoC", "JdF", "Kerm", "SWP", "SWP","Izu-Mar",
          "MohR", "EPR", "MAR", "SWP", "Okin", "EPR", "MAR", "IndR") # to add as column for re-classification of vents

sp_comm$Province <- prov # add column province
sp_comm_biog = sapply(sp_comm[-512], tapply, INDEX = sp_comm$Province, sum) # create new table with species classified in bioprovinces
sp_comm <- sp_comm[,-512] # delete column of provinces
sp_comm_biog <- decostand(sp_comm_biog, "pa") # presence/abscense transformation (1/0)

# Trait modality reduction
sp_traits_sen[c(which(sp_traits_sen$Chemosynthesis.Obligate != "No")), "Chemosynthesis.Obligate"] <- "CBE" # Chemo-obligate 2 categories
sp_traits_sen[c(which(sp_traits_sen$Habitat.Complexity != "Does not add")), "Habitat.Complexity"] <- "Add" # Habitat complexity 2 categories
sp_traits_sen[c(which(sp_traits_sen$Zonation.From.Vent != "Low")), "Zonation.From.Vent"] <- "High" # Zonation of vent 2 categories

sp_traits_sen$Relative.Adult.Mobility <- as.ordered(sp_traits_sen$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits_sen$Estimated.Max.Body.Size <- as.ordered(sp_traits_sen$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits_sen$Habitat.Complexity <- as.factor(sp_traits_sen$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits_sen$Chemosynthesis.Obligate <- ordered(sp_traits_sen$Chemosynthesis.Obligate, levels = c("No", "CBE")) # Chemo Obligate as an ordered trait and reorder levels
sp_traits_sen$Zonation.From.Vent <- ordered(sp_traits_sen$Zonation.From.Vent, levels = c("Low", "High")) # Zonation as an ordered trait and reorder levels
sp_traits_sen$Feeding.Mode <- as.factor(sp_traits_sen$Feeding.Mode) # Feeding as a factor trait

# Functional distance between species
func_dist_sp_sen <- funct.dist(sp_traits_sen, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_sen <- quality.fspaces(func_dist_sp_sen, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_sen$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_sen <- func_qual_sen$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_sen <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_sen[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_sen$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities (FE) or unique trait combinations (UTC) and derived functional indices
func_ent_sen <- sp.to.fe(sp_traits_sen, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_sen <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_sen) # compute UTC derived indices
func_ind_fe_sen$asb_fdfe # UTC functional indeces

# Summary sensivity indices
func_ind_sen <- data.frame(nb_sp= func_ind_fe_sen$asb_fdfe[,"nb_sp"]) # species richness
func_ind_sen$fdis <- func_ind_sp_sen$functional_diversity_indices$fdis # functional dispersion 
func_ind_sen$fred <- func_ind_fe_sen$asb_fdfe[,"fred"] # functional redundancy
func_ind_sen$fvuln <- func_ind_fe_sen$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_sen,2), "func_ind_sen.csv")  

# Null models for functional dispersion, redundancy and vulnerability
n_perm = 1000 # number of permutations

multif_perm_sen <- as.list(rep(NA, n_perm)) # list to save permutations results of Fdis
fe_perm_sen <- as.list(rep(NA, n_perm)) # list to save permutations results of Fred, Fvul

for(i in seq(n_perm)){
  
  sp_comm_n_sen = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_sen <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_sen[,c(1:4)],
                                    asb_sp_w = sp_comm_n_sen,
                                    ind_vect= "fdis",
                                    scaling = T,
                                    details_returned = T,
                                    verbose = F) #  Null FDis
  
  multif_perm_sen[[i]] <- multif_n_sen$functional_diversity_indices["fdis"] # store Nnull FDis results 
  
  func_ind_fe_n_sen <- alpha.fd.fe(sp_comm_n_sen, sp_to_fe = func_ent_sen) # Null UTC indices
  fe_perm_sen[[i]] <- func_ind_fe_n_sen$asb_fdfe[,c("fred","fvuln")] # store null FRed and FVuln indices result
  
} # null model of functional indices

# Wrap null model results
multif_perm_a_sen <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for null FDis results  
fe_perm_a_sen <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for UTC inices results 

for (i in c(1:n_perm)) { 
  multif_perm_a_sen <- abind(multif_perm_a_sen, multif_perm_sen[[i]]) # FDis
  fe_perm_a_sen <- abind(fe_perm_a_sen, fe_perm_sen[[i]]) # UTC indices
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_sen <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_sen), colnames(multif_perm_a_sen))) # matrix to save mean values of FDis permutations
multif_perm_sd_sen <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_sen), colnames(multif_perm_a_sen))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_sen[i,j] = mean(multif_perm_a_sen[i,j,])
    multif_perm_sd_sen[i,j] = sd(multif_perm_a_sen[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_sen <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_sen), colnames(fe_perm_a_sen))) # matrix to save mean values of UTC indices permutations
fe_perm_sd_sen <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_sen), colnames(fe_perm_a_sen))) # matrix to save sd values of UTC indices permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_sen[i,j] = mean(fe_perm_a_sen[i,j,])
    fe_perm_sd_sen[i,j] = sd(fe_perm_a_sen[i,j,])
  }
} # mean and sd values of FE indices permutations

# Standardized effect size (SES)
ses_multif_sen <- (func_ind_sp_sen$functional_diversity_indices[,"fdis"] - multif_perm_mean_sen) / multif_perm_sd_sen # SES FDis
ses_fe_sen <- (func_ind_fe_sen$asb_fdfe[,c("fred", "fvuln")] - fe_perm_mean_sen) / fe_perm_sd_sen # SES UTC indices

# Pariwise functional turnover null model 
ses_turn_sen <- matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog)))) # matrix to store SES functional turnover results
n_perm = 1000 # number of permutations

for (k in rownames(ses_turn_sen)) {
  for (i in rownames(ses_turn_sen[-(which(rownames(ses_turn_sen) == k)),])) {
    sp_comm_sub_sen <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_sen <- sp_faxes_sen[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_sen <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_sen = randomizeMatrix(sp_comm_sub_sen, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_sen <- functional.beta.pair(sp_comm_n_sen, sp_faxes_sub_sen, index.family= "jaccard")
      fb_turn_perm_sen[j, ] <- func_beta_jacc_sen$funct.beta.jtu
    }
    
    fb_turn_obs_sen <- functional.beta.pair(sp_comm_sub_sen, sp_faxes_sub_sen, index.family= "jaccard")
    ses_fb_turn_sen <- (fb_turn_obs_sen$funct.beta.jtu - mean(fb_turn_perm_sen))/ sd(fb_turn_perm_sen)
    ses_turn_sen[k, i] <- ses_fb_turn_sen
  }
} # Functional turnover SES 

#write.csv(ses_turn_sen, "ses_fb_turn_sen.csv") 


###########################################
## Sensitivity Analyses: Trait Reduction ##
###########################################

# Mobility
##########

# Import dataset and data preparation
sp_traits_mob <- sfDVent_sp[,c(2:6)] # species traits w/o Mobility 
trait_info <- data.frame(trait_name = colnames(sp_traits_size), trait_type = c("O","N","O","O","N")) # Trait info data

sp_traits_mob$Estimated.Max.Body.Size <- as.ordered(sp_traits_mob$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits_mob$Habitat.Complexity <- as.factor(sp_traits_mob$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits_mob$Chemosynthesis.Obligate <- as.ordered(sp_traits_mob$Chemosynthesis.Obligate) # Chemo Obligate as an ordered trait
sp_traits_mob$Zonation.From.Vent <- ordered(sp_traits_mob$Zonation.From.Vent, levels = c("Low", "Medium", "High")) # Zonation as an ordered trait and reorder levels
sp_traits_mob$Feeding.Mode <- as.factor(sp_traits_mob$Feeding.Mode) # Feeding as a factor trait

# Functional distance between species
func_dist_sp_mob <- funct.dist(sp_traits_mob, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_mob <- quality.fspaces(func_dist_sp_mob, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_mob$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_mob <- func_qual_mob$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_mob <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_mob[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_mob$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities or unique trait combinations (UTC) and derived functional indices
func_ent_mob <- sp.to.fe(sp_traits_mob, trait_info[-1,], fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_mob <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_mob) # compute UTC derived indices
func_ind_fe_mob$asb_fdfe # UTC indices

# Summary sensivity indices
func_ind_mob <- data.frame(nb_fe = func_ind_fe_mob$asb_fdfe[,"nb_sp"]) # sp richness
func_ind_mob$fdis <- func_ind_sp_mob$functional_diversity_indices$fdis # functional dispersion 
func_ind_mob$fred <- func_ind_fe_mob$asb_fdfe[,"fred"] # functional redundancy
func_ind_mob$fvuln <- func_ind_fe_mob$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_mob,2), "func_ind_mob.csv")  

# Null models
n_perm = 1000 # number of permutations

multif_perm_mob <- as.list(rep(NA, n_perm)) # list to save permutations results of Fdis
fe_perm_mob <- as.list(rep(NA, n_perm)) # list to save permutations results of Fred and Fvul

for(i in seq(n_perm)){
  
  sp_comm_n_mob = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_mob <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_mob[,c(1:4)],
                                    asb_sp_w = sp_comm_n_mob,
                                    ind_vect= "fdis",
                                    scaling = T,
                                    details_returned = T,
                                    verbose = F) #  Null FDis
  
  multif_perm_mob[[i]] <- multif_n_mob$functional_diversity_indices["fdis"] # store Nnull FDis results 
  
  func_ind_fe_n_mob <- alpha.fd.fe(sp_comm_n_mob, sp_to_fe = func_ent_mob) # Null UTC indices
  fe_perm_mob[[i]] <- func_ind_fe_n_mob$asb_fdfe[,c("fred","fvuln")] # store null FRed and FVuln result
  
} # null model of functional indices

# Wrap null model results
multif_perm_a_mob <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for null FDis results  
fe_perm_a_mob <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for FRed and FVuln inices results 

for (i in c(1:n_perm)) { 
  multif_perm_a_mob <- abind(multif_perm_a_mob, multif_perm_mob[[i]]) # FDis
  fe_perm_a_mob <- abind(fe_perm_a_mob, fe_perm_mob[[i]]) # FRed and FVuln 
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_mob <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_mob), colnames(multif_perm_a_mob))) # matrix to save mean values of FDis permutations
multif_perm_sd_mob <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_mob), colnames(multif_perm_a_mob))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_mob[i,j] = mean(multif_perm_a_mob[i,j,])
    multif_perm_sd_mob[i,j] = sd(multif_perm_a_mob[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_mob <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_mob), colnames(fe_perm_a_mob))) # matrix to save mean values of FRed and FVuln  permutations
fe_perm_sd_mob <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_mob), colnames(fe_perm_a_mob))) # matrix to save sd values of FRed and FVuln permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_mob[i,j] = mean(fe_perm_a_mob[i,j,])
    fe_perm_sd_mob[i,j] = sd(fe_perm_a_mob[i,j,])
  }
} # mean and sd values of FE indices permutations

# Standardized effect size (SES)
ses_multif_mob <- (func_ind_sp_mob$functional_diversity_indices[,"fdis"] - multif_perm_mean_mob) / multif_perm_sd_mob # SES FDis
ses_fe_mob <- (func_ind_fe_mob$asb_fdfe[,c("fred","fvuln")] - fe_perm_mean_mob) / fe_perm_sd_mob # SES Fred and FVuln 

# Pariwise functional turnover null model 
ses_turn_mob <-matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog))))
n_perm = 1000

for (k in rownames(ses_turn_mob)) {
  for (i in rownames(ses_turn_mob[-(which(rownames(ses_turn_mob) == k)),])) {
    sp_comm_sub_mob <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_mob <- sp_faxes_mob[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_mob <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_mob = randomizeMatrix(sp_comm_sub_mob, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_mob <- functional.beta.pair(sp_comm_n_mob, sp_faxes_sub_mob, index.family= "jaccard")
      fb_turn_perm_mob[j, ] <- func_beta_jacc_mob$funct.beta.jtu
    }
    
    fb_turn_obs_mob <- functional.beta.pair(sp_comm_sub_mob, sp_faxes_sub_mob, index.family= "jaccard")
    ses_fb_turn_mob <- (fb_turn_obs_mob$funct.beta.jtu - mean(fb_turn_perm_mob))/ sd(fb_turn_perm_mob)
    ses_turn_mob[k, i] <- ses_fb_turn_mob
  }
} # Functional turnover SES 

#write.csv(ses_turn_mob, "ses_fb_turn_mob.csv") 


# Size
######

sp_traits_size <- sfDVent_sp[,c(1,3,4,5,6)] # species traits w/o Size
trait_info <- data.frame(trait_name = colnames(sp_traits_size), trait_type = c("O","N","O","O","N")) # Trait info data

sp_traits_size$Relative.Adult.Mobility <- as.ordered(sp_traits_size$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits_size$Habitat.Complexity <- as.factor(sp_traits_size$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits_size$Chemosynthesis.Obligate <- as.ordered(sp_traits_size$Chemosynthesis.Obligate) # Chemo Obligate as an ordered trait
sp_traits_size$Zonation.From.Vent <- ordered(sp_traits_size$Zonation.From.Vent, levels = c("Low", "Medium", "High")) # Zonation as an ordered trait and reorder levels
sp_traits_size$Feeding.Mode <- as.factor(sp_traits_size$Feeding.Mode) # Feeding as a factor trait

# Functional distance between species
func_dist_sp_size <- funct.dist(sp_traits_size, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_size <- quality.fspaces(func_dist_sp_size, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_size$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_size <- func_qual_size$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_size <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_size[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_size$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities or unique trait combinations (UTC) and derived functional indices
func_ent_size <- sp.to.fe(sp_traits_size, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_size <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_size) # compute UTC derived indices
func_ind_fe_size$asb_fdfe # UTC indices

# Summary sensivity indices
func_ind_size <- data.frame(nb_fe = func_ind_fe_size$asb_fdfe[,"nb_sp"]) # sp richness
func_ind_size$fdis <- func_ind_sp_size$functional_diversity_indices$fdis # functional dispersion 
func_ind_size$fred <- func_ind_fe_size$asb_fdfe[,"fred"] # functional redundancy
func_ind_size$fvuln <- func_ind_fe_size$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_size,2), "func_ind_size.csv")  

# Null models
n_perm = 1000 # number of permutations

multif_perm_size <- as.list(rep(NA, n_perm)) # list to save permutations results of Fdis
fe_perm_size <- as.list(rep(NA, n_perm)) # list to save permutations results of UTC indices

for(i in seq(n_perm)){
  
  sp_comm_n_size = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_size <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_size[,c(1:4)],
                                     asb_sp_w = sp_comm_n_size,
                                     ind_vect = "fdis",
                                     scaling = T,
                                     details_returned = T,
                                     verbose = F) #  Null FDis
  
  multif_perm_size[[i]] <- multif_n_size$functional_diversity_indices["fdis"] # store Null FDis results 
  
  func_ind_fe_n_size <- alpha.fd.fe(sp_comm_n_size, sp_to_fe = func_ent_size) # Null UTC indices
  fe_perm_size[[i]] <- func_ind_fe_n_size$asb_fdfe[,c("fred", "fvuln")] # store null FRed and FVuln result
  
} # null model of functional indices

# Wrap null model results
multif_perm_a_size <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for FDis results  
fe_perm_a_size <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for Fred and FVuln results 

for (i in c(1:n_perm)) { 
  multif_perm_a_size <- abind(multif_perm_a_size, multif_perm_size[[i]]) # FDis
  fe_perm_a_size <- abind(fe_perm_a_size, fe_perm_size[[i]]) # FRed and Fvuln
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_size <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_size), colnames(multif_perm_a_size))) # matrix to save mean values of FDis permutations
multif_perm_sd_size <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_size), colnames(multif_perm_a_size))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_size[i,j] = mean(multif_perm_a_size[i,j,])
    multif_perm_sd_size[i,j] = sd(multif_perm_a_size[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_size <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_size), colnames(fe_perm_a_size))) # matrix to save mean values of FRed and FVuln permutations
fe_perm_sd_size <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_size), colnames(fe_perm_a_size))) # matrix to save sd values of FRed and FVuln permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_size[i,j] = mean(fe_perm_a_size[i,j,])
    fe_perm_sd_size[i,j] = sd(fe_perm_a_size[i,j,])
  }
} # mean and sd values FRed and FVuln permutations

# Standardized effect size (SES)
ses_multif_size <- (func_ind_sp_size$functional_diversity_indices[,"fdis"] - multif_perm_mean_size) / multif_perm_sd_size # SES functional multidimensional alpha diversity indices
ses_fe_size <- (func_ind_fe_size$asb_fdfe[,c("fred", "fvuln")] - fe_perm_mean_size) / fe_perm_sd_size # SES FE indices

# Functional pariwise turnover null model 
ses_turn_size <-matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog))))
n_perm = 1000 # number of permutations

for (k in rownames(ses_turn_size)) {
  for (i in rownames(ses_turn_size[-(which(rownames(ses_turn_size) == k)),])) {
    sp_comm_sub_size <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_size <- sp_faxes_size[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_size <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_size = randomizeMatrix(sp_comm_sub_size, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_size <- functional.beta.pair(sp_comm_n_size, sp_faxes_sub_size, index.family= "jaccard")
      fb_turn_perm_size[j, ] <- func_beta_jacc_size$funct.beta.jtu
    }
    
    fb_turn_obs_size <- functional.beta.pair(sp_comm_sub_size, sp_faxes_sub_size, index.family= "jaccard")
    ses_fb_turn_size <- (fb_turn_obs_size$funct.beta.jtu - mean(fb_turn_perm_size))/ sd(fb_turn_perm_size)
    ses_turn_size[k, i] <- ses_fb_turn_size
  }
} # Functional turnover SES 

#write.csv(ses_turn_size, "ses_fb_turn_size.csv")  


# Feeding
#########

sp_traits_feed <- sfDVent_sp[,c(1,2,3,4,5)] # species traits w/o Feeding
trait_info <- data.frame(trait_name = colnames(sp_traits_feed), trait_type = c("O","O","N","O","O")) # Trait info data

sp_traits_feed$Relative.Adult.Mobility <- as.ordered(sp_traits_feed$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits_feed$Estimated.Max.Body.Size <- as.ordered(sp_traits_feed$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits_feed$Habitat.Complexity <- as.factor(sp_traits_feed$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits_feed$Chemosynthesis.Obligate <- as.ordered(sp_traits_feed$Chemosynthesis.Obligate) # Chemo Obligate as an ordered trait
sp_traits_feed$Zonation.From.Vent <- ordered(sp_traits_feed$Zonation.From.Vent, levels = c("Low", "Medium", "High")) # Zonation as an ordered trait and reorder levels

# Functional distance between species
func_dist_sp_feed <- funct.dist(sp_traits_feed, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_feed <- quality.fspaces(func_dist_sp_feed, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_feed$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_feed <- func_qual_feed$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_feed <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_feed[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_feed$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities or unique trait combinations (UTC) and derived functional indices
func_ent_feed <- sp.to.fe(sp_traits_feed, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_feed <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_feed) # compute UTC derived indices
func_ind_fe_feed$asb_fdfe # UTC indices

# Summary sensivity indices
func_ind_feed <- data.frame(nb_fe = func_ind_fe_feed$asb_fdfe[,"nb_sp"]) # sp richness
func_ind_feed$fdis <- func_ind_sp_feed$functional_diversity_indices$fdis # functional dispersion 
func_ind_feed$fred <- func_ind_fe_feed$asb_fdfe[,"fred"] # functional redundancy
func_ind_feed$fvuln <- func_ind_fe_feed$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_feed,2), "func_ind_feed.csv")  

# Null models
n_perm = 1000 # number of permutations

multif_perm_feed <- as.list(rep(NA, n_perm)) # list to save permutations results of Fdis
fe_perm_feed <- as.list(rep(NA, n_perm)) # list to save permutations results of Fred and Fvul)

for(i in seq(n_perm)){
  
  sp_comm_n_feed = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_feed <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_feed[,c(1:4)],
                                     asb_sp_w = sp_comm_n_feed,
                                     ind_vect = "fdis",
                                     scaling = T,
                                     details_returned = T,
                                     verbose = F) #  Null FDis
  
  multif_perm_feed[[i]] <- multif_n_feed$functional_diversity_indices["fdis"] # store Null FDis results 
  
  func_ind_fe_n_feed <- alpha.fd.fe(sp_comm_n_feed, sp_to_fe = func_ent_feed) # Null UTC indices
  fe_perm_feed[[i]] <- func_ind_fe_n_feed$asb_fdfe[, c("fred", "fvuln")] # store null FRed and FVuln result

} # null model of functional indices

# Wrap null model results
multif_perm_a_feed <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for null FDis results  
fe_perm_a_feed <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for Fred and FVuln results 

for (i in c(1:n_perm)) { 
  multif_perm_a_feed <- abind(multif_perm_a_feed, multif_perm_feed[[i]]) # FDis
  fe_perm_a_feed <- abind(fe_perm_a_feed, fe_perm_feed[[i]]) # Fred and FVuln
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_feed <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_feed), colnames(multif_perm_a_feed))) # matrix to save mean values of FDis permutations
multif_perm_sd_feed <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_feed), colnames(multif_perm_a_feed))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_feed[i,j] = mean(multif_perm_a_feed[i,j,])
    multif_perm_sd_feed[i,j] = sd(multif_perm_a_feed[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_feed <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_feed), colnames(fe_perm_a_feed))) # matrix to save mean values of FRed and FVuln  permutations
fe_perm_sd_feed <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_feed), colnames(fe_perm_a_feed))) # matrix to save sd values of FRed and FVuln permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_feed[i,j] = mean(fe_perm_a_feed[i,j,])
    fe_perm_sd_feed[i,j] = sd(fe_perm_a_feed[i,j,])
  }
} # mean and sd values of FRed and FVuln permutations

# Standardized effect size (SES)
ses_multif_feed <- (func_ind_sp_feed$functional_diversity_indices[,"fdis"] - multif_perm_mean_feed) / multif_perm_sd_feed # SES FDis
ses_fe_feed <- (func_ind_fe_feed$asb_fdfe[,c("fred","fvuln")] - fe_perm_mean_feed) / fe_perm_sd_feed # SES FRed and FVuln

# Pariwise functional turnover null model 
ses_turn_feed <-matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog))))
n_perm = 1000

for (k in rownames(ses_turn_feed)) {
  for (i in rownames(ses_turn_feed[-(which(rownames(ses_turn_feed) == k)),])) {
    sp_comm_sub_feed <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_feed <- sp_faxes_feed[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_feed <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_feed = randomizeMatrix(sp_comm_sub_feed, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_feed <- functional.beta.pair(sp_comm_n_feed, sp_faxes_sub_feed, index.family= "jaccard")
      fb_turn_perm_feed[j, ] <- func_beta_jacc_feed$funct.beta.jtu
    }
    
    fb_turn_obs_feed <- functional.beta.pair(sp_comm_sub_feed, sp_faxes_sub_feed, index.family= "jaccard")
    ses_fb_turn_feed <- (fb_turn_obs_feed$funct.beta.jtu - mean(fb_turn_perm_feed))/ sd(fb_turn_perm_feed)
    ses_turn_feed[k, i] <- ses_fb_turn_feed
  }
} # Functional turnover SES 

#write.csv(ses_turn_feed, "ses_fb_turn_feed.csv")  


# Habitat complexity
####################

sp_traits_hc <- sfDVent_sp[,c(1,2,4,5,6)] # species traits w/o Habitat Complexity
trait_info <- data.frame(trait_name = colnames(sp_traits_hc), trait_type = c("O","O","O","O","N")) # Trait info data

sp_traits_hc$Relative.Adult.Mobility <- as.ordered(sp_traits_hc$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits_hc$Estimated.Max.Body.Size <- as.ordered(sp_traits_hc$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits_hc$Chemosynthesis.Obligate <- as.ordered(sp_traits_hc$Chemosynthesis.Obligate) # Chemo Obligate as an ordered trait
sp_traits_hc$Zonation.From.Vent <- ordered(sp_traits_hc$Zonation.From.Vent, levels = c("Low", "Medium", "High")) # Zonation as an ordered trait and reorder levels
sp_traits_hc$Feeding.Mode <- as.factor(sp_traits_hc$Feeding.Mode) # Feeding as a factor trait

# Functional distance between species
func_dist_sp_hc <- funct.dist(sp_traits_hc, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_hc <- quality.fspaces(func_dist_sp_hc, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_hc$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_hc <- func_qual_hc$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_hc <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_hc[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_hc$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities or unique trait combinations (UTC) and derived functional indices
func_ent_hc <- sp.to.fe(sp_traits_hc, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_hc <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_hc) # compute UTC derived indices
func_ind_fe_hc$asb_fdfe # UTC indices

# Summary sensivity indices
func_ind_hc <- data.frame(nb_fe = func_ind_fe_hc$asb_fdfe[,"nb_sp"]) # sp richness
func_ind_hc$fdis <- func_ind_sp_hc$functional_diversity_indices$fdis # functional dispersion 
func_ind_hc$fred <- func_ind_fe_hc$asb_fdfe[,"fred"] # functional redundancy
func_ind_hc$fvuln <- func_ind_fe_hc$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_hc,2), "func_ind_hc.csv")  

# Null models
n_perm = 1000 # number of permutations

multif_perm_hc <- as.list(rep(NA, n_perm)) # list to save permutations results of multidimensional functional alpha-diversity indices (Fric, Fdis)
fe_perm_hc <- as.list(rep(NA, n_perm)) # list to save permutations results of FE indices (FE richness, Fred, Fvul)

for(i in seq(n_perm)){
  
  sp_comm_n_hc = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_hc <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_hc[,c(1:4)],
                                   asb_sp_w = sp_comm_n_hc,
                                   ind_vect = "fdis",
                                   scaling = T,
                                   details_returned = T,
                                   verbose = F) #  Null FDis
  
  multif_perm_hc[[i]] <- multif_n_hc$functional_diversity_indices["fdis"] # store Null FDis results 
  
  func_ind_fe_n_hc <- alpha.fd.fe(sp_comm_n_hc, sp_to_fe = func_ent_hc) # Null UTC indices
  fe_perm_hc[[i]] <- func_ind_fe_n_hc$asb_fdfe[,c("fred", "fvuln")] # store null FRed and FVuln result
  
} # null model of functional indices

# Wrap null model results
multif_perm_a_hc <- array(NA, dim = c(11,1,0)) # create an array with 3 col & 11 rows for null Fric & FDis results  
fe_perm_a_hc <- array(NA, dim = c(11,2,0)) # create an array with 5 col & 11 rows for FE inices results 

for (i in c(1:n_perm)) { 
  multif_perm_a_hc <- abind(multif_perm_a_hc, multif_perm_hc[[i]]) # FDis
  fe_perm_a_hc <- abind(fe_perm_a_hc, fe_perm_hc[[i]]) # FRed and FVuln
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_hc <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_hc), colnames(multif_perm_a_hc))) # matrix to save mean values of FDis permutations
multif_perm_sd_hc <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_hc), colnames(multif_perm_a_hc))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_hc[i,j] = mean(multif_perm_a_hc[i,j,])
    multif_perm_sd_hc[i,j] = sd(multif_perm_a_hc[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_hc <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_hc), colnames(fe_perm_a_hc))) # matrix to save mean values of FRed and FVuln permutations
fe_perm_sd_hc <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_hc), colnames(fe_perm_a_hc))) # matrix to save sd values of FRed and FVuln permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_hc[i,j] = mean(fe_perm_a_hc[i,j,])
    fe_perm_sd_hc[i,j] = sd(fe_perm_a_hc[i,j,])
  }
} # mean and sd values of Fred and FVuln permutations

# Standardized effect size (SES)
ses_multif_hc <- (func_ind_sp_hc$functional_diversity_indices[,"fdis"] - multif_perm_mean_hc) / multif_perm_sd_hc # SES FDis
ses_fe_hc <- (func_ind_fe_hc$asb_fdfe[,c("fred", "fvuln")] - fe_perm_mean_hc) / fe_perm_sd_hc # SES FRed and FVuln

# Pariwise functional turnover null model 
ses_turn_hc <-matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog))))
n_perm = 1000

for (k in rownames(ses_turn_hc)) {
  for (i in rownames(ses_turn_hc[-(which(rownames(ses_turn_hc) == k)),])) {
    sp_comm_sub_hc <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_hc <- sp_faxes_hc[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_hc <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_hc = randomizeMatrix(sp_comm_sub_hc, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_hc <- functional.beta.pair(sp_comm_n_hc, sp_faxes_sub_hc, index.family= "jaccard")
      fb_turn_perm_hc[j, ] <- func_beta_jacc_hc$funct.beta.jtu
    }
    
    fb_turn_obs_hc <- functional.beta.pair(sp_comm_sub_hc, sp_faxes_sub_hc, index.family= "jaccard")
    ses_fb_turn_hc <- (fb_turn_obs_hc$funct.beta.jtu - mean(fb_turn_perm_hc))/ sd(fb_turn_perm_hc)
    ses_turn_hc[k, i] <- ses_fb_turn_hc
  }
} # Functional turnover SES 

#write.csv(ses_turn_hc, "ses_fb_turn_hc.csv") 


# Zonation from vent
####################

sp_traits_zfv <- sfDVent_sp[,c(1,2,3,4,6)] # species traits w/o Zonation from Vent
trait_info <- data.frame(trait_name = colnames(sp_traits_zfv), trait_type = c("O","O","N","O","N")) # Trait info data

sp_traits_zfv$Relative.Adult.Mobility <- as.ordered(sp_traits_zfv$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits_zfv$Estimated.Max.Body.Size <- as.ordered(sp_traits_zfv$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits_zfv$Chemosynthesis.Obligate <- as.ordered(sp_traits_zfv$Chemosynthesis.Obligate) # Chemo Obligate as an ordered trait
sp_traits_zfv$Habitat.Complexity <- as.factor(sp_traits_zfv$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits_zfv$Feeding.Mode <- as.factor(sp_traits_zfv$Feeding.Mode) # Feeding as a factor trait

# Functional distance between species
func_dist_sp_zfv <- funct.dist(sp_traits_zfv, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_zfv <- quality.fspaces(func_dist_sp_zfv, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_zfv$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_zfv <- func_qual_zfv$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_zfv <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_zfv[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_zfv$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities or unique trait combinations (UTC) and derived functional indices
func_ent_zfv <- sp.to.fe(sp_traits_zfv, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_zfv <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_zfv) # compute UTC derived indices
func_ind_fe_zfv$asb_fdfe # UTC indices

# Summary sensivity indices
func_ind_zfv <- data.frame(nb_fe = func_ind_fe_zfv$asb_fdfe[,"nb_sp"]) # sp richness
func_ind_zfv$fdis <- func_ind_sp_zfv$functional_diversity_indices$fdis # functional dispersion 
func_ind_zfv$fred <- func_ind_fe_zfv$asb_fdfe[,"fred"] # functional redundancy
func_ind_zfv$fvuln <- func_ind_fe_zfv$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_zfv,2), "func_ind_zfv.csv")  

# Null models
n_perm = 1000 # number of permutations

multif_perm_zfv <- as.list(rep(NA, n_perm)) # list to save permutations results of Fdis
fe_perm_zfv <- as.list(rep(NA, n_perm)) # list to save permutations results of Fred and Fvul

for(i in seq(n_perm)){
  
  sp_comm_n_zfv = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_zfv <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_zfv[,c(1:4)],
                                    asb_sp_w = sp_comm_n_zfv,
                                    ind_vect = "fdis",
                                    scaling = T,
                                    details_returned = T,
                                    verbose = F) #  Null FDis
  
  multif_perm_zfv[[i]] <- multif_n_zfv$functional_diversity_indices["fdis"] # store Null FDis results 
  
  func_ind_fe_n_zfv <- alpha.fd.fe(sp_comm_n_zfv, sp_to_fe = func_ent_zfv) # Null UTC indices
  fe_perm_zfv[[i]] <- func_ind_fe_n_zfv$asb_fdfe[,c("fred", "fvuln")] # store null Fred and FVuln result
  
} # null model of functional indices

# Wrap null model results
multif_perm_a_zfv <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for null Fric & FDis results  
fe_perm_a_zfv <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for FE inices results 

for (i in c(1:n_perm)) { 
  multif_perm_a_zfv <- abind(multif_perm_a_zfv, multif_perm_zfv[[i]]) # FRic & FDis
  fe_perm_a_zfv <- abind(fe_perm_a_zfv, fe_perm_zfv[[i]]) # FRed and FVuln indices
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_zfv <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_zfv), colnames(multif_perm_a_zfv))) # matrix to save mean values of FDis permutations
multif_perm_sd_zfv <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_zfv), colnames(multif_perm_a_zfv))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_zfv[i,j] = mean(multif_perm_a_zfv[i,j,])
    multif_perm_sd_zfv[i,j] = sd(multif_perm_a_zfv[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_zfv <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_zfv), colnames(fe_perm_a_zfv))) # matrix to save mean values of FRed and FVuln permutations
fe_perm_sd_zfv <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_zfv), colnames(fe_perm_a_zfv))) # matrix to save sd values of FRed and FVuln permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_zfv[i,j] = mean(fe_perm_a_zfv[i,j,])
    fe_perm_sd_zfv[i,j] = sd(fe_perm_a_zfv[i,j,])
  }
} # mean and sd values of FRed and FVuln

# Standardized effect size (SES)
ses_multif_zfv <- (func_ind_sp_zfv$functional_diversity_indices[,"fdis"] - multif_perm_mean_zfv) / multif_perm_sd_zfv # SES FDis
ses_fe_zfv <- (func_ind_fe_zfv$asb_fdfe[,c("fred", "fvuln")] - fe_perm_mean_zfv) / fe_perm_sd_zfv # SES FRed and FVuln

# Pariwise functional turnover null model 
ses_turn_zfv <-matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog))))
n_perm = 1000

for (k in rownames(ses_turn_zfv)) {
  for (i in rownames(ses_turn_zfv[-(which(rownames(ses_turn_zfv) == k)),])) {
    sp_comm_sub_zfv <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_zfv <- sp_faxes_zfv[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_zfv <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_zfv = randomizeMatrix(sp_comm_sub_zfv, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_zfv <- functional.beta.pair(sp_comm_n_zfv, sp_faxes_sub_zfv, index.family= "jaccard")
      fb_turn_perm_zfv[j, ] <- func_beta_jacc_zfv$funct.beta.jtu
    }
    
    fb_turn_obs_zfv <- functional.beta.pair(sp_comm_sub_zfv, sp_faxes_sub_zfv, index.family= "jaccard")
    ses_fb_turn_zfv <- (fb_turn_obs_zfv$funct.beta.jtu - mean(fb_turn_perm_zfv))/ sd(fb_turn_perm_zfv)
    ses_turn_zfv[k, i] <- ses_fb_turn_zfv
  }
} # Functional turnover SES 

#write.csv(ses_turn_zfv, "ses_fb_turn_zfv.csv") 


# Chemosynthetic obligate
#########################

sp_traits_chemo <- sfDVent_sp[,c(1,2,3,5,6)] # species traits w/o Chemo Obligate
trait_info <- data.frame(trait_name = colnames(sp_traits_chemo), trait_type = c("O","O","N","O","N")) # Trait info data

sp_traits_chemo$Relative.Adult.Mobility <- as.ordered(sp_traits_chemo$Relative.Adult.Mobility) # Mobility as an ordered trait
sp_traits_chemo$Estimated.Max.Body.Size <- as.ordered(sp_traits_chemo$Estimated.Max.Body.Size) # Size as an ordered trait
sp_traits_chemo$Zonation.From.Vent <- ordered(sp_traits_chemo$Zonation.From.Vent, levels = c("Low", "Medium", "High")) # Zonation as an ordered trait and reorder levels
sp_traits_chemo$Habitat.Complexity <- as.factor(sp_traits_chemo$Habitat.Complexity) # Habitat Complexity as a factor trait
sp_traits_chemo$Feeding.Mode <- as.factor(sp_traits_chemo$Feeding.Mode) # Feeding as a factor trait

# Functional distance between species
func_dist_sp_chemo <- funct.dist(sp_traits_chemo, tr_cat = trait_info, metric = "gower", ordinal_var = "metric", weight_type = "equal") # Functional dissimilarity (Gower) distance between species 

# Functional space quality
func_qual_chemo <- quality.fspaces(func_dist_sp_chemo, deviation_weighting = c("absolute", "squarred"),maxdim_pcoa = 9, fdist_scaling = c(TRUE, FALSE)) # estimating the quality of the functional space

apply(func_qual_chemo$quality_fspaces, 2, which.min) # best number of dimesions calculated with each index (4D) 
sp_faxes_chemo <- func_qual_chemo$details_fspaces$sp_pc_coord # species coordinates in the functional space

# Multidimensional alpha-diversity indices
func_ind_sp_chemo <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_chemo[,c(1:4)], asb_sp_w = sp_comm_biog, scaling = TRUE, details_returned = T, verbose = T) # compute mutidimensinal functional diversity indices in a 4D functional space
func_ind_sp_chemo$functional_diversity_indices # mutidimensinal indices of functional diversity

# Functional Entities or unique trait combinations (UTC) and derived functional indices
func_ent_chemo <- sp.to.fe(sp_traits_chemo, trait_info, fe_nm_type = "fe_rank") # compute UTC
func_ind_fe_chemo <- alpha.fd.fe(sp_comm_biog, sp_to_fe = func_ent_chemo) # compute UTC derived indices
func_ind_fe_chemo$asb_fdfe # UTC indices

# Summary sensivity indices
func_ind_chemo <- data.frame(nb_fe = func_ind_fe_chemo$asb_fdfe[,"nb_sp"]) # sp richness
func_ind_chemo$fdis <- func_ind_sp_chemo$functional_diversity_indices$fdis # functional dispersion 
func_ind_chemo$fred <- func_ind_fe_chemo$asb_fdfe[,"fred"] # functional redundancy
func_ind_chemo$fvuln <- func_ind_fe_chemo$asb_fdfe[,"fvuln"] # functional vulnerability
#write.csv(round(func_ind_chemo,2), "func_ind_chemo.csv")  

# Null models
n_perm = 1000 # number of permutations

multif_perm_chemo <- as.list(rep(NA, n_perm)) # list to save permutations results of Fdis
fe_perm_chemo <- as.list(rep(NA, n_perm)) # list to save permutations results of Fred and Fvul

for(i in seq(n_perm)){
  
  sp_comm_n_chemo = randomizeMatrix(sp_comm_biog, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
  
  multif_n_chemo <- alpha.fd.multidim(sp_faxes_coord = sp_faxes_chemo[,c(1:4)],
                                      asb_sp_w = sp_comm_n_chemo,
                                      ind_vect = "fdis",
                                      scaling = T,
                                      details_returned = T,
                                      verbose = F) #  Null Fric & FDis
  
  multif_perm_chemo[[i]] <- multif_n_chemo$functional_diversity_indices["fdis"] # store Nnull Fric & FDis results 
  
  func_ind_fe_n_chemo <- alpha.fd.fe(sp_comm_n_chemo, sp_to_fe = func_ent_chemo) # Null FRed and FVuln
  fe_perm_chemo[[i]] <- func_ind_fe_n_chemo$asb_fdfe[,c("fred", "fvuln")] # store null FRed and FVuln result
  
} # null model of functional indices

# Wrap null model results
multif_perm_a_chemo <- array(NA, dim = c(11,1,0)) # create an array with 1 col & 11 rows for null FDis results  
fe_perm_a_chemo <- array(NA, dim = c(11,2,0)) # create an array with 2 col & 11 rows for FRed and FVuln results 

for (i in c(1:n_perm)) { 
  multif_perm_a_chemo <- abind(multif_perm_a_chemo, multif_perm_chemo[[i]]) # FDis
  fe_perm_a_chemo <- abind(fe_perm_a_chemo, fe_perm_chemo[[i]]) # FRed and FVuln
} # merge each elemnt of the lists as dimensions of the arrays

multif_perm_mean_chemo <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_chemo), colnames(multif_perm_a_chemo))) # matrix to save mean values of FDis permutations
multif_perm_sd_chemo <- matrix(NA, ncol = 1, nrow= 11, dimnames = list(rownames(multif_perm_a_chemo), colnames(multif_perm_a_chemo))) # matrix to save sd values of FDis permutations 
for (i in 1:11) {
  for (j in 1:1) {
    multif_perm_mean_chemo[i,j] = mean(multif_perm_a_chemo[i,j,])
    multif_perm_sd_chemo[i,j] = sd(multif_perm_a_chemo[i,j,])
  }
} # mean and sd values of FDis permutations 

fe_perm_mean_chemo <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_chemo), colnames(fe_perm_a_chemo))) # matrix to save mean values of FRed and FVuln permutations
fe_perm_sd_chemo <- matrix(NA, ncol = 2, nrow= 11, dimnames = list(rownames(fe_perm_a_chemo), colnames(fe_perm_a_chemo))) # matrix to save sd values of FRed and FVuln permutations 
for (i in 1:11) {
  for (j in 1:2) {
    fe_perm_mean_chemo[i,j] = mean(fe_perm_a_chemo[i,j,])
    fe_perm_sd_chemo[i,j] = sd(fe_perm_a_chemo[i,j,])
  }
} # mean and sd values of FRed and FVuln permutations

# Standardized effect size (SES)
ses_multif_chemo <- (func_ind_sp_chemo$functional_diversity_indices["fdis"] - multif_perm_mean_chemo) / multif_perm_sd_chemo # SES FDis
ses_fe_chemo <- (func_ind_fe_chemo$asb_fdfe[,c("fred", "fvuln")] - fe_perm_mean_chemo) / fe_perm_sd_chemo # SES FRed and Fvuln

# Pariwise functional turnover null model 
ses_turn_chemo <-matrix(NA, nrow= length(rownames(sp_comm_biog)), ncol = length(rownames(sp_comm_biog)), dimnames = list(c(rownames(sp_comm_biog)), c(rownames(sp_comm_biog))))
n_perm = 1000

for (k in rownames(ses_turn_chemo)) {
  for (i in rownames(ses_turn_chemo[-(which(rownames(ses_turn_chemo) == k)),])) {
    sp_comm_sub_chemo <- sp_comm_biog[c(k, i), -c(which(colSums(sp_comm_biog[c(k, i),]) == 0))] # sp comm
    sp_faxes_sub_chemo <- sp_faxes_chemo[-c(which(colSums(sp_comm_biog[c(k, i),]) == 0)), c(1:4)]
    fb_turn_perm_chemo <- matrix(NA, nrow= n_perm, ncol = 1, dimnames = list(c(seq(n_perm)), "turn_jacc"))
    
    for (j in seq(n_perm)){
      sp_comm_n_chemo = randomizeMatrix(sp_comm_sub_chemo, null.model = "independentswap") # Randomizations maintaining species occurrence frequency (global) and region species richness  
      func_beta_jacc_chemo <- functional.beta.pair(sp_comm_n_chemo, sp_faxes_sub_chemo, index.family= "jaccard")
      fb_turn_perm_chemo[j, ] <- func_beta_jacc_chemo$funct.beta.jtu
    }
    
    fb_turn_obs_chemo <- functional.beta.pair(sp_comm_sub_chemo, sp_faxes_sub_chemo, index.family= "jaccard")
    ses_fb_turn_chemo <- (fb_turn_obs_chemo$funct.beta.jtu - mean(fb_turn_perm_chemo))/ sd(fb_turn_perm_chemo)
    ses_turn_chemo[k, i] <- ses_fb_turn_chemo
  }
}# Functional turnover SES 

#write.csv(ses_turn_chemo, "ses_fb_turn_chemo.csv") 

#############
## THE END ##
#############

