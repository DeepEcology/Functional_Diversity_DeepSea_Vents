
###############################################################################################
### Script of "High functional vulnerability across deep-sea hydrothermal vent communities" ###
### Function to simulate species extintion effects on functional diversity                  ###
### Joan M. Alfaro (jmalfarolucas@gmail.com)                                                ###
### Department of Biology, University of Victoria (Canada)                                  ###
### November 2023                                                                           ### 
###############################################################################################

extinct.sim <- function(comm, nperm = 10) { # comm = matrix/data.frame with communities in rows and species in columns; nperm = number of extintions simulations per community at each step   
  
  extinct_comm <- as.list(rep(NA, length(rownames(comm)))) # list for store simulated extinctions community compositions for each community 
  names(extinct_comm) = rownames(comm)
  summ_comm <- asb.sp.summary(comm) # summary of communities to retrive species richness  
  
  extinct_ind <- as.list(rep(NA, length(rownames(comm)))) # list for store the results of indices derived from extinctions simulations  
  names(extinct_ind) = rownames(comm)
  
  for (k in c(rownames(comm))) {
    
    comm_k <- summ_comm$asb_sp_nm[[k]] # retrieve commuinity k composition  
    sp <- c((length(comm_k)-1):5) # total species richness of each extinction step 
    extinct_comm[[k]] <- as.list(rep(NA, nperm)) # list to store the simulated communities after each extinction step 
    ind_ext_perm <- data.frame(comm_name = NA, perm = NA, sp_ext = NA, FRic = NA, FDis = NA)
    ind_ext_perm <- ind_ext_perm[-1,]
    
    for (j in 1:nperm) {
      
      comm_ext_perm <- as.list(rep(NA, length(sp))) # list to store each simulation extinction community
      ind_ext_perm_j <- data.frame(comm_name = NA, perm = NA, sp_ext = NA, FRic = NA, FDis = NA)
      ind_ext_perm_j <- ind_ext_perm_j[-1,]
      
      for (i in 1:length(sp)) {
        
        sp_left <- sample(x = comm_k, size = sp[i]) # simulate species left after each extintion per each community  
        comm_ext_perm[[i]] <-  matrix(rep(1, sp[i]), ncol = sp[i], nrow = 1, dimnames = list(k, names(sp_left))) # save extintions results
        
        multi <- matrix(sp_faxes[c(names(sp_left)),c(1:4)], nrow = length(names(sp_left)), ncol = 4, 
                        dimnames = list(c(names(sp_left)), c("PC1", "PC2", "PC3", "PC4")))
        
        FRic_ext <- alpha.fd.multidim(sp_faxes_coord = multi,
                                      asb_sp_w = comm_ext_perm[[i]],
                                      ind_vect = c("fric", "fdis"),
                                      scaling = F,
                                      details_returned = F,
                                      verbose = F)
        
        ind_ext_perm_i <- data.frame(comm_name = k, perm = j, sp_ext = summ_comm$asb_sp_richn[k]-sp[i], 
                                     FRic = FRic_ext$functional_diversity_indices$fric, FDis = FRic_ext$functional_diversity_indices$fdis)
        
        ind_ext_perm_j <- rbind(ind_ext_perm_j, ind_ext_perm_i)
        
      }  
      
      extinct_comm[[k]][[j]] <- comm_ext_perm
      ind_ext_perm <- rbind(ind_ext_perm, ind_ext_perm_j)
      
    }
    
    extinct_ind[[k]] <- ind_ext_perm
    
  }
  
  extinct <- as.list(rep(NA, 2)) # list to store community extinctions simulations and indices results
  names(extinct) <- c("extinct_comm", "extinct_ind")
  
  extinct[[1]] <- extinct_comm
  extinct[[2]] <- extinct_ind
  
  return(extinct)
} # Function fro simulating species extinction on functional metrics

###############
### The end ###
###############

