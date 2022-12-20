
sim_study <- function(d = 5,p=0.2, rng = NULL){
  ## perform a simulation study to measure performance of the second step algorithm.
  ##
  ## Args:
  ##    - d: dimension.
  ##    - tol: tolerance for algorithm.
  ## Returns:
  ## a tibble with time, KKT condition values, duality gap, the dual dimension 
  ## and the number of failed traingle inequalities in the first step
  
  # check arguments
  
  
  # set seed, if applicable
  if (!is.null(rng)){
    rng_sims <- rng[[1]]
    rngtools::setRNG(rng_sims)
  }
  
  # perform simulation

  time <- numeric(1)
  KKT <- numeric(3)
  ineq <- numeric(1)
  dual_gap <- numeric(1)
  start_prop <- numeric(1)
  
  ## Generate random Gamma and graph
  Sig=matrix(.99,d-1,d-1)
  diag(Sig)=rep(1,d-1)
  Y <-  mvrnorm(n=d,mu=rep(0,d-1),Sigma = Sig)
  R <-  diag(1/sqrt(apply(Y^2 , 1, sum)))%*%Y
  S <-  R%*%t(R)
  par <-  Sigma2Gamma(S)
    tri=0
    lm_check = TRUE
    while (tri==0|lm_check==TRUE) {
      graph = graphicalExtremes:::generate_random_connected_graph(d=d,p=p)
      tri=length(triangles(graph = graph))
      lm = loc_metr(par,graph)
      lm_check=lm$lm
    }
  ## Compute step 1 estimate
  par = complete_Gamma(par,graph,final_tol=1e-6,N=50000)
  
  ##Compute step 2 estimate
  ptm <- proc.time()[1]
  res = loc_metr_alg(Gamma = par, graph = graph)
  time <- proc.time()[1] - ptm
  
  G = res$Gamma_opt
  Theta <- Gamma2Theta(G)
  ineq <- length(res$mu_opt)
  
  start_prop=1-length(lm$failed_ineq)/ineq
  AtimesGvec = tAtimesGamma(G,graph)
  
  KKT[1] <- max(AtimesGvec)
  KKT[2] <- min(res$mu_opt)
  KKT[3] <- abs(sum(res$mu_opt * AtimesGvec))
  
  dual_gap <- 1/2*sum(-Gamma2Theta(par)*G) - (d-1) 
  
  tbl <- tibble(type = paste0("time"), 
              value = time) %>% 
    bind_rows(tibble(type = paste0("KKT", 1:length(KKT)), value =  KKT))%>%
    bind_rows(tibble(type = paste0("dual_gap"), value =  dual_gap))%>%
    bind_rows(tibble(type = paste0("ineq"), value = ineq ))%>%
    bind_rows(tibble(type = paste0("start_prop"), value = start_prop ))
  
  return(tbl)
}


wrapper_sim <- function(i, rowid, sim_fn, sim_fn_args){
  ## apply arguments sim_fn_args[i] to sim_fn
  ## Args:
  ##    - i (integer): row to consider from sim_fn_args
  ##    - rowid (integer): unique identifier of the current simulation row
  ##      (not necessarily equal to i)
  ##    - sim_fn (function): function to run
  ##    - sim_fn_args (tibble): tibble with arguments to pass to sim_fn
  ##
  ## Returns:
  ##    - tibble with simulation results
  
  do.call(what = sim_study, args = fun_args[i, ]) %>% 
    mutate(rowid = rowid)
}

set_rng <- function(tbl, seed){
  ## adds to tbl a column with seeds to generate independent streams of random
  ## numbers.
  ##
  ## Args:
  ##     - tbl: a tibble where the columns contain the parameter settings and the
  ##       rows contain the simulation runs.
  ##
  ## Returns:
  ##     The function returns tbl appending a column with seeds used to generate
  ##     independent streams of random numbers.
  ##
  ## Note:
  ##     This function creates ensures that the simulations are fully repeatable.
  ##     This is possible because it assigns to each simulation run a unique 
  ##     random seed (generated with L'Ecuyer RNG method, which is suitable 
  ##     for parallel processing, too).
  
  m <- n_groups(tbl)
  group_idxs <- group_indices(tbl)
  
  # create independent RNG streams with L'Ecuyer method
  rng <- RNGseq(m, seed = seed, simplify = FALSE)
  rng <- rng[group_idxs]
  
  # add RNG streams to tbl
  tbl$rng <- rng
  
  # return tibble
  return(tbl)
}

rep_tibble <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    select(-rowname)
  
} 

assign_random_seed <- function(tbl, grouping_vars, seed){
  ## tibble character_vector integer -> tibble
  ## assign random seed according to the variables in grouping_vars
  if (is.null(grouping_vars)){
    tbl <- tbl %>%
      rowwise()
  } else {
    tbl <- tbl %>% 
      group_by(across(all_of(grouping_vars)))
  }
  
  tbl %>% 
    set_rng(seed) %>% 
    ungroup()
}



#rep_tibble_new solves an issue with package intersection

rep_tibble_new <- function(tbl, m){
  ## tibble integer -> tibble
  ## replicate tibble m times and append named column rep_id = 1:m
  
  tbl <-  tbl %>% 
    rownames_to_column()
  
  expand_grid(rep_id = 1:m,
              rowname = tbl$rowname) %>% 
    left_join(tbl, by = "rowname") %>% 
    dplyr::select(-rowname)
  
}



