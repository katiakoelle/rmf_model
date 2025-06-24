


library(ggplot2)
library(dplyr)
library(tidyverse)
library(patchwork)
library(tictoc)


## PARAMETERS

params.L_P <- 267 # number of nonsynonymous phenotypic sites
params.L_PA <- 0 # number of nonsynonymous pleiotropic sites
params.L_A <- 48 # number of nonsynonymous antigenic sites
params.L_S <- 85 # number of synonymous sites
params.L_total <- params.L_P+params.L_PA+params.L_A+params.L_S

params.N <- 5000 # effective population size

proportion_maladapted <- 0.5 # proportion of sites that contain 0 in the starting genotype
params.dist_away_start <- round(proportion_maladapted*params.L_total)

params.k <- 100 # controls the maximum number to linearly transform the data

params.c <- 0.2 # controls the degree of landscape ruggedness

params.q <- 0.5 # strength of immune pressure
params.p <- 0.4 # breadth of immunity

params.mu_rate_perSite_perRepCycle <- 2.5E-5 
params.mut_prob_perRepCycle <- 1-exp(-(params.mu_rate_perSite_perRepCycle*params.L_total)) # probability of one mutation per replication cycle

params.d <- 4 # deaths per day

params.deltaT <- 1/24 # 1 hour
params.tyrs_end <- 1 # Number of Years of Infection 



## Functions Required For Simulation and Data Analysis

GetGenotypeFitnessP <- function(this_genotype) {
  D = sum(params.refSequence != this_genotype)
  eta_sigma <- rnorm(1, mean=0, sd=1)  
  
  phi <- -params.c*D + eta_sigma
  
  min_x = -params.c * (params.L_P+params.L_PA)
  max_x = 0
  
  min_y = 0
  max_y = params.k
  
  lr_b = max_y
  lr_m = (max_y - min_y) / (max_x - min_x)
  
  F_P = lr_m * phi + lr_b
  
  if(F_P < 0) {
    F_P = 0
  }
  
  return(F_P)
}

GetGenotypeFitnessA <- function(consensus_genotype, this_genotype) {
  D_A <- sum(consensus_genotype != this_genotype)
  F_A <- 1-params.q*(params.p^D_A)
  
  return(F_A)
}

CalculateNewConsensusSequence <- function(genotypes_antigenic_region, curr_N_by_genotype) {
  sumBySite <- colSums(curr_N_by_genotype*genotypes_antigenic_region)
  fractionBySite <- sumBySite/params.N
  new_consensus_genotype <- round(fractionBySite)
  
  return(new_consensus_genotype)
}


GetOffspringGenotype <- function(parent_genotype, parent_fitness_P, parent_fitnes_A, consensus_genotype) {
  
  genome_loc <- sample(1:params.L_total, 1)
  offspring_genotype <- parent_genotype
  
  if (offspring_genotype[genome_loc] == 0) {
    offspring_genotype[genome_loc] <- 1
  } else {
    offspring_genotype[genome_loc] <- 0
  }
  
  nGenotypes <- length(parentVector)
  
  if (TRUE) {
    n_ones_offspring_genotype <- sum(offspring_genotype)
    n_ones_genotypeVector <- rowSums(genotypeVector)
    locs_possible_identical <- which(n_ones_genotypeVector == n_ones_offspring_genotype)
    
    if (length(locs_possible_identical) != 0) {
      
      for (i in 1:length(locs_possible_identical)) {
        ndiff <- sum(offspring_genotype != genotypeVector[locs_possible_identical[i],])
        if (ndiff == 0) {
          boolNewGenotype <- 0
          locOldGenotype <- locs_possible_identical[i]
          
          offspring_genotype <- NaN
          offspringFitness_P <- NaN
          offspringFitness_A <- NaN
          offspringFitness <- NaN
          
          if (genome_loc > params.L_P & genome_loc <= (params.L_P+params.L_PA+params.L_A)) {
            boolAntigenicChange <- 1
          } else (
            boolAntigenicChange <- 0
          )
          return(list(boolNewGenotype = boolNewGenotype, locOldGenotype = locOldGenotype, offspring_genotype = offspring_genotype, offspringFitness_P = offspringFitness_P, offspringFitness_A = offspringFitness_A, offspringFitness = offspringFitness, boolAntigenicChange = boolAntigenicChange))
          
          }
      }
    }
  }  
  
  boolNewGenotype <- 1
  locOldGenotype <- NaN
  
  if (genome_loc <= params.L_P) {
    offspringFitness_P <- GetGenotypeFitnessP(offspring_genotype[1:(params.L_P+params.L_PA)])
    offspringFitness_A <- parent_fitnes_A
    boolAntigenicChange <- 0
  } else if (genome_loc > params.L_P & genome_loc <= (params.L_P+params.L_PA)) {
    offspringFitness_P <- GetGenotypeFitnessP(offspring_genotype[1:(params.L_P+params.L_PA)])
    offspringFitness_A <- GetGenotypeFitnessA(consensus_genotype, offspring_genotype[(params.L_P+1):(params.L_P+params.L_PA+params.L_A)]) 
    boolAntigenicChange <- 1
  } else if (genome_loc > (params.L_P+params.L_PA) & genome_loc <= (params.L_P+params.L_PA+params.L_A)) {
    offspringFitness_P <- parent_fitness_P
    offspringFitness_A <- GetGenotypeFitnessA(consensus_genotype, offspring_genotype[(params.L_P+1):(params.L_P+params.L_PA+params.L_A)])
    boolAntigenicChange <- 1
  } else {
    offspringFitness_P <- parent_fitness_P
    offspringFitness_A <- parent_fitnes_A
    boolAntigenicChange <- 0
  }
  
  offspringFitness <- offspringFitness_P*offspringFitness_A
  
  return(list(boolNewGenotype = boolNewGenotype, locOldGenotype = locOldGenotype, offspring_genotype = offspring_genotype, offspringFitness_P = offspringFitness_P, offspringFitness_A = offspringFitness_A, offspringFitness = offspringFitness, boolAntigenicChange = boolAntigenicChange))
}

GetFitnessPData <- function() {
  t_yrs <- tList/365.25
  
  nGenotypes <- nrow(nListByGenotype)
  n_timePts <- ncol(nListByGenotype)
  
  mean_fitnessP <- numeric(n_timePts) 
  
  for (i in 1:n_timePts) {
    mean_fitnessP[i] <- sum((genotypeFitness_P*nListByGenotype[ ,i])/params.N)
  }
  
  return(mean_fitnessP)
}



GetDivergenceData <- function() {
  
  nGenotypes <- nrow(nListByGenotype)
  n_timePts <- ncol(nListByGenotype)
  
  divergence <- numeric(n_timePts) 
  dist <- c()
  for (i in 1:nrow(genotypeVector)) {
    dist[i] <- sum(params$initial_genotype != genotypeVector[i,])
  }
  
  for (i in 1:n_timePts) {

    divergence[i] <- sum(dist*nListByGenotype[ ,i])/params.N
  }
  
  
  return(divergence)
}



GetAntigenicDivergenceData <- function() {
  
  nGenotypes <- nrow(nListByGenotype)
  n_timePts <- ncol(nListByGenotype)
  
  divergence <- numeric(n_timePts) 
  dist <- c()
  for (i in 1:nrow(genotypeVector)) {
    dist[i] <- sum(params$initial_genotype[(params.L_P+1):(params.L_P+params.L_PA+params.L_A)] != genotypeVector[i, (params.L_P+1):(params.L_P+params.L_PA+params.L_A)])
  }
  
  for (i in 1:n_timePts) {

    divergence[i] <- sum(dist*nListByGenotype[ ,i])/params.N
  }
  
  return(divergence)
}



## Setting up starting model framework



# Start with a genotype that has an overall hamming distance of params.dist_away_start nucleotides
# away from the reference.

init_genotype_temp <- rep(1, params.L_total) #start off with P sites mostly good. (so set to 1)
locs <- sample(1:params.L_total, params.dist_away_start, replace=FALSE) 

params.init_genotype <- init_genotype_temp
params.init_genotype[locs] <- 0

genotypeVector <- matrix(nrow=1, ncol=params.L_total) # this will be nrows = (# explored genotypes), ncols = L
genotypeVector[1,] <- params.init_genotype # add starting gneoytpe to the first row


parentVector <- NaN
params.refSequence <- rep(1, (params.L_P+params.L_PA)) # (highest fitness genotype - which is all 1's)


# only the part of the genome that codes for P and PA gets handed in to calculate the phenotypic fitness of a genotype:
# (organizing the genome as: P, PA, A, S)

genotypeFitness_P <- c()
genotypeFitness_P[1] <- GetGenotypeFitnessP(genotypeVector[1, (1:(params.L_P+params.L_PA))])



### Adding Antigenicity
if (params.L_PA+params.L_A > 0) {
  consensus_genotype <- genotypeVector[1, (params.L_P+1):(params.L_P+params.L_PA+params.L_A)]
} else if (params.L_PA+params.L_A == 0) {
  consensus_genotype <- c()
} else (stop("params.L_PA + params.L_A is less than zero"))

genotypeFitness_A <- c()
genotypeFitness_A <- GetGenotypeFitnessA(consensus_genotype, genotypeVector[1, (params.L_P+1):(params.L_P+params.L_PA+params.L_A)])


genotypeFitness <- c(genotypeFitness_P * genotypeFitness_A)


params <- list(L_P = params.L_P, L_PA = params.L_PA, L_A = params.L_A, L_S = params.L_S, L = params.L_total, N = params.N, dist_away_start = params.dist_away_start, k = params.k, c = params.c, q=params.q, p=params.p, mu_perSite_perRepCycle = params.mu_rate_perSite_perRepCycle, d = params.d, deltaT = params.deltaT, tyrs_end = params.tyrs_end, t_end = params.t_end, mut_prob_per_rep = params.mut_prob_perRepCycle, initial_genotype = params.init_genotype, init_parent = parentVector, ref_seq = params.refSequence)


tic()

tList <- c(0)

curr_N_by_genotype <- params.N

nListByGenotype <- matrix(nrow=length(curr_N_by_genotype), ncol=1) # rows equate to the number of genotypes and columns equate to the time points
nListByGenotype[,1] <- curr_N_by_genotype

cntr <- 1

n_repeated_draws_tau_leap <- 0



## Run Model

for (t in seq(from = 0, to = params.t_end, by = params.deltaT)) {
  
  
  
  while (TRUE) {
    
    #determining the number of birth/death events in the deltaT interval:
    temp_curr_N_by_genotype <- curr_N_by_genotype
    
    n_deaths_by_genotype <- rpois(length(temp_curr_N_by_genotype), (params.d*temp_curr_N_by_genotype*params.deltaT))
    n_births <- sum(n_deaths_by_genotype)
    
    
    prob_birth <- (genotypeFitness*curr_N_by_genotype) / sum(genotypeFitness*curr_N_by_genotype)
    n_births_by_genotype <- as.numeric(rmultinom(1, n_births, prob_birth))
    
    # figure out which of the parental genotypes will have mutations:
    mutations <- rbinom(length(n_births_by_genotype), n_births_by_genotype, params.mut_prob_perRepCycle)
    
    
    # update population size vector by removing deaths and adding in non-mutated births:
    temp_curr_N_by_genotype <- temp_curr_N_by_genotype - n_deaths_by_genotype + (n_births_by_genotype - mutations)
    
    # check to make sure tau leap resulted in no negative numbers of genotypes due to too large step size:
    if (min(temp_curr_N_by_genotype) >=0) {
      curr_N_by_genotype <- temp_curr_N_by_genotype
      
      break
    }
    
    # otherwise, do this again!
    n_repeated_draws_tau_leap <- n_repeated_draws_tau_leap + 1
    
  }
  
  
  
  
  ### Adding MUTATED GENOTYPES
  
  
  n_existing_genotypes <- length(curr_N_by_genotype)
  
  antigenicChange <- 0
  
  if (sum(mutations) > 0) {
    
    for (i in 1:n_existing_genotypes) {
      n_mut_from_this_genotype <- mutations[i]
      if (n_mut_from_this_genotype > 0) {
        for (j in 1:n_mut_from_this_genotype) {
          bool_results <- GetOffspringGenotype(genotypeVector[i, ], genotypeFitness_P[i], genotypeFitness_A[i], consensus_genotype)
           if (!bool_results$boolNewGenotype) {
            curr_N_by_genotype[bool_results$locOldGenotype] <- curr_N_by_genotype[bool_results$locOldGenotype] + 1
          } else {
            genotypeVector <- rbind(genotypeVector, bool_results$offspring_genotype)
            genotypeFitness_P <- c(genotypeFitness_P, bool_results$offspringFitness_P)
            genotypeFitness_A <- c(genotypeFitness_A, bool_results$offspringFitness_A)
            genotypeFitness <- c(genotypeFitness, bool_results$offspringFitness)
            parentVector <- c(parentVector, i)
            curr_N_by_genotype <- c(curr_N_by_genotype, 1)
          }
          antigenicChange <- antigenicChange + bool_results$boolAntigenicChange
          
        }
      }
    }
  }
  
  if (antigenicChange > 0) { # then there was at least one mutation that occurred in either PA or A sites.
    new_consensus_genotype <- CalculateNewConsensusSequence(genotypeVector[, (params.L_P+1):(params.L_P+params.L_PA+params.L_A)], curr_N_by_genotype)
    
    if (sum(new_consensus_genotype != consensus_genotype) != 0) {# then consensus has changed
      # then update genotypeFitness_A and genotypeFitness vectors
      
      for (this_g in 1:length(curr_N_by_genotype)) {
        genotypeFitness_A[this_g] <- GetGenotypeFitnessA(new_consensus_genotype, genotypeVector[this_g, (params.L_P+1):(params.L_P+params.L_PA+params.L_A)])
      }
      genotypeFitness <- genotypeFitness_P * genotypeFitness_A
      consensus_genotype <- new_consensus_genotype
    }
  }
  
  
  ### saving data after mutations are added - updating results dataframe every 24 hours and saving data every 1 year
  cntr <- cntr + 1
  
  if (cntr %% 24 == 0) {
    tList <- c(tList, t)
    
    nrows_to_add <- length(curr_N_by_genotype) - nrow(nListByGenotype)
    zeros_to_add <- matrix(0, nrow = nrows_to_add, ncol = ncol(nListByGenotype))
    nListByGenotype <- rbind(nListByGenotype, zeros_to_add)
    nListByGenotype <- cbind(nListByGenotype, curr_N_by_genotype)
    
  }
  
  years_simulated <- t/365.25
    print(years_simulated) # print the current time
  
  
}

### Updating the results data after simulation has completed
tList <- c(tList, t)
nrows_to_add <- length(curr_N_by_genotype) - nrow(nListByGenotype)
zeros_to_add <- matrix(0, nrow = nrows_to_add, ncol = ncol(nListByGenotype))
nListByGenotype <- rbind(nListByGenotype, zeros_to_add)
nListByGenotype <- cbind(nListByGenotype, curr_N_by_genotype)

mean_fitnessP <- GetFitnessPData()
divergence <- GetDivergenceData()
Antigenic_divergence <- GetAntigenicDivergenceData()

warnings()


toc()



results <- list(tList=tList, genotypeVector=genotypeVector, parentVector=parentVector, nListByGenotype=nListByGenotype, genotypeFitness=genotypeFitness, genotypeFitness_P=genotypeFitness_P, genotypeFitness_A=genotypeFitness_A, mean_fitnessP=mean_fitnessP, divergence=divergence, Antigenic_divergence=Antigenic_divergence)


saveRDS(results, "results.RData")
saveRDS(params, "params.RData")






