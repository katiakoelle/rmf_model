# Immune Pressure is Key to Understanding Observed Patterns of Respiratory Virus Evolution in Prolonged Infections

This repository contains the code that we used to simulate within-host evolution and analyze the evolutionary patterns. 

> Immune Pressure is Key to Understanding Observed Patterns of Respiratory Virus Evolution in Prolonged Infections
Amber Coats, Yintong Rita Wang, Katia Koelle
bioRxiv 2025.04.30; doi: https://doi.org/10.1101/2025.04.28.651059


## Accessing Model Code
The models that were developed here can be found in the directories: 

1. RMF_New_Landscape: contains the code for a single stochastic realization. Running this code twice with the same parameter settings will result in two independent replicates of virus evolution on different fitness landscapes that were generated from the same parameter settings. 
2. RMF_Same_Landscape: contains the code for running multiple stochastic realizations (n_reps) on the exact same fitness landscape. To change the number of stochastic realizations, you would simply change the value of the parameter n_reps.   



## Post-Simulation Analyses
Following model simulation completion, several analyses were performed. The code for those analyses can be found in the Post_Simulation_Analysis directory. If an analysis in our paper is not included here, that analysis was performed at the time of model simulation, and the data can be found in the model output data. A


