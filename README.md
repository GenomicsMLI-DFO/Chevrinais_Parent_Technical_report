# Scripts of analyse and figure for pilot study in estuarine environmental DNA

__Main author:__  Marion Chevrinais  
__Affiliation:__  Fisheries and Oceans Canada (DFO)   
__Group:__        Laboratory of genomics   
__Location:__     Maurice Lamontagne Institute, Mont-Joli, Québec, Canada  
__Affiliated publication:__ Chevrinais and Parent. Improving estuarine water environmental DNA detections using tank experiment. Canadian Technical Report of Fisheries and Aquatic Sciences.       
__Contact:__      e-mail: marion.chevrinais@dfo-mpo.gc.ca 

- [Description](#description)
- [Statistics](#statistics)
- [References](#references)

## Description 

Welcome to this repository with the raw data and code necessary to make the analyses of the pilot study in estuarine environmental DNA paper. 

## Statistics

The code **statistics.R** includes all the code to make statistics from raw data using generalised linear mixed models followed by least-square means pairwise comparisons for each experiment. We estimated the effect of the following categorical variables on eDNA copies detected. 
- Experiment 1: sample preservation treatments and duration
- Experiment 2: filter treatments
- Experiment 3: filter preservation treatments and duration and their interactions
- Experiment 4: extraction treatments
Models include a mixed effect nesting qPCR replicates into DNA extracts ro avoid pseudo-replication. 
