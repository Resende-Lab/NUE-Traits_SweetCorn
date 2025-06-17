# NUEtraits_SweetCorn
Repository with the data and scripts used in the paper: "Genomic strategies for nitrogen use efficiency traits breeding in a sweet corn breeding program"


<br>

--------

## data folder

- Folder containing:
    - BLUEs for lines and hybrids
    - Markers for lines and hybrids


## File GenomicModels.R
- File containing the genomic prediction model for NUE traits.
    - BayesB, GBLUP (ST and MT) and Spike-and-Slab models implemented in BGLR.
    - CV0 for the inbred line population for R1.LN, R3.LN, and R6.LN and CV1 (Hybrid prediction) for R3.LN trait. 


## Folder Simulations/

- Multi trait simulation pipeline. We simulated the Sweet corn breeding from University of Florida.
- Three scenarios were implemented:
    - Conv Scenario: phenotypic truncted selection.
    - GS scenario: genomic selection implemented targetting the prediction of individuals.
    - OCS scenario: genomic selection implemented targetting the prediction of individuals and cross prediction and optimization via **SimpleMating** package.

## File CrossPred.R
- Cross prediction for all possible crosses in sweet corn population.
- Prediction and optimation was implemented via **SimpleMating** package.
