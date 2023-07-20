# Personalized inflammatory scores Pipeline
This pipeline provides four personalized inflammatory states: antimicrobial resistance (R), Disease Tolerance (T), Metabolic syndrome (MetS) and Systemic Inflammation (SI). The calculation of personalized MetS/SI states are based on the MetS/SI gene map (Frishberg, A. et al., 2021), whereas the calculation of the R/T states are based on the R/T gene map (Cohn, O. et al., 2022). The output also includes the personal balance between R and SI (the “R/SI-balance score”, Brandes et al. 2023).
## Pipeline Overview
The pipeline, illustrated in Figures A,B, is implemented in states_pipeline.py. The input to the full_inflammation_scores_pipeline(expression_df, control_individuals) method is (i) a gene expression matrix with individuals as columns and genes as the index, (ii) a list of columns representing control individuals in the gene expression matrix.
In addition to these inputs, the method utilizes the following files, which should be placed in the code folder: (i)The R/T gene map, (ii) The MetS/SI gene map.
The method provides the following output: The four inflammatory states (R, T, MetS, and SI) and their p-values, as well as the R/SI-balance score, for each individual (see details for the calculations in Brandes et al., 2023). 

The control individuals are used for standardization of the states and facilitating comparisons between the states derived from the MetS/SI map and those derived from the R/T map, thus enabling the calculation of the R/SI-balance score. Additionally, this allows the comparison between different cohorts. If the dataset does not include control individuals, it is possible to use all individuals in the cohort for the standardization step (i.e., the input list of control individuals will contain all individuals in the cohort).

Note: It is also possible to calculate the state for any other input gene map, as long as the index represents genes and the columns represent the x and y positions of these genes. This enables analysis of any gene map of interest. Please refer to the running example file to see how to perform this calculation.

## File Descriptions
•	***states_pipeline.py***: Contains the implementation of the pipeline methods.

•	***running_example.ipynb***: Jupyter Notebook file demonstrating the usage of the pipeline with a running example.

•	***resistance_tolerance_map.csv***: The map for the R/T.

•	***MetS_systemic_inflammation_map.csv***: The map for MetS/SI.

## Usage
Make sure ***'resistance_tolerance_map.csv'*** and ***'MetS_systemic_inflammation_map.csv'*** files in the project folder.
Example code to run the pipeline in python:
```python 
from states_pipeline import full_R_SI_score_pipeline

# Load the gene expression matrix and control individuals
expression_df = pd.read_csv('expression_matrix.csv',index = 'genes')
control_individuals = ['control1', 'control2', 'control3']#subset of expression_df columns

# Run the full pipeline
results = full_inflammation_scores_pipeline(expression_df, control_individuals)
```
The ***full_inflammation_scores_pipeline*** function returns a data frame with the T, R, MetS, SI states, their p-values, and the R/SI balance scores (columns), for each individual (index).

## References
1.	Cohn, O. et al. Distinct gene programs underpinning disease tolerance and resistance in influenza virus infection. Cell Syst 13, 1002-1015.e9 (2022).
2.	Frishberg, A. et al. An integrative model of cardiometabolic traits identifies two types of metabolic syndrome. eLife 10, e61710 (2021).
3.	Brandes R. et al. Manuscript in revision.
 

![pipeline v2](https://github.com/rachelbl2/Personalized-inflammatory-scores-Pipeline/assets/81696220/44467dd1-b63e-445e-877e-37507c86a53a)
