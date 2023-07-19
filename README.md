#Personalized inflammatory scores Pipeline
This pipeline provides four personalized inflammatory scores: antimicrobial resistance (R), Disease Tolerance (T), Metabolic syndrome (MetS) and Systemic Inflammation (SI). The calculation of personalized SI/MetS states are based on the MetS/SI map (Frishberg, A. et al., 2021), whereas the calculation of the R/T states are based on the R/T map (Cohn, O. et al., 2022). The output also includes the personal balance between R and SI (R/SI balance score). 
##Pipeline Overview
The pipeline, illustrated in Figure B, comprises several methods implemented in 'states_pipeline.py'. It enables the calculation of the four inflammatory states (R, T, SI, MetS and SI) as well as the R/SI balance score, by integrating the input R/T and MetS/SI gene expression maps. To obtain the complete pipeline results from the input matrix, use the 'full_inflammation_scores_pipeline (expression_df, control_individuals)' method. The input is a gene expression matrix with individuals as columns and genes as the index, and a list of columns representing control individuals.
The control individuals are used for standardization of the states and facilitating comparisons between the states derived from the MetS/SI map and those derived from the T/R map, thus enabling the calculation of the R/SI state. Additionally, this allows the comparison between different cohorts. In case your data does not include control individuals, provide a list of all individuals in the cohort instead.Note: It is also possible to calculate the state for any other gene map, as long as the index represents genes and the columns represent x and y positions. This allows you to analyze additional gene maps of interest. Please refer to the running example file to see how to perform this calculation.

##File Descriptions
•	states_pipeline.py: Contains the implementation of the pipeline methods
•	running_example.ipynb: Jupyter Notebook file demonstrating the usage of the pipeline with a running example.
•	resistance_tolerance_map.csv: Gene expression map for resistance/tolerance.
•	MetS_systemic_inflammation_map.csv: Gene expression map for MetS/systemic inflammation
##Usage
Make sure 'resistance_tolerance_map.csv' and 'MetS_systemic_inflammation_map.csv' files in the project folder
Example code to run the pipeline in python:
```python 
from states_pipeline import full_R_SI_score_pipeline

# Load the gene expression matrix and control individuals
expression_df = pd.read_csv('expression_matrix.csv',index = 'genes')
control_individuals = ['control1', 'control2', 'control3']#subset of expression_df columns

# Run the full pipeline
results = full_inflammation_scores_pipeline(expression_df, control_individuals)
```
The full_inflammation_scores_pipeline function returns a data frame with the T, R, MetS, SI, and the R/SI balance scores with the corresponding p values (columns) for each individual (index). 

##References
1.	Cohn, O. et al. Distinct gene programs underpinning disease tolerance and resistance in influenza virus infection. Cell Syst 13, 1002-1015.e9 (2022).
2.	Frishberg, A. et al. An integrative model of cardiometabolic traits identifies two types of metabolic syndrome. eLife 10, e61710 (2021).
 
