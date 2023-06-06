# R-SI-balance-score Pipeline
This pipeline calculates the R/SI balance score based on two gene expression maps (Figure A): Resistance(R)/Tolerance(T) (Cohn, O. et al., 2022) and MetS/Systemic Inflammation(SI) (Frishberg, A. et al., 2021).
## Pipeline Overview
The pipeline, illustrated in Figure B, comprises several methods implemented in ***'states_pipeline.py'***. It enables the calculation of the full R/SI balance score by integrating the R/T and MetS/SI gene expression maps. To obtain the complete pipeline results from the input matrix, use the ***'full_R_SI_score_pipeline(expression_df, control_individuals)'*** method. 
The input is a gene expression matrix with individuals as columns and genes as the index, and a list of columns representing control individuals.

Note: It is also possible to calculate the state for any other gene map, as long as the index represents genes and the columns represent x and y positions. This allows you to analyze additional gene maps of interest. Please refer to the running example file to see how to perform this calculation.
## File Descriptions
- ***states_pipeline.py***: Contains the implementation of the pipeline methods
- ***running_example.ipynb***: Jupyter Notebook file demonstrating the usage of the pipeline with a running example.
- ***resistance_tolerance_map.csv***: Gene expression map for resistance/tolerance.
- ***MetS_systemic_inflammation_map.csv***: Gene expression map for MetS/systemic inflammation
## Usage
Make sure ***'resistance_tolerance_map.csv'*** and ***'MetS_systemic_inflammation_map.csv'*** files in the project folder

Example code to run the pipeline in python:
```python 
from states_pipeline import full_R_SI_score_pipeline

# Load the gene expression matrix and control individuals
expression_df = pd.read_csv('expression_matrix.csv',index = 'genes')
control_individuals = ['control1', 'control2', 'control3']#subset of expression_df columns

# Run the full pipeline
results = full_R_SI_score_pipeline(expression_df, control_individuals)
```

## References
1.	Cohn, O. et al. Distinct gene programs underpinning disease tolerance and resistance in influenza virus infection. Cell Syst 13, 1002-1015.e9 (2022).
2.	Frishberg, A. et al. An integrative model of cardiometabolic traits identifies two types of metabolic syndrome. eLife 10, e61710 (2021).
 



![pipeline](https://github.com/rachelbl2/R-SI-balance-score/assets/81696220/4461b7d6-71ab-4614-b4d4-912ba1c7763b)
