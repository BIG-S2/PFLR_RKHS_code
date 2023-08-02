# CODE and DATA for “A partially functional linear regression framework for integrating genetic imaging and clinical data”
This is the ReadMe document for running the simulation and real data analysis presented in the paper.

## 1. Real Data Demo
Democode for the ADNI real data analysis in the paper.

### Data
Data used in the paper.
- final_clinical_m12_used.csv: the demographic information
  - the first column is the ID for each subject.
  - the second and the third columns denote Gender and Handedness.
  - the fourth to the sixth columns denote Education length, Retirement status and Age respectively.
  - The seventh to eighth columns denote the disease status (if MCI, the corresponding element would be 1; if AD, the corresponding element would be 1; if the two columns are all zero, it corresponds to CN). 
  - The "APOE4" represents the number of AOPE gene. 
  - PC1  to PC5 contain the leading 5 principal components of all the genetic data.
  - MMSE is the score assessing the cognitive function, with lower values indicating impairment.

- Hippocampus_left.txt and Hippocampus_right.txt: two matrix from left and right hippocampus
  - each is a n* 15000 matrix, where n is the sample size. 
  - For each row, you need to reshape it as a 100*150 matrix (simply use matrix(Hippocampus_left[i:, ], 100, 150) in R). 

- SNP_select_combine.csv
  - genetic data after sure independence screening while controlling the demographic variables and the top 5 PCs.

- SNP_select_label_combine.csv
  - Label informaton of the genetic data. The "chr" column represents the number of chromosome, the "name" is the name of the snp, the "location" represents the position of the SNP. 

### R functions
R functions for real data analysis in the paper.
- PFLR_RKHS_image.R: the proposed estimation method for imaging variable.
- RealData_demo.R: real data analysis demo codes.

## 2. Simulation Demo
Demo code for simulation studies.

### R functions
R functions for simulation.
- Gen_data.R: generate the simulation data.
- PFLR_RKHS_estimate.R: the proposed estimation method for functional variable.
- Simu_demo.R: simulation demo codes.







