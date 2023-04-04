CODE for paper “A partially functional linear regression framework for integrating genetic imaging and clinical data”

****************

Demo code for simulation is in the R file "Simu_demo" of the "Simulation demo folder"

****************
Democode for the real data analysis is in R file "RealData_demo.R" in the folder “Rea Data demo”


********************************************************************************
The "Data" folder in the folder “Rea Data demo” contains the real data:
********************************************************************************

final_clinical_m12_used.csv: it contains the demographic information, the first column is the ID for each subject, the second and the third columns denote Gender and Handedness, the fourth to the sixth columns denote Education length, Retirement status and Age respectively. The seventh to eighth columns denote the disease status (if MCI, the corresponding element would be 1; if AD, the corresponding element would be 1; if the two columns are all zero, it corresponds to CN). The "APOE4" represents the number of AOPE gene.
PC-1  to PC-5 contain the leading 5 principal components of all the genetic data.
MMSE is the score assessing the cognitive function, with lower values indicating impairment.


Hippocampus_left.txt and Hippocampus_right.txt contains two matrix from left and right hippocampus, each is a n*15000 matrix, where n is the sample size. For each row, you need to reshape it as a 100*150 matrix (simply use matrix(yalingeneticsright19[i,:],100,150) in R) 

SNP_select_combine.csv: it contains genetic data after sure independence screening while controlling the demographic variables and the top 5 PCs.
SNP_select_label_combine.csv: it contains the label informaton of the genetic data. The "chr" column represents the number of chromosome, the "name" is the name of the snp, the "location" represents the position of the SNP. 




