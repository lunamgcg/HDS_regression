# Estimation of conversion rates and selectivity in the hydrodesulfurization of dibenzothiophene (DBT) using supervised machine learning.

The hydrodesulfurization (HDS) process is commonly used to remove sulfur compounds from fuels, but eliminating dibenzothiophene (DBT) and its derivatives is challenging. 
This study examines factors affecting catalyst efficiency in HDS of DBT using machine learning (ML) algorithms.
It is applied Lasso, Ridge, and Random Forest regression techniques to estimate DBT conversion and selectivity. 
Random Forest and Lasso provided reliable conversion predictions, while regularized regression methods were effective for selectivity prediction.
Key structural parameters, such as pore size and slab length, significantly influence selectivity by affecting active site availability. 

# Objetive
There exist several kinds of catalysts for HDS. However, removing DFT is a challenge for these materials. This contribution considers data on MoS2-derived catalysts from scientific articles in Scopus, Web of Science, and SciFinder. 
The main objective is to estimate the objective variables, namely the conversion rate and selectivity of catalysts in the HDS of DBT. Moreover, the features are identified as influencing the objective variables. 
For the estimation of objective variables, the following methods are used: 
- Lasso regression 
- Ridge regression
- Random forest regression
The identification of features relevant to objective variables is performed by analyzing the regression coefficient.

# Dataset

The database compiles data from scientific articles in Scopus, Web of Science, and SciFinder, focusing on publications from 2005 to 2022. This period was selected due to the development of the Nebula catalyst in the early 2000s, which initiated new research in deep hydrodesulfurization (HDS). Articles were included if they provided essential information about the HDS process, encompassing catalyst composition, structural parameters, and reaction conditions, with 30 descriptors considered.

Out of 81 papers discussing the conversion of dibenzothiophene (DBT), selectivity, and catalyst material derived from molybdenum disulfide (MoS2), only 19.7% were used to build the database due to the absence of crucial structural parameters in many articles. This subset identified 86 different catalysts. Outliers with high values were adjusted to the maximum value of the feature to maintain data integrity.

# Methodology

-Preprocessing
To clean the data, any records with null values were removed. Out of 30 predictors, 21 are related to composition, 3 pertain to experimental conditions, and 5 relate to structural parameters. Outliers with high values were replaced with the maximum value of the respective feature. The composition predictors were transformed into nominal values. Finally, all data were standardized so that each variable has a minimum value of 0 and a maximum value of 1.


<img width="442" height="647" alt="image" src="https://github.com/user-attachments/assets/294f656d-66ea-497e-851f-4e3c86774eaf" />

