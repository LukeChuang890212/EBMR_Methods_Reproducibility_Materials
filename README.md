# EBMR_Methods_Reproducibility_Materials

This repository contains reproducibility materials for “Ensemble-Based Multiply Robust Methods for Missing Not at Random Data with Model Misspecification.” It includes three main folders:
-	EBMRalgorithm – R package implementing the proposed method.
-	Simulation_Studies – Scripts for running simulation experiments.
-	Mental_Health_Data_Application – Application of the method to real-world mental health data.
	
**EBMRalgorithm**

The EBMRalgorithm folder contains an R package that can be installed using:
```r
devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")".
```
Inside EBMRalgorithm/R, there are five essential R files: EBMRAlgorithm.r, WangShaoKim2014.r, Methods.r, Preprocessor.r and Fool_proofing.r. In EBMRAlgorithm.r, we defined a R6 class which encapsulate all the necessary functions to implement the proposed method. WangShaoKim2014.r contains the function for implementing the proposed method by Wang, Shao and Kim (2014) to estimate the candidate propensity score models, which corresponds to the first step of the proposed ensemble framework. Methods.r contains the functions for implementing the second step of the proposed ensemble framework, including the one that estimate $\bf{\nu}$ and that for computing the proposed inverse probability weighting (IPW) estimator and its asymptotic variance. Finally, Preprocessor.r and Fool_proofing.r contain the technical functions for preprocessing the input arguments and data to make sure the input arguments and data structure are appropriate for the implementation.

Key Files in EBMRalgorithm/R
-	EBMRAlgorithm.r – Defines an R6 class encapsulating all essential functions for the proposed method.
-	WangShaoKim2014.r – Implements the method by Wang, Shao, and Kim (2014) for estimating candidate propensity score models (first step of the ensemble framework).
-	Methods.r – Implements the second step of the ensemble framework, including estimation of ν and computation of the inverse probability weighting (IPW) estimator with its asymptotic variance.
-	Preprocessor.r & Fool_proofing.r – Handle input validation and data preprocessing to ensure correct argument and data structures for implementation.
	
## Prerequisites

Ensure you have the following software and dependencies installed:

- **R (>=4.0)**
- Required R packages:
  ```r
  install.packages(c("ggplot2", "dplyr", "tidyr", "MASS", "Matrix"))
  ```

You may also need additional dependencies as specified in individual scripts.

## Installation & Setup

Clone this repository to your local machine:
```sh
git clone https://github.com/LukeChuang890212/EBMR_Methods_Reproducibility_Materials.git
cd EBMR_Methods_Reproducibility_Materials
```

## Usage

Run the main analysis script:
```sh
Rscript scripts/main_analysis.R
```

For simulations:
```sh
Rscript scripts/simulation_study.R
```

Generated results will be saved in the `results/` directory.

## Reproducibility

To ensure full reproducibility, check the session information:
```r
sessionInfo()
```

If you encounter any package version issues, consider using [renv](https://rstudio.github.io/renv/) for dependency management.
```r
install.packages("renv")
renv::init()
```

## Contributing

If you'd like to contribute, please fork the repository and submit a pull request.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

For any questions, feel free to reach out to [Luke Chuang](https://github.com/LukeChuang890212).

