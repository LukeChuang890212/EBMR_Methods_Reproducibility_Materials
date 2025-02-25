# EBMR_Methods_Reproducibility_Materials
This repository contains the reproducibility materials for "Ensemble-Based Multiply Robust Methods for Missing Not at Random Data with Model Misspecification". There are three folders: EBMRalgorithm, Simulation_Studies and Mental_Health_Data_Application. Below introduces the details of these three folders.

**EBMRalgorithm**

This folder is a r package that can be loaded into R by the syntax: 
```r
devtools::install_github("LukeChuang890212/EBMR_Methods_Reproducibility_Materials/EBMRalgorithm")".
```
Inside EBMRalgorithm/R, there are five essential R files: EBMRAlgorithm.r, WangShaoKim2014.r, Methods.r, Preprocessor.r and Fool_proofing.r. In EBMRAlgorithm.r, we defined a R6 class which encapsulate all the necessary functions to implement the proposed method. WangShaoKim2014.r contains the function for implementing the proposed method by Wang, Shao and Kim (2014) to estimate the candidate propensity score models, which corresponds to the first step of the proposed ensemble framework. Methods.r contains the functions for implementing the second step of the proposed ensemble framework, including the one that estimate $\bm{\nu}$ and that for computing the proposed inverse probability weighting (IPW) estimator and its asymptotic variance. Finally, Preprocessor.r and Fool_proofing.r contain the technical functions for preprocessing the input arguments and data to make sure the input arguments and data structure are appropriate for the implementation.

# ## Repository Structure
# 
# ```
# EBMR_Methods_Reproducibility_Materials/
# ├── data/                 # Example datasets used in the analysis
# ├── scripts/              # R scripts for model estimation and simulations
# ├── results/              # Generated results from simulations
# ├── figures/              # Plots and visualizations
# ├── docs/                 # Additional documentation
# └── README.md             # Project overview and instructions
# ```

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

