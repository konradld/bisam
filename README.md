# BISAM: Estimating Bayesian Indicator Saturated Models

This repository contains R functions for Bayesian structural break detection and analysis, with comparison tools for evaluating performance against Gets (General-to-Specific) estimation.

## Overview

BISAM provides a Bayesian approach to detecting structural breaks in time series data using Bayesian model averaging and selection. The methodology uses Posterior Inclusion Probabilities (PIPs) to identify step-breaks and outliers and select break points based on configurable thresholds.

## Repository Structure

### Core Functions

- **`estimate_bisam_fun.R`**  
  Main estimation function for Bayesian structural break analysis. This is the key function of the repository.

- **`pip_window_fun.R`**  
  Calculates Joint Posterior Inclusion Probabilities (PIPs) and selects break points according to a specified threshold.

- **`contr_sim_breaks_fun.R`**  
  Data simulation function that generates panel data with structural breaks and outliers of specified sizes and location.

### Installation

- **`install_bespoke_mombf_3.5.4.R`**  
  Installs a bespoke version of the `mombf` package (based on version 3.5.4) required for the estimation functions. This modified version is included in the repository and must be installed before using `estimate_bisam_fun.R`.

### Examples and Simulations

- **`00_bisam_simulation.R`**  
  Example script demonstrating the basic workflow: data generation using `contr_sim_breaks_fun.R` and estimation using `estimate_bisam_fun.R`.

- **`01_estimation_gets_bisam_comparison.R`**  
  Simulation study comparing BISAM performance against Gets (General-to-Specific) estimation.

- **`02_analysis_gets_bisam_comparison.R`**  
  Analysis and visualization of results from the output of `01_estimation_gets_bisam_comparison.R`.

## Installation

1. Clone this repository:
```r
# In your terminal
git clone https://github.com/konradld/bisam.git
cd bisam
```

2. Install the bespoke mombf package:
```r
source("install_bespoke_mombf_3.5.4.R")
```

3. Load required functions:
```r
source("estimate_bisam_fun.R")
source("pip_window_fun.R")
source("contr_sim_breaks_fun.R")
```

## Quick Start

See the example simulation for a complete workflow in **`00_bisam_simulation.R`**

## References

tba

## Contact

<lucas.konrad@wu.ac.at>

## Citation

If you use this code in your research, please cite:

```
tba
```
