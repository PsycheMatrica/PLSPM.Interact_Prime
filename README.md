# PLSPM.Interact_Prime

### Author:
Heungsun Hwang and Gyeongcheol Cho

## Description:
- The **PLSPM.Interact_Prime** package enables users to estimate and evaluate basic PLSPM models.

## Features:
- Estimate the parameters of PLSPM with construct interactions and calculate their standard errors (SE) along with 95% confidence intervals (CI).
- Enable parallel computing for bootstrap sampling.
- Allow users to specify a sign-fixing indicator for each component.

## Installation:
To use this package in MATLAB:
1. Clone or download the repository:
   ```bash
   git clone https://github.com/PsycheMatrica/PLSPM.Interact_Prime.git
   ```
2. Add the package to your MATLAB path:
   ```matlab
    addpath(genpath('PLSPM.Interact_Prime'))
   ```

## Usage:
- For examples on how to use the package, refer to the `Run_Example_InteractPLSPM.m` file. This file demonstrates the implementation of `InteractPLSPM()` using the ACSI dataset.

## Compatibility:
- Tested on MATLAB R2023b.
- Likely compatible with earlier MATLAB versions.

### Citation (APA):
- If you use **PLSPM.Interact_Prime** in your research or publications, please cite it in APA format as follows:

```plaintext
Hwang, H. & Cho, G. (2024). PLSPM.Interact_Prime: A package for partial least squares path modeling with construct interactions [Computer software]. GitHub. https://github.com/PsycheMatrica/PLSPM.Interact_Prime
```
