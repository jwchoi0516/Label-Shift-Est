# Simulation Study for Estimation under Label Shift

This repository contains research code for a simulation study based on **Simulation Study for Estimation under Label Shift**.

## Overview

The project considers a **label-shift** setting where labeled observations are available from a **source population**, while only unlabeled observations are available from a **target population** with a different label distribution.

The main goal of this repository is to implement the **Label Shift** procedure and investigate its finite-sample performance through simulation experiments.

The implementation also includes benchmark estimators to compare the performance of **Label Shift** under different simulation settings.

## Methods

The simulation study mainly considers the following estimation methods:

**● Label Shift estimator**

**● Importance-weighted estimator**

**● Naive estimator**

These methods are compared to evaluate the effect of accounting for label shift in **target-population estimation**.

## Simulation Settings

| Parameter | Value / Distribution | Description |
|---|---|---|
| **Replications** | `300` | Number of independent simulation replications |
| **Total Sample Size** | `1000` | Total number of source and target observations |
| **Source Sample Size** | `500` | Number of labeled source observations |
| **Target Sample Size** | `500` | Number of unlabeled target observations |
| **Source Proportion** | `0.5` | Proportion of source observations in the combined sample |
| **Source Y** | `Normal(0, 2)` | Source outcome distribution |
| **Target Y** | `Normal(1, 1)` | Target outcome distribution |
| **X1** | `-0.5Y + ε` | First covariate, with ε ~ Normal(0, 1) |
| **X2** | `0.5Y + ε` | Second covariate, with ε ~ Normal(0, 1) |
| **X3** | `Y + ε` | Third covariate, with ε ~ Normal(0, 1) |
| **Target Parameter** | `θ = (b0, b1, b2, b3)ᵀ` | Target linear regression parameter |
| **Compared Methods** | `Naive, IPW, ARLSA, Oracle` | Estimators compared in the simulation |

## Simulation Results

| Method | beta0 | beta1 | beta2 | beta3 | RMSE | MAE |
|:------|------:|------:|------:|------:|------:|------:|
| **Naive**  | 0.0015 | -0.2510 | 0.2492 | 0.5004 | 0.2105 | 0.1503 |
| **IPW**    | 0.3473 | -0.2143 | 0.2173 | 0.4323 | 0.0528 | 0.0409 |
| **Label Shift**  | 0.4183 | -0.2043 | 0.2063 | 0.4108 | 0.1763 | 0.0659 |
| **Oracle** | 0.3971 | -0.2011 | 0.1985 | 0.4016 | 0.0290 | 0.0233 |

## Repository Structure
01_setup.R
    Basic simulation settings, data-generating mechanisms,
    model parameters, and functions used in the experiments.

02_simulation.R
    Main simulation code for generating source and target data
    and computing the estimators across repeated experiments.

03_results.R
    Code for summarizing simulation results
    and producing tables and graphical outputs.
## Environment

The project is implemented in **R**.

Required packages are specified in the source files.

## Research Status

This repository is provided for research and educational purposes to document the implementation and simulation study of the **ARLSA** methodology.

The code represents an independent implementation based on the methodology described in the original **ARLSA** paper and is not an official implementation provided by the authors.

## Conclusion
<img width="1165" height="783" alt="sim300_boxplot" src="https://github.com/user-attachments/assets/022824a3-92a3-4013-a906-9783e95b64ef" />
<img width="1109" height="783" alt="sim300_rho" src="https://github.com/user-attachments/assets/410ded6a-eb04-4023-8953-dab50695f5a3" />
<img width="1109" height="783" alt="sim300_xdist" src="https://github.com/user-attachments/assets/f5514f7b-285d-408b-8651-2ce3631f7cc3" />



# Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 
- http://dsplus.uos.ac.kr/
