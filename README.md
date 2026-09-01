# ARLSA: Analytic and Robust Semi-Supervised Adaptation under Label Shift

This repository contains research code for a simulation study based on **Analytic and Robust Semi-Supervised Adaptation under Label Shift (ARLSA)**.

## Overview

The project considers a label-shift setting where labeled observations are available from a source population, while only unlabeled observations are available from a target population with a different label distribution.

The main goal of this repository is to implement the **ARLSA** estimation procedure and investigate its finite-sample performance through simulation experiments.

The implementation also includes benchmark estimators to compare the performance of **ARLSA** under different simulation settings.

## Methods

The simulation study mainly considers the following estimation methods:

**·ARLSA estimator**

**·Importance-weighted estimator**

**·Naive estimator**

These methods are compared to evaluate the effect of accounting for label shift in **target-population estimation**.

## Simulation Study

Simulation experiments are conducted by repeatedly generating labeled source data and unlabeled target data under controlled label-shift settings.

The experiments are used to evaluate the finite-sample behavior of the estimators.

Performance is summarized using measures such as:

·**Bias**

·**Standard Error**

·**Mean Squared Error (MSE)**

·**Root Mean Squared Error (RMSE)**

Additional graphical summaries are used to examine the simulated data distributions and estimator performance.

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

The project is implemented in R.

Required packages are specified in the source files.

## Research Status

This repository is provided for research and educational purposes to document the implementation and simulation study of the **ARLSA** methodology.

The code represents an independent implementation based on the methodology described in the original **ARLSA** paper and is not an official implementation provided by the authors.

# Acknowledgement

This repository was developed with support from the 서울시립대학교 데이터 사이언스 플러스 차세대 융합인재 양성사업단 
- http://dsplus.uos.ac.kr/
