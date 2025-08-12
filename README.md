# LargeCoheEsti

Supporting Code for "Reliable Estimation of Coherence from Scarce Data: Theory and Experiment"

## Introduction

This repository contains the official source code for the paper "Reliable Estimation of Coherence from Scarce Data: Theory and Experiment".

It includes all the necessary code to reproduce the figures, generate the experimental data, and validate the algorithms presented in the paper.

## Code Structure

This repository is organized into two main parts: MATLAB and Python code. These are located in separate directories and are used to support the various experiments and analyses in our work.

## Requirements

Please ensure that your environment has the following prerequisites installed.

### MATLAB
To successfully run the MATLAB scripts, you will need the following toolboxes:

CVX: https://cvxr.com/cvx/ A modeling system for constructing and solving disciplined convex programs.

CVXQUAD: https://github.com/hfawzi/cvxquad An extension for CVX that adds support for functions involving the matrix logarithm, exponential, and entropy.

Solver: We recommend using MOSEK as the solver for CVX to achieve the best performance and stability. However, other CVX-compatible solvers (e.g., SDPT3, SeDuMi) are also supported.

### Python
The Python environment requires the following standard scientific computing libraries:

pandas numpy matplotlib torch scipy

You can quickly install these dependencies using pip:

pip install pandas numpy matplotlib torch scipy

## Usage
Data Generation & Algorithm Validation: These tasks are primarily handled by the MATLAB scripts. Please refer to the comments within the respective files for instructions on how to run them.

Figure Plotting: The figures are generated using Python scripts. These scripts read the data files (in .mat or .csv format) produced by MATLAB and create the plots shown in the paper.

## Citation
If you use this code or our methods in your research, please cite our paper:
