# GIR-Combat
GIR-Combat: Batch effect correction method for scRNA-seq data with large variance based on global information reference batches

## Overview
This repository contains the code for the paper "GIR-Combat: Batch effect correction method for scRNA-seq data with large variance based on global information reference batches".

## Regenerating the figures
### Usage
1. The implementation of the code is based on R language.
2. All the R packages needed for the experiment are included in the code.
3. Our project consists of five experiments, three of which are simulation experiments and two are real data experiments. Each experiment is comprised of four .R files, named part1 to part4. For each experiment, in order to reproduce the results, you need to perform Parts 1 through 4 sequentially.
4. A Gaussian smoothing kernel function is implemented in dusmooth_gaussian_kernel.cpp, which is suitable for data analysis in R.
5. For each experiment:
   - part1 is used for SPKM clustering and forming reference batch.
   - part2 and part3 are used to  merge Correction Matrix.
   - part4 is used for empirical bayesian estimation and removing reference batch. And in part4, the figures of the experiments are generated.

### Simulations
To run the simulations, enter the sim1/sim2/sim3 directory and run part1 to part4 sequentially.
1. In sim1, a balanced dataset is utilized, characterized by approximately equal variances across all sample label types within the simulated data.
2. In sim2, an equilibrium dataset featuring large variance is utilized, in which the variance of the normal distribution fitted to one specific sample label in the simulated data was markedly greater than that of the other sample label. 
3. In sim3, a large variance unbalanced dataset is utilized, whereby the variance of a normal distribution fitted to one of the sample labels in the simulated data not only displayed large variance but also revealed that a batch of the simulated data lacked the presence of all sample labels.
### Real data

## Maintainers
@BioStatAnalysis1

## Contributing
Feel free to dive in! Open an issue or submit PRs.
