# vRSA_MEEG
Variational RSA for M/EEG

## Setup

1. Download the code and put it somewhere on the MATLAB path.
2. In MATLAB, go into the eeg_example directory and open vRSA_demo_run_all.m . This will download the data and run the analyses.

## Overview
This toolbox implements variational Representational Similarity Analysis (vRSA), a Bayesian extension of traditional RSA that models trial-by-trial covariance in neural
data as a linear mixture of theoretical or empirically derived model covariance components. It is especially suited for time-resolved analyses, where neural responses are 
expressed in terms of temporal basis functions.

## Key features
**Approximate Bayesian Approach**: in vRSA, the trial by trial covariance matrix of the data is modelled as a linear mixture of model covariance 
components (equivalent to model matrices in traditional RSA) using a generalized linear model (GLM). Specifically, each model component is
scaled by a log-space parameter (lamba), ensuring positive contributions. Given its bayesian nature, this approach requires specifying priors
expectation (i.e. mean) and variance for each beta parameters (see below how these are selected). Using variational Laplace, the posterior 
distribution of each parameter is computed, as well as the (approximate) evidence in the form of free energy for each of the parameters. 

**Modelling of temporal dynamics underlying neural representation**
With time-resolved signal, the question is not only whether experimental conditions are reflected in the second order statistics of the data, but also with which temporal
dynamics. With variational RSA, the temporal dynamics of the neural responses can be modelled using within trial basis function, encoding pre-defined temporal dynamics of neural activation
In the example described in this toolbox, we use a **Finite Impulse Response** (FIR) basis set, splitting the analysis in time window to investigate the contribution of experimental contrast over time (see [here](https://github.com/pzeidman/vRSA_MEEG/blob/main/eeg_example/vRSA_demo_run_vRSA_FIR.m)). 
The FIR approach is especially useful when there is no strong prior hypothesis regarding the time course of neural representation and the goal is to investigate the temporal dynamics of the 
experimental contrasts.
However, it is also possible to investigate effects in pre-spefied time window by specifying different basis function. For example, if the aim of the analysis is to test whether a given experimental
contrast is represented in the P100 and P300, one can simply use 2 basis function each with ones in the time windows of each component and 0s elsewhere. 
The within trial basis set can be used to test any predictions regarding the temporal dynamics underlying representations, not only binary (i.e. "windowed") approaches, by using for example neural activation 
components derived from models of neural activation (normative baseline model, neural mass models...). In addition, it is also possible to investigate representations in the frequency domain by modelling responses using a Fourier basis set. This approach allows characterization of representational dynamics in terms of oscillatory components rather than time windows. We provide a tutorial illustrating this frequency-domain analysis as well (see [here](https://github.com/pzeidman/vRSA_MEEG/blob/main/eeg_example/vRSA_demo_run_vRSA_FourierSet.m)).

**Hierarchical (Group) Analysis**. Single-subject vRSA models are estimated to obtain posterior distributions over lambda parameters. These 
individual posteriors then feed into a parametric empirical Bayes (PEB) step, which combines them for group-level inference. Since PEB 
accounts for uncertainty in each subject’s estimates, it provides more robust results compared to naive correlation-based group analyses in
traditional RSA.

## Advantages Over Traditional RSA

- vRSA models all covariance components simultaneously within a single generative model, instead of correlating data with each model 
separately (i.e. correlating the observed distance matrix separately with each model matrix). 
- The Bayesian framework yields for each parameter a posterior distribution, capturing the mean but also variance of the estimates. This enables
accurate weighting of each subject's data at the group level, to account for differences in the uncertainty across participants. 
- The time-resolved framework, meanwhile, avoids the assumption that similarity emerges uniformly over a single time window, making 
vRSA especially useful when the temporal profile of neural responses is unknown or complex.

