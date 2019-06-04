# Cpp code for Simulation

Bayesian Latent Factor on Image Regression with Nonignorable Missing Data

## Getting Started

- Prerequisites:
    - Eigen 3.3.7
    - cmake 3.3.1
    - C++ 14
- Set eigen path in `CmakeLists.txt` by revising `include_directories(<your-path>/eigen3 .)`
- Compile and run as follows.
``` shell
mkdir build
cd build
cmake ..
make
./exe
```
- Outputs
    - Bayesian estimation of replications are saved in `out/rep_*.txt`.
    - Posterior samples of the latest replication are saved in `out/samples_*.txt`

## Documentation

- For changing simulation setting please refer to `Class Constant` including:
    - `setting` (model setup) 
        - 1: Simulation 1 (with image)
        - 0: Simulation 2 
    - `ssl` (Spike-and-Slab prior choice)
        - 1: Adjust SSL
        - 0: SSL
    - `model_type` (missing mechanism)
        - 1: Missing completely at random
        - 2: Missing at random
        - 3: Missing not at random
    - Sample size ($N$), number of eigenimages ($P$), number of repliations, true value, tunning parameters, hyperparameters etc.
- For data generation please refer to `Class Generate`.
- For MCMC initial values, please refer to `Class Initial`.
- All the C++ functions related to posterior sampling are contained in `Class MCMC`.

## Details
- For Simulation 1
    1. Preparation (code is not provided)
        1. Matlab: read preprocessed ADNI MRIs, save as .mat.
        2. Matlab: apply `svd` function for singular value decomposition and acquire eigenscores and eigenimages, save as .txt.
        3. Python: create true image parameter, specify other true values, use the image (.mat) in step a to generate the dateset for simulatin 1, save as .txt.
    2. Cpp
        1. Copy and past the simulated dataset (.txt) into Cpp folder.
        2. Set `setting` to 1 in `Constants.h`. 
        3. Specify simulation setting, such as sample size, number of MCMC iterfation, and true value in `Constants.h` and `Constants.cpp`.
        4. Compile and run.
- For Simulation 2
    - Set `setting` to 0 in `Constants.h`.
    - Specify simulation setting, such as sample size, number of MCMC iterfation, and true value in `Constants.h` and `Constants.cpp`.
    - Compile and run.
      
    

PS: Please feel free to contact the authors if you have any questions!






