# BNP-WMAR Example Code

This folder contains code to fit the models and perform simulation studies as described in the associated article by Heiner and Kottas. All MCMC is run in [Julia](https://julialang.org/) v1.1+, as implemented in the [BNP_WMReg_Joint](https://) package. One can work directly with the .jl scripts, or call them with the .sh scripts.

## Packages

Before running the .jl scripts, it is necessary to install the packages that appear after all instances of `using` in the .jl scripts. To do this, open the Julia REPL and enter package mode by pressing `]`. In package mode, enter `add` followed by each package name, separated by spaces. The primary package is unregistered and can be installed with

```julia
pkg> add https://github.com/mheiner/BNP_WMReg_Joint.git
```

Post processing is run with R scripts, which require the coda R package.

## Data

The data folder contains the pink salmon data set (in the public domain, see <https://inport.nmfs.noaa.gov/inport/item/17256>), as well as code to generate the simulation data (by running the shell script `simulate_Ricker.sh`).

## Models

The fastest way to fit models illustrated in the article is to source the shell scripts `runjobs_Ricker.sh` and `runjobs_pinkSalmon.sh`, which call `runMCMC.jl`. Alternatively, `runMCMC.jl` can be run interactively in Julia (in which case, the user would manually assign the global variables normally passed by `ARGS`). Results in the article are based on longer MCMC runs than the current settings in `runMCMC.jl`.

Following MCMC, run 
```bash
source runpostproc.sh 112119
```
(possibly replacing 112119 with an integer identifying specific runs), which will call `postsim_to_R.jl` and `postProcess_loop.R` to generate basic posterior summaries.

The script `postProcess.jl` contains several useful code snippets for exploring posterior inferences interactively. The script `KL.jl` supports the simulation study from the article, and should be run by sourcing `runKL_Ricker.sh` after `runpostproc.sh`. WARNING: `KL.jl` employs Monte Carlo approximation using 2,000 replicates for each of 1,000 validation time points and each posterior sample; it is prohibitively slow for long MCMC chains.

## Folders

- data: Contains data.
- KL: Density estimation performance metrics from the Ricker simulation study are collected in this folder.
- logs: Shell output from `runjobs` scripts is redirected to this folder.
- plots: The script `postProcess_loop.R` saves plots to the trace/ sub-directory.
- postsim: Posterior simulations are saved to this folder.
- postsimB: The script `postsim_to_R.jl` moves poster simulation .bson files to this folder.
- postsim_progress: Files reporting MCMC progress and statistics are saved in this folder.
