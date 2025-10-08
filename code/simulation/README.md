# FSI-Simulation

This folder contains the scripts for the fluid–structure interaction (FSI) simulation presented in the paper.
In contrast to the other numerical experiments, this case was performed using an updated version
of TrixiParticles.jl (v0.4.1) and a more recent Julia release (v1.11.7).
To ensure reproducibility, a dedicated Julia project environment is provided in this directory.

## Set up Julia
Install Julia following the instructions at https://julialang.org/downloads/. We have used Julia
v1.11.7.

To install all necessary Julia packages, execute the following statement from within the folder that
contains the `README.md` file you are currently reading:
```shell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```
This will recreate the exact Julia environment we used to obtain our results for full
reproducibility.

## Run the experiments
To re-create the results, run one of the following commands.
Make sure the [`../data`](../data/) folder with all required input data is present.

### Create all initial conditions
```bash
julia --project=. --threads=32 run_setup_initial_condition.jl
```

### Run the simulations
```bash
julia --project=. --threads=32 run_simulation.jl
```

### Run all steps together
```bash
julia --project=. --threads=32 create_data.jl
```
