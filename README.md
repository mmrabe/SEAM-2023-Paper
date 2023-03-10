## SEAM (Rabe et al., 2023) Supplementary Material

### General

This repository contains the analysis code as used in the Rabe et al. (2023) SEAM manuscript.

The RMarkdown script SEAM.Rmd can be run with cached intermediate results without having to compile the SEAM and SWIFT models. However, for a full replication, you should empty the cache directory and compile the models yourself so that you can perform the necessary simulations.

### Compiling SWIFT and SEAM

This repository contains the computational models SEAM and SWIFT as used in the Rabe et al. (2023) SEAM manuscript.

Each model resides in its own subdirectory (`SEAM/` or `SWIFT/`, respectively). The models are coded in C but some supplementary files to prepare compilation are coded in Python. In order for the models to run properly, they must first be compiled to binary executables, which is different from machine to machine.

On UNIX systems (such as Linux or MacOS), compilers are usually pre-installed but we recommend the latest GNU compiler (GCC) to compile the files. You can use the `./build.sh` build script inside each directory for compiling all executables for the respective model, e.g.

```
> cd SEAM
> ./build.sh
> cd ../SWIFT
> ./build.sh
```

Upon compiling, a `*/SIM/` subdirectory will be created for the respective model, which contains all binary executables. The libraries for use with the scripting languages R and Python are moreover copied into the `*/MCMC/` subdirectories. There, you will also find example scripts in R and Python, and a shell script calling the command-line interface. In the examples, the corpus and an example fixation sequence are loaded, parameters are manually updated, and the likelihood is evaluated.

Each model comes with an example corpus and dataset in `*/DATA/`. The corpus represents the items used in the experimental Mertzen et al. (2023) data set. The test data set is a simulated fixation sequence.

### Running the MCMC simulations

The script swift_recovery.R performs the recovery study (fit SEAM to 100 simulated datasets with known parameters. The script swift_fit.R performs the parameter estimation for experimental data.
