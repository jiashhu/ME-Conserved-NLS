This is a reproducibility repository for the paper "High-order mass- and energy-conserving methods for the nonlinear Schrödinger equation" submitted to the SIAM Journal for Scientific Computing (SISC) by Genming Bai, Jiashun Hu, and Buyang Li in 2023. The repository serves two purposes:

* Establishing and maintaining the sources of all data, and figures presented in the paper. In particular, the repository archives all raw data and their sources for the visualizations presented in the paper. In most cases, simulations involve many time steps before reaching the final solution, and computational measurements only monitor a small portion of the simulation domain, it is sufficient to archive the actual raw data used to generate the figures themselves.

* Assisting and automating reproducibility. This repository describes the hardware/software environment and software used for each experimental result reported in the paper. 

The automation of this repository relies on Python (including ngsolve/netgen, numpy, matplotlib, and other packages). 

For the installation of ngsolve/netgen, we recommend readers to follow the tutorial at https://ngsolve.org/downloads.

Finally, the REPRODUCED file lists the readers who have been able to reproduce the experiments in this repository and notes any discrepancies found.

If you reproduce any data from this repository, please inform me or the other authors of this paper so that we can add your results to the REPRODUCED file.

## Structure of the repository

* Implementation of a fully discrete format for high-order mass and energy-conserving methods for solving the nonlinear Schrödinger equation. The core classes can be found in ./Collocation/BaseMethod.py, which utilize implicit Newton iteration methods for solving the nonlinear system.

* To run the code, add the path of the basic codes and AuxiliaryPackages to PYTHONPATH.

```
export PYTHONPATH=$PYTHONPATH:../AuxiliaryPackages
```