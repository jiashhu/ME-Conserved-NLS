This is a reproducibility repository for the paper "High-order mass- and energy-conserving methods for the nonlinear Schrödinger equation" submitted to the SIAM Journal for Scientific Computing (SISC) by Genming Bai, Jiashun Hu, and Buyang Li in 2023. The repository serves two purposes:

* Establishing and maintaining the sources of all data, and figures presented in the paper. In particular, the repository archives all raw data and their sources for the visualizations presented in the paper. 

* Assisting and automating reproducibility. This repository describes the hardware/software environment and software used for each experimental result reported in the paper. 

The automation of this repository relies on Python (including ngsolve/netgen, numpy, matplotlib, and other packages). 

For the installation of ngsolve/netgen, we recommend readers to follow the tutorial at https://ngsolve.org/downloads. There are also automatic scripts, you can follow the instructions based on your operating system:
* For Windows: (make sure conda is installed)
  - Open a terminal.
  - Run the following command: `ngsolve_setup.bat`

* For macOS/Linux: (conda will be automatically installed if not existed)
  - Open a terminal.
  - Run the following command: `source ngsolve_setup.sh`

## Structure of the repository

The repository is structured as follows:

* The implementation of a fully discrete format for high-order mass and energy-conserving methods for solving the nonlinear Schrödinger equation can be found in `./Packages/Collocation/BaseMethod.py`. This module contains the core classes that utilize implicit Newton iteration methods for solving the nonlinear system.

* The repository includes a collection of numerical examples located in the `./Numerical_Tests/` directory. These examples are organized as follows:

    - `1d_Standard_Soliton`:
        - 1D convergence tests (both spatial and temporal)
        - CPU time tests
        - Comparison of the proposed method with the standard Gauss collocation method
        - Tests on mass and energy conservation
        - Data collection on the number of Newton iterations and parameter iterations

    - `1d_Bi_Soliton`:
        - Long-term performance tests of the proposed method and comparison with the standard Gauss collocation method

    - `2d_Soliton`:
        - 2D convergence tests (both spatial and temporal)
        - Tests on mass and energy conservation
        - Data collection on the number of Newton iterations and parameter iterations

## How to reproduce or use the code

To reproduce the figures presented in the paper, you can follow the instructions in the `AutoFigGen` directory. All the required data are also included in the repository.

To regenerate all the data for the numerical experiments, navigate to the `Numerical_Tests` directory. These numerical examples can be performed by running the main functions `NLS_Collo_1d.py` and `NLS_Collo_2d.py` in terminal with the parameters provided in `.json` format. You can find the running commands collected in the scripts `experiments.bat` (for Windows) or `experiments.sh` (for Linux/Mac).

Please note that special attention should be given to the numerical examples in dimension 2, as they require a significant amount of system memory due to the large number of Degrees of Freedom (DoFs) involved. If you encounter errors related to insufficient memory (`err=not enough memory` or `NgException: PardisoInverse`), you can try reducing the parameters or decreasing the mesh resolution to test the code within the memory limits of your system.