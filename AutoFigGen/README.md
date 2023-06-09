## Automatic Generation of Figures

The figures presented in this paper are generated using scripts that utilize the data collected in the `../Numerical_Tests` directory.

To visualize all figures, follow the instructions based on your operating system:

* For Windows:
  - Open a terminal.
  - Run the following command: `AutoFig.bat`

* For macOS/Linux:
  - Open a terminal.
  - Run the following command: `source AutoFig.sh`

You can also generate specific figures by referring to the following list, along with their corresponding explanations.

## List of files and brief comment

### Numerical example 1:
In this example, we test the evolution of the standard soliton solution for T in [0,1]. Here is a description of each plotting script and their purposes:

* `Std1d_CPU_dt.py`: Plots the CPU time for the standard Gauss collocation part and the post-process correction part for different time steps. 
* `Std1d_Sol.py`: Plots the propagation of the modulus of the solution.
* `Std1d_En-Mass_Cp.py`: Plots the energy loss versus different time steps for two methods: the standard Gauss collocation method and the proposed post-processing correction method. This allows for comparison between the two methods.
* `Std1d_Spat_Conv.py`: Plots the spatial convergence of the proposed method.
* `Std1d_Temp_Conv.py`: Plots the temporal convergence of the proposed method.
* `Std1d_ME_Conserv.py`: Plots the mass and energy loss for the proposed post-processing correction method at the end of each time step, along with the number of Newton iterations and parameter iterations.

These scripts provide visualizations and comparisons for various aspects of the standard soliton solution.

### Numerical Example 2:

In this example, we test the long-term performance of the proposed method. We consider the bi-soliton solution, whose modulus is periodic in time.

* `BiSoli1d_LongT_ME_Err_Compare.py`: Plots the energy loss and the H1 error of the numerical solutions computed by the standard Gauss collocation method and the proposed method over a long time interval T in [0, 128].
* `BiSoli1d_LongTimeSol.py`: Plots the modulus of the solution for the time interval T in [0, 128].
* `BiSoli1d_ShortTimeSol.py`: Plots the modulus of the solution for some small time intervals.

### Numerical Example 3:

In this example, we test the convergence of the proposed method on a 2-dimensional problem.

* `2dSoli_ME_Conserv.py`: Plots the mass and energy loss for the proposed post-processing correction method at the end of each time step, along with the number of Newton iterations and parameter iterations.
* `2dSoli_Spat_Conv.py`: Plots the spatial convergence of the proposed method.
* `2dSoli_Temp_Conv.py`: Plots the temporal convergence of the proposed method.

