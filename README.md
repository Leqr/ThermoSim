# ThermoSim
Thermocline Thermal Energy Storage Simulation

Project for the graduate course "Fundamentals of Computational Fluid Dynamics Methods" at ETH Zürich.

In order to run the program you need a recent version of C++, CMake and Python 3.x with numpy and matplotlib.

In order to run the program the following need to be done :

1. Open a terminal and go inside the build folder (mkdir build then cd build if the folder is not present).
2. cmake ..
3. make

Now the executable CFDSim is built. The program can be run in different modes by passing a different argument (integer) when calling the executable. The C++ code calls python for plotting automatically, only the integer that specify the mode needs to be changed :

./CFDSim 0 runs the program in MMS mode by reaching the steady state of the MMS problem with $n=1$ and plots the cosine against the numerical solution.

./CFDSim 1 runs the uncoupled OVS and plots the convergence graphs, by default $Pe = 1000$ and $n = 1$, these values can be changed in line 142 of main.cpp.

./CFDSim 2 runs the coupled OVS and plots the convergence graphs,by default $Pe = 1000$, this value can be changed in line 150 of main.cpp.

./CFDSim 3 compares the numerical solution with the analytical result by fetching the analytic solution in sol-exact-5.00000E+03.dat.

./CFDSim 4 perform the design study with diameter = 8 (can be changed in setup.txt, needs also to change the values of $u_f$, height, $h_{v,s}$ and $h_{v,f}$ with the one corresponding to the wanted diameter written at the bottom of the file) and prints the calculated exergy factor, capacity factor and temperature values, also plots a time dependent simulation of the numerical solution.

All the parameters of the simulation are defined in the setup files which are called when launching the program. Each file is relevant for a specific simulation.

setup.txt is used for the design study (mode 4).

setuppr5.txt for comparing the numerical and analytical solution (mode 3).

setupOVS.txt for the OVS (mode 1 and 2).

setupMMS.txt for comparison with the cosine (mode 0).

Overall the main.cpp file performs the initialisation of a Simulator and Exporter object as well as fetching the parameters and mode of the simulation.
The Simulator is launched by the main either with Simulator.simulate(...) which performs a standard cycle based simulation or with Simulator.OVS(...) which performs the order verification study.
The Exporter object is used as an interface between the .txt output files (which are created and stored at runtime in the build folder) and the program.
