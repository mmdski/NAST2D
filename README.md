# NAST2D
# A Finite Volume Code for Fluid Flow

Downloaded from [https://people.math.sc.edu/Burkardt/cpp_src/nast2d/nast2d.html](https://people.math.sc.edu/Burkardt/cpp_src/nast2d/nast2d.html)

NAST2D is a C++ program which uses the finite volume method to model the behavior of an incompressible fluid in a 2D flow region.

The program treats the incompressible time-dependent Navier Stokes equations (velocity and pressure) as well as the heat equation. The program can handle free boundaries.

The program uses finite differences on a structured equidistant staggered grid for discretization of derivatives, with central and upwind (donor-cell) discretization of the convective term and an explicit time stepping scheme, Chorin's projection method. Free boundaries are treated with the MAC technique.

NAST2D writes a (binary) graphics dump file, typically with an extension of ".uvp", which contains the values of the solution. This file can be read in and plotted by IDL.

NAST2D can write out a "streak" file, typically with the extension ".str", containing a sequence of locations of marker particles, to give a striking illustration of the flow. The streak files can be very large, and their creation occasionally causes the program to crash!

NAST2D writes out a file containing the stream function values at the last step. This file is always called "psi_data.txt". The file contains ASCII data.

Usage:

nast2d input.par
where
input.par is a "parameter file" that specifies the problem to be solved and associated parameter settings.

## Reference:

1. ftp://ftp.lrz-muenchen.de/pub/science/fluiddynamics/cfd/NaSt2D the web site
1. Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
Numerische Simulation in der Stroemungsmechanik,
Vieweg 1995.
3. Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
Numerical Simulation in Fluid Dynamics,
SIAM, 1998,
This book can be ordered from the SIAM web site.
