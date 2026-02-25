# Angular_Correlations

# This program is under construction in order to replace a fortran 77 angular distribution/ mixing ratio determination calculator used a FSU's Nuclear Accelerator laboratory. 

# Reading in detector width, radius, distance from the radiative source, and the gamma-ray energy of interest, the program calculates the QDk geometric attenuation coefficients. 

# The program also asks for the magnetic substate probability distribution - or sigma - and the j1 and j2 values of the states from the gamma-ray decay from state 2 to 1. Using this information the coupling of states via clebsh - gordon coefficients in which the quantum mechanical restrictions given the angular momentum of the system can be used to understand the states decay probabilities following the reaction of interest. 

# This program then reads in a [].csv, the format must match the test file in the directory. This file is a file of intensities as a function of angles or different detectors for a given decay energy. The data is then used for the legendre polinomial f(x) = A0 P0(x) + A2 P2(x) + A4 P4(x), gaussian elimination is used to return A0, A2, A4. These values are used to weight the chi-squared minimization comparison from theory to experiment. The minimization is key to locating the electromagnetic mixing ratio, as a function of angular momentum and coupling quantum mechanics.

# Finally reading from Rose and Brink, the program reads the first 2 Racah coefficient sets, 6 total, depending on j2 and j1. The program can be modified to use more than 3 angles, and for  higher multipoles of the angular distribution. Using these coefficients, theoretical A2 and A4 values can be compared, illuminating the mixing ratio. 

/*
======================== README (inline) ========================

Overview
--------
This program computes and visualizes gamma-ray angular distributions and a chi-squared scan over mixing ratio δ
in a single ROOT GUI window.

GUI Inputs (user-editable)
--------------------------
- j1, j2        : spins (integer or half-integer)
- Eγ (keV)      : gamma-ray energy in keV
- Sigma         : magnetic substate width σ (0 for perfect alignment)
- Angular file  : experimental data file (Browse)

Buttons
-------
- Run / Compute : runs all physics and updates both plots
- 1             : show angular distribution view
- 2             : show chi-squared scan view
- 3             : close GUI and exit program
- Redraw        : redraw view
- Reset Zoom    : reset autoscale

Input File Format
-----------------
Experimental angular file must have at least 3 numeric columns:

  theta(deg)   Y   Yerr

Delimiters allowed: comma, spaces, tabs (mixed is OK).
Lines starting with '#' are comments.

Plot Behavior
-------------
- Angular view always displays theta from 0 to π (full domain), regardless of data range.
- Fit overlay (red) is computed from A0/A2/A4 and drawn across [0, π].
- Chi2 view displays log(χ²) vs atan(δ).

Terminal Output
---------------
Each Run prints:
- A0, A2, A4 and normalized a2=A2/A0 and a4=A4/A0
- QD2 and QD4 with the geometry (R,D,T) and energy used

Racah Coefficients
------------------
Rk coefficients are computed internally using Wigner 3j/6j (no lookup table is used).
Selection rule enforced: |j1 - j2| must be an integer (electromagnetic transitions).

Build
-----
g++ -std=c++17 AD3.cxx `root-config --cflags --glibs` -lGui -o AD3

If you see ROOT attribute warnings, they are benign; to silence:
g++ -std=c++17 AD3.cxx `root-config --cflags --glibs` -lGui -Wno-attributes -o AD3

===============================================================
*/
