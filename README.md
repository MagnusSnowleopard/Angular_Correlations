# Angular_Correlations

# This program is under construction in order to replace a fortran 77 angular distribution/ mixing ratio determination calculator used a FSU's Nuclear Accelerator laboratory. 

# Reading in detector width, radius, distance from the radiative source, and the gamma-ray energy of interest, the program calculates the QDk geometric attenuation coefficients. 

# The program also asks for the magnetic substate probability distribution - or sigma - and the j1 and j2 values of the states from the gamma-ray decay from state 2 to 1. Using this information the coupling of states via clebsh - gordon coefficients in which the quantum mechanical restrictions given the angular momentum of the system can be used to understand the states decay probabilities following the reaction of interest. 

# This program then reads in a [].csv, the format must match the test file in the directory. This file is a file of intensities as a function of angles or different detectors for a given decay energy. The data is then used for the legendre polinomial f(x) = A0 P0(x) + A2 P2(x) + A4 P4(x), gaussian elimination is used to return A0, A2, A4. These values are used to weight the chi-squared minimization comparison from theory to experiment. The minimization is key to locating the electromagnetic mixing ratio, as a function of angular momentum and coupling quantum mechanics.

# Finally reading from Rose and Brink, the program reads the first 2 Racah coefficient sets, 6 total, depending on j2 and j1. The program can be modified to use more than 3 angles, and for  higher multipoles of the angular distribution. Using these coefficients, theoretical A2 and A4 values can be compared, illuminating the mixing ratio. 

# This program has its own graphical interface. Using X11 graphics, there is a functional terminal based menu that will promt and execute the 2 plots of interest. Option 1 displays the chi-squared minimization. Option 2 displays the Angular distribution overlayed with a legendre polinomial fit and respective error bars. Option 3 closes the program. 

# The GUI can zoom in by dragging and letting go, you can draw with left click, with right click the coordinates of the point is printed to the terminal. With space-bar you can unzoom the veiw. 
