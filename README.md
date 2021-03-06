The aim of this repository is to share the code of the reluctance network model presented in the Energies journal article: **Open Circuit Performance of Axial Air Gap Flux Switching Permanent Magnet Synchronous Machine for Wind Energy Conversion: Modeling and Experimental Study.**

Files description:

`main`: The main project file, run it to generate the results needed.

`move-rotor`: Complementary to *main* and can't be run alone.

`swap_matrix_x_720`: Complementary script that manipulates the reluctance values accroding to the new rotor position. Without it the *move-rotor* script is useless.

`CoggingForceCalculationMeanRadius`: Calculate the electromagnetic torque through a numerical approximation based on Maxwell's stress tensor. The user is free to use whatever method he wants.

Because the project may take some time to run for a large number of iterations, the previously generated *.mat* files are results ready for viewing and post processing.
