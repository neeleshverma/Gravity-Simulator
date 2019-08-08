# Gravity-Simulator
Graphical simulation of particles in 2-D space 

This is a graphical simulation of particles in 2-D space using **Barnes-Hut algorithm**.  It uses quad tree data structure to recursively divide the 2D space into 4 equal parts (4 rectangles) till each part contains at max one particle. 

main.rkt is the main file that contains the code to handle simulation and other things.
helper_functions.rkt contains the helper functions that are being used by main.rkt 
drawing-routing.rkt contains the code to draw (present) the simulation.
