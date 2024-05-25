# identify_grain
# example results on page 4 of 220318 NiO nanoparticle calculation.pdf
This script finds centers and sizes of nano-particles in an atomic force microscopy (AFM) image
The Matlab script named AFM0_ellipse.m takes in an AFM data file, such as 'La_0.csv'
then substracts off a non-physical polynomial background 
The script assigns indicators to small grains in the image 
and identifies the center and 2D area of each grain assuming an ellipital fit
combined with the height, the volume of each grain is obtained
to obtain the grain number counts as a function of grain volume

