# Double-Pendulum
Researching the motion of a double pendulum with no air resistance.

This project was initiated and completed in 2017, and presented during the 2017 Bremen Summer School of Mathematics.

## Abstract
The abstract contains all necessary information to understand this project, and the highest resolution picture of the pretty graph we get as a result. 

## Images
The images are in ASCII ppm format. You can use GIMP or Photoshop to view them. I also recommend the website http://paulcuth.me.uk/netpbm-viewer/, which is a simple netBPM viewer.
The filename format is: image[x-res]x[y-res]\_[max number of steps]\_[enery error in sci. notation]\_

So, for example, image300x300_250_1e-1_ is a 300x300 ppm image, where for each pixel the pendulum motion is generated until either 250 steps or a pendulum-flip occurs. The step size is 1\*10^-1, or .1 seconds. 

## Solution method
There are two parameters in the double pendulum system - theta and phi, the two angles of the rods. Each of these are governed by a second-order differential equation (acceleration is due to gravity only). So, we have two paired DE's. These two DE's are split into two first-order DE's, resulting in a total of 4 paired first-order DE's. Each of these are solved using the Runge-Kutta numerical method to give us the curve the pendulum follows. The resulting curve is traced until either a) the maximum time for simulation is reached, or b) either angles have done a greater-than-360 degree rotation (a flip). We record the time it took for a flip to occur, and plot the result over the domain of choice (default pi\*pi region), with color representing the length of time. The resulting function looks smooth in some areas, but poorly-behaved and fractal-like in othe areas. An elliptical shape can be derived such that within it, no initial states of the system will lead to a flip because there is not enough energy to do so. Outside of this shape, predicting flipping time is no longer analytically feasible.  



