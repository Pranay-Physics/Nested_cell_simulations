<h1> Nested herriott cell concentric ring simulation GUI </h1>

This project is to simulate the formation of a spot pattern in a nested herriot cell. 
For this simulation I have chosen 3 physical parameters that we can control, Radius of cuvature of outer mirror R1, Radius of curvature of inner mirror R2, and the coaxial mirror seperation d.

For use of this program, just run the nested_cell_slider program, it contains all of the functions and GUI code. The horizontal sliders control the value of the 3 parameters.
At the bottom left you can see the entry and exit state which are of the form [x, x', y, y'], As well as the value of B/A, time in cavity and the final spot size.
note that spot size displayed aren't the actual spot sizes that will be observed, only the relative size between the points will be accurate. this is due to slow plotting of circles compared to scatter plot in python

You may change the range of values each of the parameters can take and the resolution by which they change in the slider by adjusting the values in the self.{parameter}_slider 
lines in the init function of the class (starting line 116). 

You may wish to simulate a different number of concentric rings, to do that you will need to change the value of n (line 94), initial distance d (line 102) to a reasonable value and the range of values of d
(line 116). Reasonable values of d for different n are as follows. 

| n   | d   |
|------------|------------|
| 2      | 1160   |
| 3     | 535     |
| 4      | 320     |
