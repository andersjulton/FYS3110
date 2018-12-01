

# Notes 

- I think the one dimensional impementation works now 

- In main1D.cpp: the exact solution is set to one (line 33), since it is not implemented yet

- How can we get rid of pragma once? (What can we replace it with?)

- something is strange with time steps in main1D! See output: 
----------------------------------
t = 0.50, t = 0.05
dt = 0.000050, t = 0.500000 
101
time steps = 10001
dt = 0.000050, t = 0.050000 
101
time steps = 1001
dt = 0.005000, t = 0.500000 
11
time steps = 100
dt = 0.005000, t = 0.050000 
11
time steps = 10
---------------------------------
One too few time steps when dt = 0.005 ??

