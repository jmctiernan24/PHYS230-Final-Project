# PHYS230-Final-Project
The repository has the following documents (excluding codes):

1. Phys_230_final_presentation.pdf which holds my non-modified presentation, along with several additional slides at the end to showcase new results.

2. Phys_230_project.pdf which is a combination of my HW 3 submission (without modifications) and a brief summary describing the state of the project at the end with corresponding math.

Each .py file in the repository has the following purpose/instructions:

1. rk_fe_numba_combined.py 

A two component liquid-liquid phase separation model represented in Section 1 of the Phys_230_project.pdf file. It consists of two different time evolution methods: Forward Euler and 4th order Runge Kutta.

The instructions to run the code consist of: 

a. Choosing preferred model parameters on lines 21-23 

b. Choosing preferred system parameters on lines 25-30

c. Deciding whether you want graphs to be shown during the simulation ('yes' shows, everything else doesn't); replace line 32 with graph='yes' (can choose steps between each graph with line 33). Also possible to save each shown graph by uncommenting line 107.

d. Deciding if you want to keep track of simulation time throughout ('yes' keeps track, everything else doesn't); replace line 34 with tcls = 'yes' (steps between each calculation on line 36).

e. Pick which time evolution scheme you want, if line 96(97) is uncommented you are using Forward Euler(Runge Kutta). Can comment/uncomment respectively to change this.

f. Decide if you want to save the final figure, which can be done by uncommenting line 127.


