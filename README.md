# PHYS230-Final-Project
The repository has the following documents (excluding codes):

1. Phys_230_final_presentation.pdf which holds my non-modified presentation, along with several additional slides at the end to showcase new results.

2. Phys_230_project.pdf which is a combination of my HW 3 submission (without modifications) and a brief summary describing the state of the project at the end with corresponding math.

Each .py file in the repository has the following purpose/instructions:

Note: Using Runge Kutta means the simulation can handle higher dt.

1. rk_fe_numba_combined.py 

A two component liquid-liquid phase separation model represented in Section 1 of the Phys_230_project.pdf file. It consists of two different time evolution methods: Forward Euler and 4th order Runge Kutta. Note: density array ranges from -1 to 1.

The instructions to run the code consist of: 

a. Choosing preferred model parameters on lines 21-30. Be careful deciding initial density and variance, if every component of the array is not within -1 or 1 the code will run infinitely.

b. Deciding whether you want graphs to be shown during the simulation ('yes' shows, everything else doesn't); replace line 32 with graph='yes' (can choose steps between each graph with line 33). Also possible to save each shown graph by uncommenting line 107.

c. Deciding if you want to keep track of simulation time throughout ('yes' keeps track, everything else doesn't); replace line 34 with tcls = 'yes' (steps between each calculation on line 36).

d. Pick which time evolution scheme you want, if line 96(97) is uncommented you are using Forward Euler(Runge Kutta). Can comment/uncomment respectively to change this.

e. Decide if you want to save the final figure, which can be done by uncommenting line 127.



2. modified2_final.py

This code is very similar to 1., except for the slightly modified system being simulated: it models a two component mixture based off the Flory-Huggins free energy density (explained more in Section 6.2 of Phys_230_project.pdf). Note: density array ranges from 0 to 1.

Instructions to run code:

a. Choosing preferred model parameters on lines 23-32. Be careful deciding initial density and variance, if every component of the array is not within -1 or 1 the code will run infinitely.

b. Deciding whether you want graphs to be shown during the simulation ('yes' shows, everything else doesn't); replace line 34 with graph='yes' (can choose steps between each graph with line 35). Also possible to save each shown graph by uncommenting line 105.

c. Deciding if you want to keep track of simulation time throughout ('yes' keeps track, everything else doesn't); replace line 36 with tcls = 'yes' (steps between each calculation on line 38).

d. Pick which time evolution scheme you want, if line 94(95) is uncommented you are using Forward Euler(Runge Kutta). Can comment/uncomment respectively to change this.

e. Decide if you want to save the final figure, which can be done by uncommenting line 123.



3. three_component.py

This code models the mixture of three components, as described in Section 6.3 of Phys_230_project.pdf. Compared to the previous two, this code is functionally simpler only outputting a final graph for two of the components. Note: this code takes much longer than 1 and 2.

Instructions to run code:

a. Choosing preferred model parameters on lines 21-31. Be extra careful choosing the initial densities and variances; every component of the two arrays added together must be between 0 and 1, while each array ranges from 0 and 1. A smaller time value is required in comparison to previous models.

b. Decide if you want to save the final figures, which can be done by uncommenting lines 106 and 112.


4. old_two_comp.py

An old version of code 1 that uses numpy exclusively (and only includes Runge Kutta evolution), leading to a very slow code. It is included in the repository to show how much more efficient it is to use numba as opposed to numpy.roll.

Instructions to run code:

a. Choose parameters for system on lines 50-57. 

b. Decide whether you want to keep final image by uncommenting line 74.

5. time_graph_gen.py

The only use for this code is to plot the results when determining the simulation time from code 1: it can be used to compare the simulation time between the Runge Kutta and Forward Euler methods.
