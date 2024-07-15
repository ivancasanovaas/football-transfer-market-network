# Codes

This folder contains the necessary codes to reproduce the results shown in the [report](../football-transfer-market-network.pdf). The programing language is FORTRAN95, as it is more efficient to compute properties and simulate network dynamics, but the plots are generated with the [MatPlotLib library](https://matplotlib.org/), from [Python](https://www.python.org/).

## Modules

4 modules were created to complete the requirements: 

0. **Random number generator:**  (`random_number_generator.f95`): This module generates a random number by initializing the system wwith a seed. It is necessary to include the file `r1279block.h` in the same folder of this module.

1. **Structural properties** (`structural_properties.f95`): This module contains functions & subroutines to read an undirected and unweighted network (given in edge list format), and calculate its structural properties (number of edges and nodes, list of degrees, list of neighbors and pointer, degree distributions, average nearest neighbors degree and clustering coefficient).

2. **Equilibrium models:**  (`equilibrium_models.f95`): This module contains two functions to maximally randomize the original undirected andunweighted network with two methods, mantaining the degree distribution. The random rewiring process (RW) & the configuration model (CM) are both presented.

3. **SIS dynamics:**  (`SIS_dynamics.f95`): This module reads a undirected and unweighted network given in edge list format and simulates the Susceptible Infected-Susceptible model with the Gillepsie algorithm. An initial number of infected nodes and an infection rate are required as input parameters.

## Running the main program

The functions and subroutines of all these modules are executed sequentially in the main program `main.f95`. It receives two inputs in the command line: the file name containing the edge list (`<edgelist_filepath>`) and the path for the output files (`<output_filespath>`). 
Running `main` is quite straightforward:
```
# Command line
gfortran -c random_number_generator.f95
gfortran -c structural_properties.f95
gfortran -c equilibrium_models.f95
gfortran -c SIS_dynamics.f95
gfortran -c main.f95
gfortran -o main main.o random_number_generator.o structural_properties.o equilibrium_models.o SIS_dynamics.o
./main <edgelist_filepath> <output_filespath>
```
The script `exe.sh`  contains these command lines, with the paths already specified. To execute it:
```
# Command line
./exe.sh
```
## Format of input files

The structure of the network to be embedded is passed to the program via a file containing its edgelist (one link per line). The edgelist file consists in a simple text file with the following convention

```
[name of node1]  [name of node2]  [remaining information will be ignored]
[name of node2]  [name of node3]  [remaining information will be ignored]
[name of node4]  [name of node5]  [remaining information will be ignored]
# comments can be inserted between links
[name of node5]  [name of node6]  [remaining information will be ignored]
[name of node7]  [name of node6]  [remaining information will be ignored]
...
```

Note that the nodes' name will be imported as `std::string` and can therefore be virtually anything as long as they do not include white spaces (i.e., there is not need for the nodes to be identified by contiguous integers).

**IMPORTANT**: Lines beginning with "#" (comments)are not ignored, so drop them out. File must contain solely the edge list.

**IMPORTANT**: this class only considers **simple undirected** networks **without self-loops**. Any multiple edges (e.g., if the graph is originally directed) or self-loops will be ignored.

**IMPORTANT**: in the actual version of the code, the network must have **only one component**.

## Figures

To reproduce the figures in the report, open the Jupyter Notebook `plots.ipynb` and run it. (Results text files must be in the specified output files path.)