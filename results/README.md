# Results

This folder contains the results obtained running the script `exe.sh` from the codes folder. It calculates the structural properties of the original network, compared them with equilibrium models, and finally runs the SIS dynamics. The files with the results are in the `dat` folder, and the plots of these results are in the `figs`folder. 


## Output text files (`dat`)

When executing the main program `main.f95`, it outputs several files divided in two folders: `properties` and `SIS`. On the one hand, the extracted network structural properties are:

- Printing to screen the number of nodes `N`, the number of edges `E` and the average node degree `<k>`.
- `degrees.txt`: list of degrees for each node
- `neighbors.txt`: list of neighbors for each node
- `pointer.txt`: list with the positions, in V, of the first and last neighbors for each node
- `ddd.txt`, `cdd.txt` and `ccdd.txt`: direct, cumulative and complementary cumulative degree distributions 
- `knn.txt` and `c.txt`: average nearest neighbors degree and clustering coefficient as a function of the degree for the original network
- `knn_rw1.txt` and `c_rw1.txt`: average nearest neighbors degree and clustering coefficient as a function of the degree for 1 RW network
- `knn_rw100.txt` and `c_rw100.txt`: average nearest neighbors degree and clustering coefficient as a function of the degree for 100 RW networks (averaged)
- `knn_cm1.txt` and `c_cm1.txt`: average nearest neighbors degree and clustering coefficient as a function of the degree for 1 CM network
- `knn_cm100.txt` and `c_cm100.txt`: average nearest neighbors degree and clustering coefficient as a function of the degree for 100 CM networks (averaged)

On the other hand, the Susceptible-Infected-Susceptible dynamics provides 3 interesting magnitudes to measure as a function of 2 parameters (initial number of infected nodes, and the infection rate):

- `SIS-life_time`: 
- `SIS-rho_st`: average steady-state prevalence
- `SIS-rho_t`: empirical prevalence as a function of time

## Output figures (`figs`)

The output figures can be modelled from the Jupyter Notebook `plots.ipynb` in a simple manner.