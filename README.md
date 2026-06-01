# MATLAB Spin-dependent Fermi-Hubbard
Matlab package to compute thermodynamic observables of spin-dependent Fermi-Hubbard model using numerical linked-cluster expansion (NLCE).

&copy; Sohail Dasgupta, Haotian Wei and Kaden R. A. Hazzard

If you use this code, please cite Mongkolkiattichai et al. Quantum gas microscopy of three-flavor Hubbard systems, (In final stages of preparation for Science), 2026.

# Acknowledgements
This work was supported in part by the NOTS cluster operated by Rice University's Center for Research Computing (CRC).

# Contents
1. [What is this package?](#what-is-this-package)
2. [Usage](#usage)
3. [Code structure](#code-structure)
4. [License](#license)
5. [Contributors](#contributors)

# What is this package?
This package generates the values of the thermodynamic observables of the three-flavor Fermi-Hubbard model with spin-dependent interaction potentials for any temperatures and (spin-dependent) chemical potentials for up to a given order of the site-expansion numerical linked-cluster expansion (NLCE).

# Usage

## Running the script
`matlab NLCE_add.m`

Note : Update **Parameters** section of NLCE_add.m for generating results for different parameters.

## Input
## Experimental data
Directory : `../experimental_data_U13_<u13>_U12_<u12>_U23_<23>`
Relevant filenames : `*<observable>*.txt` - Observable as function of $r$ 


## NLCE graph data
NLCE data may be generated from the NLCE algorithm of [Tang et al., Comp. Phys. Comm., 18, 3 (2013)](https://www.sciencedirect.com/science/article/pii/S0010465512003414). 

### path for graphs
 `/NLCE/NLCE 2D graphs/graphsSimplified<n>.txt` 
 Stores the graphs of order n and all their subgraphs as a list of tuples $\{\{v_1,v_2\},\{v_3,v_4\},\cdots\}$ representing edges between verices, $v_i, v_j$.
### path for the coefficients
`/NLCE/NLCE 2D coefficients/coefficientsOfGraphs<n>.txt`
Stores the corresponding coefficients which includes contributions of the graph only to order n. 

For a graph, $c_n$ with n vertices, the coefficients are $L_{\mathcal{L}}(c_n)$, which counts the number of graphs topologically equivalent to $c_n$ in lattice $\mathcal{L}$ (square lattice in our case).

For a subgraph $s$, the coefficient is $\sum_{c_n,\ s.t.\ s\subset c_n} L_{\mathcal{L}}(c_n) \times$ (contributions of $s$ to $W_P(c_n)$). The weight, $W_P(c_n)= P(c_n) + \sum_{s\subset c_n}L_{c_n}(s)P(s)$, where $L_{c_n}(s)$ is computed iteratively from the definition, $W_P(c_n) = P(c_n) - \sum_{s\subset c_n}W_P(s)$.

## Output
Directory: `../data/csv_files/N=3/` (*change in NLCE_add.m*)
### Fit parameters
`u=<U values>_fit_vals.txt` : The atomic limit and order 7 fitted temperature and center-of-trap chemical potential - $(T, \vec{\mu})/t$ with their respective fit errors. 

### Other files
`u=<U vals>_mus.csv` : Stores $\vec{\mu}(r) = \vec{\mu} - \frac{1}{2}m\omega^2r^2$, the chemical potentials away from the trap center with $m$ and $\omega$ being the mass of $^6\text{Li}$ and trap frequency respectively.

`u=<U vals>_T.csv` : The temperature fitted by NLCE order 7.

`NLCE_order=<order val>_u=<U val>_<observable_name>_<flavor info>.csv`: Stores the values of the observable computed at the corresponding NLCE order and $U$. If the values are flavor-dependent then files are named with the flavor info.

# Code Structure
## Key functions
ED_solver : Exact diagonalization on a graph at $\vec{\mu}=(\mu_1,\mu_2,\mu_3)=(0,0,0)$. Outputs the matrix elements of the observables in the eigenbasis. 

thermal_average : Given the matrix elements, performs the thermal average.

NLCE_sum : Performs NLCE sum on themal observables up to a given order, for a given list of $T$ and $\vec{\mu}$.

obs_vs_r_fitting_function : Using a given experimental data set of an observable \(Eg. density, doublon\) as a function of distance from the trap center $r$, finds the best fit values of $T$ and $\vec{\mu}(r=0)$ for a given order of NLCE.

# License

# Contributors
Sohail Dasgupta, Haotian Wei and Kaden Hazzard