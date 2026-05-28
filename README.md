# Generating the data for the paper
matlab NLCE_add.m

Update **Parameters** section of NLCE_add.m for generating results for different parameters.

# System parameters:
    $t,\ U=[U_{12}, U_{13}, U_{23}],\ \mu=[\mu_1, \mu_2, \mu_3],\ T$ .

# Generating observable expectation value for a set of system parameters
1. ED on every graph up to the given order. Stores the matrix elements of observables.
2. NLCE sum up to a given order on the thermal averages (computed using the matrix elements) of observables.

# NLCE graph data
NLCE data is generate using the standard NLCE algorithm of [Tang et al., Comp. Phys. Comm., 18, 3 (2013)](https://www.sciencedirect.com/science/article/pii/S0010465512003414).
### path for graphs
 `/NLCE/NLCE 2D graphs/graphsSimplified<n>.txt` 
 Stores the graphs of order n and all their subgraphs as a list of tuples $\{\{v_1,v_2\},\{v_3,v_4\},\cdots\}$ representing edges between verices, $v_i, v_j$.
### path for the coefficients
`/NLCE/NLCE 2D coefficients/coefficientsOfGraphs<n>.txt`
Stores the corresponding coefficients which includes contributions of the graph only to order n. 

For a graph, $c_n$ with n vertices, the coefficients are $L_{\mathcal{L}}(c_n)$, which counts the number of graphs symmetrically and topologically equivalent to $c_n$ in lattice $\mathcal{L}$ (square lattice in our case).

For a subgraph $s$, the coefficient is $\sum_{c_n,\ s.t.\ s\subset c_n} L_{\mathcal{L}}(c_n) \times$ (contributions of $s$ to $W_P(c_n)$). The weight, $W_P(c_n)= P(c_n) + \sum_{s\subset c_n}L_{c_n}(s)P(s)$, where $L_{c_n}(s)$ is computed iteratively from the definition, $W_P(c_n) = P(c_n) - \sum_{s\subset c_n}W_P(s)$.


# Fit parameters
`u=<U values>_fit_vals.txt` : The atomic limit and order 7 fitted temperature and center-of-trap chemical potential - $(T, \vec{\mu})/t$ with their respective fit errors. 

# Other files
`u=<U vals>_mus.csv` : Stores $\vec{\mu}(r) = \vec{\mu} - \frac{1}{2}m\omega^2r^2$, the chemical potentials away from the trap center with $m$ and $\omega$ being the mass of $^6\text{Li}$ and trap frequency respectively.

`u=<U vals>_T.csv` : The temperature fitted by NLCE order 7.

`NLCE_order=<order val>_u=<U val>_<observable_name>_<flavor info>.csv`: Stores the values of the observable computed at the corresponding NLCE order and U. If the values are flavor-dependent then files are named with the flavor info.
