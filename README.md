# MATLAB_spin-dependent_Fermi-Hubbard
Matlab package to compute thermodynamic observables of spin-dependent Fermi-Hubbard model using numerical linked-cluster expansion (NLCE) on a 2D square lattice , and fit to experimental data.

&copy; Sohail Dasgupta, Haotian Wei and Kaden R. A. Hazzard

If you use this code, please cite Mongkolkiattichai et al. *Quantum gas microscopy of three-flavor Hubbard systems*, (In final stages of preparation for Science), 2026.

# Acknowledgements
This work was supported in part by the NOTS cluster operated by Rice University's Center for Research Computing (CRC).

# Contents
1. [What is this package?](#what-is-this-package)
2. [Getting Started](#getting-started)
3. [Usage](#usage)
4. [Key functions](#key-functions)
5. [License](#license)
6. [Contributors](#contributors)

# What is this package?
This package generates the values of the thermodynamic observables such as particle density, density of pairs, and nearest-neighbor density-density correlation functions of the three-flavor spin-dependent Fermi-Hubbard model for any temperatures and (spin-dependent) chemical potentials up to a given order of the site-expansion numerical linked-cluster expansion (NLCE). 

The inputs required are the optical lattice experimental data (observables as a function of distance from the optical trap center), NLCE graphs (provided as a cell of edges $\{\{v_1,v_2\},\{v_3,v_4\},\cdots\}$) and the coefficients for up to the maximum order of interest, and the Fermi-Hubbard parameters. 

The outputs are the best-fit value of the experimental temperature and the center-of-trap chemical potentials for every spin flavor, and the observable values as a function of the local chemical potential across a radial cut of the trap.

# Getting Started

## Software Dependencies
The following software and its corresponding toolboxes need to be installed prior to using this repository.

 - `MATLAB` : R2025b or higher. 
    - `Optimization Toolbox` : Version 25.2
    - `Symbolic Math Toolbox` : Version 25.2 
    - `Parallel Computing Toolbox`: Version 25.2. (Optional. Needed for `ED_solver_parallel` only) 

Older versions have not been tested by the authors.

## Directory Structure
 **All paths in this document are relative to the directory of the local repository.**

Make sure the following directories exist before running the script.
- `../data/csv_files/N=3/` : The final outputs are generated here as .csv files.
- `../mat_files/ED_mat_files/N=3/`: The intermedidate matrix elements are stored as .mat files.
- `../mat_files/NLCE_mat_files/N=3/` : The thermal averages for every graph is stored as .mat files.
- `../experimental_data_U13_<u13>_U12_<u12>_U23_<23>` : The experimental data. 

These dependencies can be changed from `NLCE_add.m`.

# Usage

## Running the script
Make sure the above software dependencies and directory structure is met, NLCE graph data and experimental data are provided.

Then use run the NLCE_add script either from command line or from MATLAB GUI
### Command Line
```
matlab -batch NLCE_add.m
```
Note that this may require adding MATLAB to system path.
### MATALB GUI
Open the NLCE_add.m using MATLAB and run. 

## Input
### Experimental data
| Filepath | Description |
| -------- | ----------- |
| `../experimental_data_U13_<u13>_U12_<u12>_U23_<23>/*<observable>*.txt` | Observable as function of $r$ (radial distance from the center of trap). The code ignores the first line as header. The first column is read as the $r$-values. The next three columns are the spin-dependent values of the observables. The last three columns are the corresponding experimental error. |
| `../experimental_data_U13_<u13>_U12_<u12>_U23_<23>/Parameters.txt`| Experimental parameters-- $s$ in the units of the recoil energy, $t$ in Hz, Lattice confinement in horizontal direction Hz, Atom number, $U/t$-- and zeroth order high-T series fit values.|

### NLCE graph parameters
NLCE data may be generated from the NLCE algorithm of [Tang et al., Comp. Phys. Comm., 18, 3 (2013)](https://www.sciencedirect.com/science/article/pii/S0010465512003414). 

A thermodynamic property of a translatinally invariant graph, $\mathcal{L}$ can be written as $P(\mathcal{L})/N = \sum_{n=1}^\infty\sum_{c_n} l(c_n) W_P(c_n)$, where $N$ is the lattice-size and $l(c_n)$ is the number of ways the cluster $c_n$ of size $n$ can be embedded in the lattice upto translational invariance, and $W_P(c_n)$ is the weight of the cluster, defined iteratively as $W_P(c_n) = P(c_n) - \sum_{s\subset c_n}W_P(s)$.

This is re-written as $P(\mathcal{L})/N = \sum_n \left(\sum_{c_n}k(c_n)P(c_n) + \sum_{m=1}^{n-1}\sum_{c_m} k^{(n)}(c_m) P(c_n)\right)$. Here, $P(c_n)$ is the property value on graph $c_n$ computed using exact digonaliztion, $k(c_n)$ is the number of ways to embed the graph $c_n$ in the lattice and $k^{(n)}(c_m)$ counts the contribution of $c_m$ as a subgraph to order $n$ graphs. 

For this repository, data for **2D square lattice** provided. 

| parameter | Path| Description |
| ------| ---- | ---------- |
| graphs | `/NLCE/NLCE 2D graphs/graphsSimplified<n>.txt` | Stores the graphs of order $n$, $c_n$ and all their subgraphs $c_m$ for $m=1,\cdots,n-1$ and $k^n(c_m)\neq 0$ as a list of tuples $\{\{v_1,v_2\},\{v_3,v_4\},\cdots\}$ representing edges between verices, $v_i, v_j$. |
| coefficients | `/NLCE/NLCE 2D coefficients/coefficientsOfGraphs<n>.txt` |Stores the corresponding non-zero coefficients of the graphs, $k(c_n)$ and their subgraphs, $k^{(n)}(c_m)$ for $m=1, \cdots, n-1$).

## Simulation parameters
The simulation parameters can be modified in the `Parameters` section of `NLCE_add.m`
| Variables | Description |
| --------- | ----------- |
| t | Hubabrd $t$. This is set to $1$ to normalize other values to it. |
| u | Spin-dependent Hubbard $U$. For three-flavors $U = [U_{12},U_{13}, U_{23}]$.| 
| m | Number of spin flavors.| 
| mfermion | The mass of the fermion used in the experiment in kg.|
| trap_laser | The frequency of the trap laser in m.|
| h | The Planck constant in SI units.| 
| order_max | The maximum order of NLCE to be computes.|

## Output
All the following files are stored under `../data/csv_files/N=3/` (*change in NLCE_add.m*).

| Files | Description |
| ----- | ----------- |
|`u=<U values>_fit_vals.txt`  |  The atomic limit and order 7 fitted temperature and center-of-trap chemical potential - $(T, \vec{\mu})/t$ with their respective fit errors. |
| `u=<U vals>_mus.csv` | Stores $\vec{\mu}(r) = \vec{\mu} - \frac{1}{2}m\omega^2r^2$, the chemical potentials away from the trap center with $m$ and $\omega$ being the mass of $^6\text{Li}$ and trap frequency respectively. |
| `NLCE_order=<order val>_u=<U val>_<observable_name>_<flavor info>.csv` | Stores the values of the observable computed at the corresponding NLCE order and $U$. If the values are flavor-dependent then files are named with the flavor info. |

## Example use
Fit the density vs $r$ experimental data in `sample_exp_data` to NLCE order 5.

- Open the NLCE_add.m and make sure the `Parameters` section (the first section)  has the following values - 
```
%% Parameters
% Fermi-Hubbard parameters
t = 1; u = [7.9, 13.7, 1.8]; % U12, U13, U23
m = 3;

% Experimental parameters
mfermion = 6*1.66054e-27; %mass of a Li6 atom in kg
trap_laser = 752e-9; % Frequency of the trapping laser (in nm)
h = 6.62607015e-34; % Planck's constant in J-s

% Maximum NLCE order (site-expansion)
order_max = 5;
```
- Open a command line from the local repository directory and exeute the following commands.
```
cp sample_exp_data/. ../.
matlab -batch NLCE_add
```
- Check the best-fit values. 
```
vi ../data/csv_files/N=3/u=7p9,\ 13p7,\ 1p8_fit_vals.txt
```
It should read 
```
Atomic Limit Fit
T/t = 2.1052 +/- 0.0664
mu0/t = 5.3669 +/- 0.1169	4.2317 +/- 0.0993	4.3785 +/- 0.1001	



NLCE order: 5
T/t = 1.8000 +/- 0.0768
mu0/t = 5.4797 +/- 0.1201	4.1914 +/- 0.1026	4.4203 +/- 0.1040	
```

# Key functions
The following tables show the key functions that for the NLCE algorithm.

| Functions | Description |
| --------- | ----------- |
| ED_solver | Exact diagonalization on a graph at $\vec{\mu}=(\mu_1,\mu_2,\mu_3)=(0,0,0)$. Outputs the matrix elements of the observables in the eigenbasis. |
| thermal_average |  Given the matrix elements, performs the thermal average. |
| NLCE_sum | Performs NLCE sum on themal observables up to a given order, for a given list of $T$ and $\vec{\mu}$. |
| obs_vs_r_fitting_function | Using a given experimental data set of an observable \(Eg. density, doublon\) as a function of distance from the trap center $r$, finds the best fit values of $T$ and $\vec{\mu}(r=0)$ for a given order of NLCE. |  
| celltuples | Helps generate basis states from spin configurations.|
| key_gen | Generates unique key for every graph. |

# License

The code in this repository is licensed under the [MIT License](LICENSE).

This repository contains code supporting a forthcoming research publication.


# Contributors
Sohail Dasgupta, Haotian Wei and Kaden Hazzard