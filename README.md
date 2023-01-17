# markPH-2missing

MATLAB code to analyze a simulated data set based on the CYD14 efficacy trial (Section 3.3 of the Web-based Supplementary Material for the JRSSC paper: Yanqing Sun (yasun@uncc.edu), Li Qi, Fei Heng, and Peter B. Gilbert (2021) Hybrid Approach for the Stratified Mark-Specific Proportional Hazards Model  with Missing Covariates and Missing Marks, with Applications to Vaccine Efficacy Trials

The folder contains:
	* sim_data.mat: a simulated dataset
	* MATLAB functions:
		+ main.m: main function
		+ aipw_imp1_strt.m: compute using the hybrid-MIEE method
		+ weight.m: calculate the inverse probability weight
		+ ksrmv.m: multivariate kernel smoothing regression

In main.m, we use load() to load the simulated dataset: load("sim_data.mat").
It includes input variables as shown below. 
	
INPUT:
------
	X: min{T,C}, where T is the failure time and C is the censoring time
	Zm: covariates subject to missingness
	Zc: covariates observed in all subjects
	V: the mark, also subject to missingness
	delta: I(T<=C), censoring indicator
	Zaux: the auxiliary variable predictive of Zm
	Vaux: the auxiliary variable predictive of missing marks
	strata: strata indicator. Strata jj has value jj-1. 
	
PARAMETER SETTINGS:
-------------------
	h: bandwidth
	L_n: number of nearest neighbors
	rep_n: number of imputations	
	Az_index: indicator for A_z, 1: with Az; 2: without Az
	Nbindex: indicator for A_v, 1: H=(T, Z1, Z2, A_v), 2: H=(T, Z1, Z2)	
	Parameters for tests: a, b, and aa;
	tt: the value of t used for estimating CIF rate
	zpd2: the value of Z1 used for estimating CIF rate
	
OUTPUT:
-------
* Estimates of covariate coefficients: eff_beta_hat
* Estimated standard deviation of estimators: sig_eff_beta
* Estimates of the CIF rates: 
	+ F10n: estimate of the CIF rate at the 10th percentile of Zm, 
	+ F50n: estimate of the CIF rate at the 50th percentile of Zm, 
	+ F90n: estimate of the CIF rate at the 90th percentile of Zm, 
* p-values: 
    + p_value_a: for the alternative H1a
	+ p_value_m: for the alternative H1m
	+ p_value_t2_a: for the alternative H2a
	+ p_value_t2_m: for the alternative H2m
	
COMPUTATION TIME:
-----------------
The required computation time for sample size 5000 is about 5 minutes running on the High Performance Computing Cluster at UNC Charlotte.
