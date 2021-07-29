# Solver stability tests

This set of MATLAB scripts let's you test the numerical stability of three ice-dynamics solvers (DIVA, L1L2-SIA, Hybrid) for a range of grid resolutions. The main script to run is `sample_stability.m`. The user specifies a solver to test and the test case of interest ('strong' or 'weak'). The script will then calculate the analytical 1D solution for the given solver and test case, followed by numerical verfication using the 1D model. If `run_sample=false`, then the script will load results that were previously saved to a .mat file. Finally, a figure will be produced that compares the numerical results to the analytical solution. 