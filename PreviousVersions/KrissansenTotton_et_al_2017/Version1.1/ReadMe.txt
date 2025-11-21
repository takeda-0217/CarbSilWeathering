Version 1.1

This set of python scripts runs our geological carbon cycle model for the last 100 Ma and plots selected outputs alongside proxy data from the literature. The model is described in  J. Krissansen-Totton and D. C. Catling (2017) "Constraining climate sensitivity and continental versus seafloor weathering using an inverse geological carbon cycle model", Nature Communications, DOI: 10.1038/ncomms15423. With the appropriate choice of parameters, this code can reproduce Fig. 2, 3, and 4 in the main text.

As a matter of courtesy, we request that people using this code please cite Krissansen-Totton et al. (2017). In the interest of an "open source" approach, we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author (joshkt@uw.edu)

Updates this version: Code for performing inverse analyses has been added. 

REQUIREMENTS: Python, including numpy, pylab, and scipy modules. Inverse model calculations require emcee and corner.

HOW TO RUN CODE FORWARD MODEL:
(1) Put all the python scripts and the folder geochemical proxies in the same directory, and ensure python is working in this directory.
(2) Open Main_code_forward_modeling.py and input desired parameter ranges and number of iterations (default parameters reproduce Fig. 2 in main text)
(3) Run Main_code_forward_modeling.py. Code will output model confidence interval alongside proxy data. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXPLANATION OF CODE STRUCTURE:

%% Main_code_forward_modeling.py
This script provides the shell to repeatedly call the forward model and plot the output distributions. 
Ranges for uncertain parameters and number of forward model calls (ÒiterationsÓ) can be modified by the user. 
Modifying anything else may produce errors in the code. Once parameter ranges have been defined, 
the script calls the forward model once and uses the dimensions of the outputs to define an output array to store all future outputs. 
The forward model is then called many times (equal to ÒiterationsÓ) and parameter ranges are randomly sampled for each forward model call. 
Forward model calls resulting in errors or non-physical outputs are discarded (the code may still produce error messages from these discarded outputs). 
The remainder of the script calculates 90% confidence intervals for model outputs based on the distribution of outputs, 
and plots these alongside binned proxy data (see below).

%% model_functions.py
This script contains the following functions which, taken together, define and solve the forward model:

% forward_model: 
    Given parameter inputs, forward_model calculates the initial conditions for the carbon cycle e.g. equilibrium ocean chemistry and modern fluxes. 
    Proportionality constants for carbon cycle functions are also calculated from initial conditions. 
    Next, the ODE solver is called, and the system of equations describing the carbon cycle are solved. 
    The ODE solver only returns DIC and ALK as a function of time for both the ocean and the pore space. 
    These outputs are fed back into Cretaceous_cc to obtain the time evolution for carbon cycle fluxes and ocean chemistry. 
    Selected outputs are returned to Main_code_forward_modeling.py. 
    The function has additional lines of code that can be used to check mass balance, but these have been commented out for simplicity.

% system_of_equations:
    Contains the ODEs that describe the time evolution of the carbon cycle (equation 6 in manuscript). 
    The function takes the current state of the carbon cycle (calculated using the Cretaceous_cc function), 
    and returns the time derivatives of DIC and ALK in the ocean and the pore space. 
    This function is fed into the ODE solver to compute the time evolution of DIC and ALK for the ocean and pore space.

% weatherability_func: 
    Calculates weatherability as a function of time using equation 4 main text

% Cretaceous_cc: 
    This function computes equilibrium chemistry and carbon cycle fluxes from DIC and ALK, for both the ocean and the pore space. 
    Cretaceous_cc first calculates equilibrium chemistry for the ocean and the pore space using equations 20-25 in the main text, 
    then calculates surface and pore space temperatures from pCO2 using our climate model (equation 25) and surface-deep ocean relationship (equation 11). 
    Finally, carbon cycle fluxes are calculated from this information and prescribed [Ca] and shelf area through time. 
    Both carbon cycle fluxes and equilibrium chemistry are returned to system_of_equations or forward_model.

%% thermodynamic_variables.py:
    Contains the function Sol_prod that calculates carbonate solubility product 
    as a function of temperature from Pilson, M. E. (1998) "An Introduction to the Chemistry of the Sea", Prentice-Hall, Inc. 
    This function reproduces table G.1 for salinity=35.

%% mcmc_run_new.py:
    This script runs Markov Chain Monte Carlo inverse calculations as described in the paper. 
    This script defines the priors, calculates the likelihood function, and actually runs the MCMC analysis. 
    Then it calls some external plotting scripts "plotting_everthing" to plot everything. 
    To run MCMC calculations, python packages, "emcee" for the MCMC and "corner" for plotting are required. 
    Note that you may want to change the variable "threads" 
    This tells your computer how many cores to distribute the calculation over. 
    The variables nwalk and nsteps control how many forward model calls to do. 
    nwalk = 1000 and nsteps = 10000 are appropriate for a final calculation, but this could take a very long time. 
    For a quick test of functionality, nwalk = 50 and nsteps = 500 is adequate. 
    If you change the number of walkers and steps, be careful how you're defining the variable "production". 
    The default is to throw out the first 1000 model steps because it takes the walkers a while to converge on the posterior 
    (i.e. the first few hundred steps) are usually garbage. 
    But if you only have nsteps = 500 total then you won't be left with anything and will get errors!

%% plotting_everything.py:
    Called to plot MCMC outputs
    Fig

END EXPLANATION OF CODE STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
BINNED PROXY DATA EXPLANATION

The folder geochemical_proxies contains binned proxy data for atmospheric pCO2, ocean saturation state, surface and deep ocean temperature, 
and ocean pH. 
These text files are called by Main_code_forward_modeling.py which plots them alongside model outputs. T
Here is no text file for seafloor weathering since there is only a single mean Cretaceous data point that is defined in Main_code_forward_modeling.py script.

obs_CO2.txt: 
    binned atmospheric pCO2 proxies, plotted in Fig. S9 in the supplementary materials. 
    Rows refer to time (in Ma), pCO2 (in ppm), and uncertainty in pCO2 (in ppm). 
    See supplementary materials for references.

obs_omega_calc.txt: 
    binned calcite saturation state proxies, calculated from CCD proxies in Fig. S11 and equation S10 in supplementary materials. 
    Rows refer to time (in Ma), calcite saturation state, and uncertainty in calcite saturation state. 
    See supplementary materials for references.

obs_pH.txt: 
    binned ocean pH proxies, plotted in Fig. S8 in the supplementary materials. 
    Rows refer to time (in Ma), ocean pH, and uncertainty in ocean pH. 
    See supplementary materials for references.

obs_T.txt: 
    binned mean surface temperature proxies, plotted in Fig. S10 in the supplementary materials. 
    Rows refer to time (in Ma), mean surface temperature (in deg C), and uncertainty in mean surface temperature (deg C). 
    See supplementary materials for references.

obs_Td.txt: 
    binned deep ocean temperature proxies, plotted in Fig. S10 in the supplementary materials. 
    Rows refer to time (in Ma), mean deep ocean temperature (in deg C), and uncertainty in deep ocean temperature (deg C). 
    See supplementary materials for references.


END EXPLANATION OF DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

-------------------------
Contact e-mail: joshkt@uw.edu

