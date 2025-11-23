# This set of python scripts runs a geological carbon cycle model over the last 100 Ma and plots the outputs alongside proxy data from the literature.
# The model is described in J. Krissansen-Totton and D.C. Catling (2017) "Constraining climate sensitivity and continental versus seafloor weathering 
# using an inverse geological carbon cycle model", Nature Communications, in revision. 

# As a matter of courtesy, we request that people using this code please cite Krissansen-Totton et al. (2017). In the interest of an "open source" approach,
# we also request that authors who use and modify the code, please send a copy of papers and modified code to the lead author  (joshkt@uw.edu)

#! /usr/local/bin/python
import numpy
import pylab
import scipy.stats

################################################################################
################################################################################
# Define uniform range for 17 uncertain variables.
# Compare to Table 1, column 2 in the main text
# Set lower limits and upper limits to the same value for point estimates
################################################################################
# CO2 dependence of continental weathering
[CO2_dep_low,CO2_dep_high]=[0.2,0.5]

# e-folding temperature continental weathering (K)
[Te_low,Te_high]=[5.,15.] # [30.,40.] to reproduce Fig. 3 and 4.

#Relative Cretaceous weatherability
[W_plus_1_low,W_plus_1_high]=[1.0,1.0] #[0.4,0.6] to reproduce Fig. 4

#Climate sensitivity (K per CO2 doubling)
[cl_sens_low,cl_sens_high]=[1.5,8.0] 

# Change in outgassing, V
[out_low,out_high]=[0.2,1.5]

#Carbonate weatherability factor
[C_WF_low,C_WF_high]=[-0.9,1.5]

#Modern outgassing flux (mol C/yr)
[F_outgass_low,F_outgass_high]=[4e12,10e12
                                ]
#Modern carbonate weathering flux (mol/yr)
[F_carbw_low,F_carbw_high]=[7e12,14e12]

#Pore space circulation time (yr)
[Tau_low,Tau_high]=[2e4,1e6]

#Carbonate precipiataion exponent
[n_low,n_high]=[1.0,2.5]

#Modern seafloor dissolution relative to precipitation
[alt_frac_low,alt_frac_high]=[0.5,1.5]

#Surface-deep temperature gradient, a_grad
[deep_grad_low,deep_grad_high]=[0.8,1.4]

#pH dependence seafloor dissolution, gamma 
[coef_for_diss_low,coef_for_diss_high]=[0,0.5]

#Temperature dependence seafloor weathering (J/mol)
[Ebas_low,Ebas_high]=[40000,110000]

#Modern pelagic fraction
[frac_pel_low,frac_pel_high]=[0.4,0.6]

#Spreading rate dependence
[beta_low,beta_high]=[0.0,1.0]

#Paleogeography climate parameter (K)
[PG_low,PG_high]=[0,5.0]

iterations=100 # Number of forward model calls, typically 1000-10000 to build up reliable confidence intervals
################################################################################
################################################################################

from model_functions import forward_model #import scripts that execute the forward model
#imbalance_array=[] # for checking mass balance

#Fixed parameters
Mp_frac=0.01 #pore space volume as fraction of whole ocean  
mod_sea=.45e12 # modern seafloor dissolution flux (mol/yr)
ppCO2=10**-6  #CO2 conversion factor from ppm to bar

# Run forward model once with arbitrary inputs and use outputs to define size of output array
[fmod1,imbal]=forward_model(6.0e4,6e12,1.7,0.3,13.7,.45e12,0.7,0.01,0.7,1.0,1.0,11e12,0.5,0.0,1.0,0.0,0.0,92000,0.0)
all_output=numpy.zeros([len(fmod1),len(fmod1[4]),iterations]) #has dimensions outputs X timesteps X iterations
all_output[:,:,0]=fmod1

ii=0
while ii<iterations:
    print (ii," of ",iterations)

    ### Sample uniform distribution for unknown parameters given ranges defined above
    CO2_dep=numpy.random.uniform(CO2_dep_low,CO2_dep_high) # CO2-depenence continental weathering, alpha
    Te=numpy.random.uniform(Te_low,Te_high) # e-folding temperature continental weathering
    W_plus_1=numpy.random.uniform(W_plus_1_low,W_plus_1_high) #Relative Cretaceous weatherability
    cl_sens=numpy.random.uniform(cl_sens_low,cl_sens_high) #Climate sensitivity
    change_out=numpy.random.uniform(out_low,out_high) # Change in outgassing, V 
    C_WF=numpy.random.uniform(C_WF_low,C_WF_high) #Carbonate weatherability factor
    F_outgass=numpy.random.uniform(F_outgass_low,F_outgass_high) #Modern outgassing flux (mol C/yr)
    F_carbw=numpy.random.uniform(F_carbw_low,F_carbw_high) #Modern carbonate weathering flux (mol/yr)
    Tau=numpy.random.uniform(Tau_low,Tau_high) #Pore space circulation time
    n=numpy.random.uniform(n_low,n_high) #Carbonate precipiataion exponent
    alt_frac=numpy.random.uniform(alt_frac_low,alt_frac_high) #Modern seafloor dissolution relative to precipitation
    deep_grad=numpy.random.uniform(deep_grad_low,deep_grad_high) #Surface-deep temperature gradient, a_grad
    coef_for_diss=numpy.random.uniform(coef_for_diss_low,coef_for_diss_high) #pH dependence seafloor dissolution, gamma 
    Ebas=numpy.random.uniform(Ebas_low,Ebas_high) #Temperature dependence seafloor weathering (J/mol)
    frac_pel=numpy.random.uniform(frac_pel_low,frac_pel_high) #Modern pelagic fraction
    beta=numpy.random.uniform(beta_low,beta_high) #Spreading rate dependence
    PG=numpy.random.uniform(PG_low,PG_high) #Paleogeography climate parameter (K)
    ############################################################################
    
    ## Attempt to run forward model
    try:
        [all_output[:,:,ii],imbal]=forward_model(Tau,F_outgass,n,CO2_dep,Te,mod_sea,alt_frac,Mp_frac,W_plus_1,cl_sens,change_out,F_carbw,frac_pel,C_WF,deep_grad,coef_for_diss,beta,Ebas,PG)   
        # Check if model parameters produce valid outputs:
        if (numpy.isnan(all_output[7,98,ii]))or(all_output[14,98,ii]<0.0):
            print ("error") #repeat forward model iteration
        else: # Model produce valid outputs:
            ii=ii+1
            #### for checking mass balance
            #imbalance_array.append(imbal) 
            #if abs(imbal)>0.5: 
            #    print "LARGE IMBALANCE"
            #    print Tau,F_outgass,n,CO2_dep,tdep_weath,mod_sea,alt_frac,Mp_frac,W_plus_1,cl_sens,change_out,F_carbw,frac_pel,C_WF,deep_grad,coef_for_diss,beta
    
    # Rerun this iteration if combination of parameters produces errors 
    except:
        print ("error, try again")

## Define 90% confidence intervals for key outputs:
# pH ocean and pore space      
confidence_pH_o=scipy.stats.scoreatpercentile(all_output[5,:,:],[5,50,95], interpolation_method='fraction',axis=1)
confidence_pH_p=scipy.stats.scoreatpercentile(all_output[7,:,:],[5,50,95], interpolation_method='fraction',axis=1)
# atmospheric pCO2
confidence_CO2o=scipy.stats.scoreatpercentile(all_output[6,:,:],[5,50,95], interpolation_method='fraction',axis=1)  
# Saturation state ocean and pore space
confidence_omega_o=scipy.stats.scoreatpercentile(all_output[15,:,:],[5,50,95], interpolation_method='fraction',axis=1)
confidence_omega_p=scipy.stats.scoreatpercentile(all_output[16,:,:],[5,50,95], interpolation_method='fraction',axis=1)
# Surface and deep ocean temperature
confidence_Tsurf=scipy.stats.scoreatpercentile(all_output[17,:,:],[5,50,95], interpolation_method='fraction',axis=1)
confidence_Tdeep=scipy.stats.scoreatpercentile(all_output[18,:,:],[5,50,95], interpolation_method='fraction',axis=1)
# Seafloor dissolution and continental weathering fluxes
confidence_Fd=scipy.stats.scoreatpercentile(all_output[19,:,:],[5,50,95], interpolation_method='fraction',axis=1)
confidence_Fs=scipy.stats.scoreatpercentile(all_output[20,:,:],[5,50,95], interpolation_method='fraction',axis=1)
# Carbonate precipitation ocean and pore space
confidence_Prec_o=scipy.stats.scoreatpercentile(all_output[21,:,:],[5,50,95], interpolation_method='fraction',axis=1)
confidence_Prec_p=scipy.stats.scoreatpercentile(all_output[22,:,:],[5,50,95], interpolation_method='fraction',axis=1)

# Create distribution for 'observed' Cretaceous dissolution flux from observed carbon content
# This is calculated from outgassing change (V), and spreading rate dependence on outgassing change (beta)
spread_dist=[] # create array to store distribution of dissolution rate changes
for kk in range(0,10000): # sample a large number to obtain uncertainty in observed dissolution change
    rate=1+numpy.random.uniform(out_low,out_high) #sample 1+V
    beta_plot=numpy.random.uniform(beta_low,beta_high) #sample beta
    precip_e=(numpy.random.uniform(-.75,0.75)+2.35)*10**12 #Observed change in oceanic crust carbon content
    spread_dist.append(rate**beta_plot*precip_e) #calculate change in seafloor dissolution and add to array
spread_dist=numpy.array(spread_dist) #converty to numpy array
[low,med,high]=scipy.stats.scoreatpercentile(spread_dist,[5,50,95], interpolation_method='fraction',axis=0) #90% confidence interval

########################################################
########################################################
### Plot figure for publication
########################################################
all_output[4,:,0]=all_output[4,:,0]/1e6 #convert times to millions of years for plotting
strt_lim=-.01e2 #start time for plotting (Ma)
fin_lim=1.01e2 #end time for plotting (Ma)

# Create 2x3 subplot to recreate Fig. 2-4 in  main text
pylab.figure(figsize=(24,12))
pylab.subplot(2, 3, 1)  # First subplot shows ocean pH through time
pylab.plot(all_output[4,:,0],confidence_pH_o[1],'k',label='ocean') #median model output
pylab.fill_between(all_output[4,:,0], confidence_pH_o[0], confidence_pH_o[2], color='grey', alpha=0.4) #90% confidence model output
observ_pH_saved=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_pH.txt',delimiter=',') #load proxy pH data from text file
pylab.errorbar(observ_pH_saved[0,:],observ_pH_saved[1,:],observ_pH_saved[2,:],color='k',marker='o',linestyle="None") #plot proxy data with uncertainty
## Define axes, labels, and limits of x-axis:
pylab.xlabel('Time (Ma)')
pylab.ylabel('Ocean pH')
pylab.xlim([strt_lim,fin_lim])
pylab.text(-15, 8.35, 'A', fontsize=16, fontweight='bold', va='top')

pylab.subplot(2, 3, 2) # Second subplot shows atmospheric pCO2 through time
pylab.plot(all_output[4,:,0],confidence_CO2o[1]/ppCO2,'k',label='RCO2') # median model output
pylab.fill_between(all_output[4,:,0], confidence_CO2o[0]/ppCO2, confidence_CO2o[2]/ppCO2, color='grey', alpha=0.4) #90% confidence model output
observ_CO2=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_CO2.txt',delimiter=',') #load proxy CO2 data from text file
pylab.errorbar(observ_CO2[0,:],observ_CO2[1,:],observ_CO2[2,:],color='k',marker='o',linestyle="None") #plot proxy data with uncertainty
## Define axes, labels, and limits of x-axis:
pylab.xlabel('Time (Ma)')
pylab.ylabel('Atmospheric CO$_2$ (ppm)')
pylab.xlim([strt_lim,fin_lim])
pylab.text(-15, 2400, 'B', fontsize=16, fontweight='bold', va='top')

pylab.subplot(2,3,3) # Third subplot shows ocean saturation state through time
pylab.plot(all_output[4,:,0],confidence_omega_o[1],'k',label='ocean') # median model output
pylab.fill_between(all_output[4,:,0], confidence_omega_o[0], confidence_omega_o[2], color='grey', alpha=0.4) #90% confidence model output
observ_omega=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_omega_calc.txt',delimiter=',') #load proxy CO2 data from text file
pylab.errorbar(observ_omega[0,:],observ_omega[1,:],observ_omega[2,:],color='k',marker='o',linestyle="None") #plot proxy data with uncertainty
## Define axes, labels, and limits of x-axis:
pylab.ylabel('Saturation state')
pylab.xlabel('Time (Ma)')
pylab.xlim([strt_lim,fin_lim])
pylab.text(-15, 4.9, 'C', fontsize=16, fontweight='bold', va='top')    

pylab.subplot(2, 3, 4) # Fourth subplot shows surface and deep ocean temperature through time
pylab.plot(all_output[4,:,0],confidence_Tsurf[1],'k',label='Surface') # median surface temp model output
pylab.fill_between(all_output[4,:,0], confidence_Tsurf[0], confidence_Tsurf[2], color='grey', alpha=0.4) #90% confidence model output
pylab.plot(all_output[4,:,0],confidence_Tdeep[1],'r',label='Deep ocean') # median deep ocean temp model output
pylab.fill_between(all_output[4,:,0], confidence_Tdeep[0], confidence_Tdeep[2], color='red', alpha=0.4) #90% confidence model output
observ_Td=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_Td.txt',delimiter=',') #load deep ocean temperature data from text file
pylab.errorbar(observ_Td[0,:],observ_Td[1,:]+273.15,observ_Td[2,:],color='r',marker='o',linestyle="None") #plot proxy data with uncertainty
observ_T=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_T.txt',delimiter=',')  #load surface temperature data from text file
pylab.errorbar(observ_T[0,:],observ_T[1,:]+273.15,observ_T[2,:],color='k',marker='o',linestyle="None") #plot proxy data with uncertainty
## Define axes, labels, and limits of x-axis:
pylab.ylabel('Temperature (K)')
pylab.xlabel('Time (Ma)')
pylab.legend(loc=2)
pylab.xlim([strt_lim,fin_lim])
pylab.text(-15, 305, 'D', fontsize=16, fontweight='bold', va='top')    

pylab.subplot(2, 3, 5) # Fifth subplot shows continental weathering flux and ocean carbonate precipitation through time
pylab.plot(all_output[4,:,0],confidence_Fs[1]/1e12,'r',label='Cont. weathering') # median cont. weathering model output, converted to Tmol/yr
pylab.fill_between(all_output[4,:,0], confidence_Fs[0]/1e12, confidence_Fs[2]/1e12, color='red', alpha=0.4)  #90% confidence model output
pylab.plot(all_output[4,:,0],confidence_Prec_o[1]/1e12,'k',label='ocean precip.') # median ocean precipitation model output, converted to Tmol/yr
pylab.fill_between(all_output[4,:,0], confidence_Prec_o[0]/1e12, confidence_Prec_o[2]/1e12, color='grey', alpha=0.4)  #90% confidence model output
## Define axes, labels, and limits of x-axis:
pylab.ylabel('Fluxes (Tmol/yr)')
pylab.xlabel('Time (Ma)')  
pylab.xlim([strt_lim,fin_lim])
pylab.legend(loc=2)
pylab.text(-14, 60, 'E', fontsize=16, fontweight='bold', va='top')        

pylab.subplot(2, 3, 6)  # Sixth subplot shows seafloor dissolution flux and seafloor carbonate precipitation through time
pylab.plot(all_output[4,:,0],confidence_Fd[1]/1e12,'r',label='Seafloor dissolution') # median seafloor dissolution model output, converted to Tmol/yr
pylab.fill_between(all_output[4,:,0], confidence_Fd[0]/1e12, confidence_Fd[2]/1e12, color='red', alpha=0.4) #90% confidence model output
pylab.plot(all_output[4,:,0],confidence_Prec_p[1]/1e12,'k',label='pore precip.') # median pore space carbonate precp. model output, converted to Tmol/yr
pylab.fill_between(all_output[4,:,0], confidence_Prec_p[0]/1e12, confidence_Prec_p[2]/1e12, color='grey', alpha=0.4) #90% confidence model output
## Define axes, labels, and limits of x-axis:
pylab.legend(loc=2)     
pylab.xlim([strt_lim,fin_lim])
pylab.ylabel('Fluxes (Tmol/yr)')
pylab.xlabel('Time (Ma)')  
pylab.errorbar(98.,med/1e12,yerr=numpy.array([[med-low],[high-med]])/1e12,color='k',marker='o',linestyle="None") #Plot observed Cretaceous seafloor precip. flux (Tmol/yr)
pylab.legend(loc=2)  
pylab.text(-12, 6.0, 'F', fontsize=16, fontweight='bold', va='top')    
pylab.show()
pylab.tight_layout()

## Plot distribution mass imbalance
#pylab.figure()
#pylab.hist(imbalance_array,500)
#pylab.show()

########################################################
## End figure plotting
########################################################
########################################################

