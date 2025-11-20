#! /usr/local/bin/python
import numpy
import pylab
import scipy.stats
#import function to calculate temperature dependent solubility product
from thermodynamic_variables import Sol_prod

# Given parameter inputs, forward_model calculates the initial conditions for the carbon cycle e.g. equilibrium ocean chemistry and modern fluxes.
# Proportionality constants for carbon cycle functions are also calculated from initial conditions.
# Next, the ODE solver is called, and the system of equations describing the carbon cycle are solved.
# The time evolution for selected outputs are returned to Main_code.py
def forward_model(W,F_outgass,n,CO2_dep,Te,mod_sea,alt_frac,Mp_frac,W_plus_1,cl_sens,change_out,F_carbw,frac_pel,CWF,deep_grad,coef_for_diss,beta,Ebas,PG):
    #import pdb
    #pdb.set_trace()
    # Make carbon chemistry constants and other parameters global variables:    
    global pK1,pK2,H_CO2,Pore_mod,k_d,k_w,k_c,k_o,Mp,ppCO2_o,Ea,Mo,ko1,ko2
    
    # Constants for carbon chemisty
    T=18+273 #Temperature held constant for the purposes of calculating equilibrium constants (see manuscript)
    pK1=17.788 - .073104 *T - .0051087*35 + 1.1463*10**-4*T**2
    pK2=20.919 - .064209 *T - .011887*35 + 8.7313*10**-5*T**2
    H_CO2=pylab.exp(9345.17/T - 167.8108 + 23.3585 * pylab.log(T)+(.023517 - 2.3656*10**-4*T+4.7036*10**-7*T**2)*35)
   
    Mo=1.35e21 #Mass of ocean (kg)
    W=Mo/W #convert timescale of circulation (yr) to water flux through pore space (kg/yr)
    
    Mp=Mp_frac*Mo #Calculate mass of ocean in pore space        
    partition=mod_sea/F_outgass #seafloor precipitation as a fraction of total carbon sink (used to initialize)

    ############################################################################
    ## Compute initial conditions (i.e. fluxes and ocean chemistry for modern carbon cycle):

    F_diss0=alt_frac*partition*F_outgass # Modern seafloor dissolution flux in pore space
    F_silic0=(1-partition)*F_outgass+(1-alt_frac)*partition*F_outgass #Modern continental silicaet weathering
    Precip_pore0=partition*F_outgass # Modern carbonate precipitation in pore-space
    
    ### Initial conditions for modern atmosphere and ocean
    pH_o=8.2 #pH of modern ocean
    ppCO2_o=.000280 #preindustrial pCO2 (ppm)
    CO2_aq_o=H_CO2 * ppCO2_o #modern ocean aqueous pCO2 molality
    hco3_o=(10**-pK1) * CO2_aq_o / (10**-pH_o) #modern ocean bicarbonate molality
    co3_o=(10**-pK2) * hco3_o / (10**-pH_o) # modern ocean carbonate molality 
    DIC_o=co3_o+hco3_o+CO2_aq_o #modern ocean dissolved inorganic carbon
    ALK_o=2*co3_o+hco3_o # modern ocean alkalinity
    Ca_o=0.0100278 # modern oceanic calcium molality (mol/kg)
    salt = ALK_o-2*Ca_o ## set salt to ensure ALK is consistent with Ca
    
    T_surface_diff=0 
    T_surface=T_surface_diff+285 # Modern average surface temperature (K)
    interc=274.037-deep_grad*285 # Calculate intercept for surface vs. deep ocean temperature relationship
    buffer_T=deep_grad*T_surface+interc # Calculate modern deep ocean temperature (K)
    Pore_mod=9.0 # Constant difference between deep ocean temperature and pore space temperature (K)
    
    ## Use these constants and steady state assumption to calculate remaining initial conditions:
    b1=2*F_diss0+W*ALK_o-2*W*DIC_o
    b2=2*F_silic0-W*ALK_o-2*F_outgass+2*W*DIC_o
    a1=W
    a2=-2*W
    a3=-W
    a4=2*W
    
    # Solving system of equations for pore space properties
    DIC_p=DIC_o-Precip_pore0/W #Modern dissolved inorganic carbon in pore space
    ALK_p=(b1-a2*DIC_p)/a1 #Modern alkalinity in pore space
    Precip_ocean0=(F_outgass+F_carbw)-(DIC_o-DIC_p)*W #Carbonate precipitation in ocean
    
    if DIC_p<0:
        print ("WARNING: NO PHYSIAL STEADY STATE EXISTS! DIC is negative")

    #Solve quadratic ignoring H in ALK but including CO2aq
    [vv,xx]=pylab.roots([ALK_p/(10**-pK1*10**-pK2),(ALK_p-DIC_p)/(10**-pK2),ALK_p-2*DIC_p])
    H_p=numpy.max([vv,xx]) # Modern hydrogen molality in pore space
    pH_p=-pylab.log10(H_p)  # Modern pH of pore space
    
    #Find remaining carbon chemistry from equilibrium conditions:
    CO2aq_p=DIC_p/((10**-pK1)*(10**-pK2)/H_p**2+10**-pK1/H_p+1) # Modern aqueous CO2 in pore space
    co3_p=ALK_p-DIC_p-H_p+CO2aq_p # Modern carbonate alkalinity in pore space
    hco3_p=co3_p*10**-pH_p/(10**-pK2) #Modern bicarbonate alkalinity in pore space (not used)
    
    Ca_p = 0.5*(ALK_p - salt) #Modern calcium molality in pore space
    omega_p=Ca_p *co3_p/Sol_prod(buffer_T+Pore_mod) #Modern saturation state pore space
    omega_o=Ca_o *co3_o/Sol_prod(T_surface) #Modern saturatino state of ocean
    
    ## All initital conditions have been calculated
    ############################################################################
    
    ############ Calculate proportionality constants ###########################
    # Given modern precipitation fluxes and saturation states, calculate proportionality constants
    k_c=Precip_pore0/(omega_p - 1)**n 
    k_o=Precip_ocean0/(omega_o - 1)**n

    ## Partition ocean carbonate sink into shelf precipitation and carbonate precipitation
    ko1=(1-frac_pel)*Precip_ocean0/(omega_o - 1)**n #Modern shelf precipitation
    ko2=frac_pel*Precip_ocean0/(omega_o**2.84) #Modern pelagic precipitation

    # Calculate remaining proportionality constants
    k_d=F_diss0/(2.88*10**-14*10**(-coef_for_diss*pH_p)*numpy.exp(-Ebas/(8.314*(buffer_T+Pore_mod)))) 
    k_w=F_silic0
    ##### all proportionality constants have been calculated ###################
              
    time=numpy.linspace(0,1e8,100) # Time array for outputs, 100 timesteps between 0 and 100 Ma
    #print ('IC new',DIC_o+1.8e20/Mo*ppCO2_o,ALK_o,DIC_p,ALK_p)    
    
    # Given initial conditions, unknown parameters, and system of ODEs (below), solve for time evolution of DIC and ALK of ocean and pore space:
    [out,mes]=scipy.integrate.odeint(system_of_equations, [DIC_o+1.8e20/Mo*ppCO2_o,ALK_o,DIC_p,ALK_p], time, args=(W,F_outgass,n,CO2_dep,Te,mod_sea,alt_frac,Mp_frac,W_plus_1,cl_sens,change_out,F_carbw,CWF,deep_grad,coef_for_diss,beta,Ebas,PG),full_output=1)
 
    # Create empty arrays for storing model outputs
    pH_array_o=0*time # pH of ocean
    CO2_array_o=0*time # pCO2 of ocean
    pH_array_p=0*time # pH of pore space
    CO2_array_p=0*time # pCO2 pore space
    Ca_array_o=0*time # Ca molality ocean
    Ca_array_p=0*time # Ca molality pore space
    CO3_array_o=0*time # Carbonate molality ocean
    CO3_array_p=0*time # Carbonate molality pore space
    HCO3_array_o=0*time # Bicarbonate molaity ocean
    HCO3_array_p=0*time # Bicarbonate molality pore space
    omega_o=0*time # Saturation state ocean
    omega_p=0*time # Saturation state pore space
    Tsurf_array=0*time # Surface temperature
    Tdeep_array=0*time # Pore space temperature
    Fd_array=0*time # Seafloor dissolution flux
    Fs_array=0*time # Continental silicate weathering flux
    Precip_ocean_array=0*time # Carbonate precipitaion ocean flux
    Precip_pore_array=0*time # Carbonate precipitation pore space flux
    volc_array=0*time # Volcanic outgassing flux 
    carbw_ar=0*time # Carbonate weathering flux
    co2_aq_o=0*time # aqueous CO2 ocean molality
    co2_aq_p=0*time # aqueous CO2 pore space molality
    
    ### Step through time, and fill output arrays with solution calculated above
    for i in range(0,len(time)):
        
        # Given time evolution of DIC and ALK for ocean and pore space, calculate time evolution for other carbon cycle variables (fluxes and equilibrium chemistry):
        [Ca,EM_H_o,EM_pH_o,EM_co3_o,EM_hco3_o,EM_co2aq_o,EM_ppCO2_o,EM_Ca_o,EM_H_p,EM_pH_p,EM_co3_p,EM_hco3_p,EM_co2aq_p,EM_ppCO2_p,EM_Ca_p,area,T_surface,T_surface_diff,buffer_T,EM_omega_o,EM_omega_p,Relative_outgassing,weatherability,F_outg,carb_weath,Precip_ocean,Precip_pore,F_diss,F_silic] = \
            Cretaceous_cc([out[i,0],out[i,1],out[i,2],out[i,3]],time[i],cl_sens,deep_grad,PG,change_out,F_outgass,W_plus_1,F_carbw,CWF,CO2_dep,Te,n,coef_for_diss,beta,Ebas)
        
        array_of_outputs=[EM_pH_o,EM_co3_o,EM_hco3_o, EM_co2aq_o,EM_ppCO2_o,EM_Ca_o,EM_omega_o,T_surface,buffer_T,F_diss,F_silic,Precip_ocean,Precip_pore,EM_omega_p,carb_weath,EM_co3_p,EM_hco3_p, EM_co2aq_p,EM_pH_p,EM_Ca_p, EM_co2aq_p]
        
        # Fill output array with appropriate model outputs:
        co2_aq_o[i]=array_of_outputs[3]
        pH_array_o[i]=array_of_outputs[0]
        CO2_array_o[i]=array_of_outputs[4] 
        co2_aq_p[i]=array_of_outputs[17]        
        pH_array_p[i]=array_of_outputs[18] 
        CO2_array_p[i]=array_of_outputs[20] 
        Ca_array_o[i]=array_of_outputs[5]
        Ca_array_p[i]=array_of_outputs[19]
        CO3_array_o[i]=array_of_outputs[1]
        CO3_array_p[i]=array_of_outputs[15]
        HCO3_array_o[i]=array_of_outputs[2]
        HCO3_array_p[i]=array_of_outputs[16]
        omega_o[i]=array_of_outputs[6]
        omega_p[i]=array_of_outputs[13]
        carbw_ar[i]=array_of_outputs[14]
        Tsurf_array[i]=array_of_outputs[7] 
        Tdeep_array[i]=array_of_outputs[8]
        Fd_array[i]=array_of_outputs[9]  
        Fs_array[i]=array_of_outputs[10]
        Precip_ocean_array[i]=array_of_outputs[11]
        Precip_pore_array[i]=array_of_outputs[12]
        volc_array[i]=F_outgass+change_out*F_outgass*(time[i]/1e8)
    
    ############################################################################
    #### Checking for mass balance
    ############################################################################
    #max_el=numpy.size(time)-1
    #tstep=numpy.max(time)/max_el #duration of timestep
    #sources=numpy.sum(volc_array[2:]+carbw_ar[2:])*tstep
    #sinks=numpy.sum(Precip_ocean_array[2:]+Precip_pore_array[2:])*tstep
    #flux_bal=sources-sinks
    #res_change=1.8e20*(CO2_array_o[max_el]-CO2_array_o[2])+Mo*(HCO3_array_o[max_el]+CO3_array_o[max_el]+co2_aq_o[max_el]-co2_aq_o[2]-CO3_array_o[2]-HCO3_array_o[2])+Mp*(HCO3_array_p[max_el]+CO3_array_p[max_el]+co2_aq_p[max_el]-co2_aq_p[2]-CO3_array_p[2]-HCO3_array_p[2])
    #res_change2=(out[max_el,0]-out[2,0])*Mo+(out[max_el,2]-out[2,2])*Mp
    #print flux_bal,res_change,res_change2,Mo*(HCO3_array_o[max_el]+CO3_array_o[max_el]+co2_aq_o[max_el]-co2_aq_o[2]-CO3_array_o[2]-HCO3_array_o[2])+Mp*(HCO3_array_p[max_el]+CO3_array_p[max_el]+co2_aq_p[max_el]-co2_aq_p[2]-CO3_array_p[2]-HCO3_array_p[2])
    #imbalance=(flux_bal-res_change)/(out[max_el,0]*Mo)
    #print imbalance
    imbalance = 0 ## Comment out if wish to calculate mass imbalance
    ############################################################################
    
    # Return selected outputs to Main_code.py for plotting
    # out[:,0] = DIC_o+1.8e20/Mo*ppCO2_o; out[:,1] = ALK_o; out[:,2] = DIC_p; out[:,3] = ALK_p
    return [out[:,0],out[:,1],out[:,2],out[:,3],time,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array],imbalance
   
    
# system_of_equations contains the ODEs that describe the time evolution of the carbon cycle (equation 6 in manuscript). 
# The function takes the current state of the carbon cycle, and returns the time derivatives of DIC and ALK in the ocean and the pore space
def system_of_equations (y,t0,W,F_outgass,n,CO2_dep,Te,mod_sea,alt_frac,Mp_frac,W_plus_1,cl_sens,change_out,F_carbw,CWF,deep_grad,coef_for_diss,beta,Ebas,PG):
    
    global pK1,pK2,H_CO2,Pore_mod,k_d,k_w,k_c,k_o,Mp,ppCO2_o,Ea,Mo,ko1,ko2
    
    s=1.8e20/Mo #correction factor for mass balance (see manuscript) 
    
    # Call Cretaceous carbon cycle forward model to obtain current state of carbon cycle, namely ocean chemistry and carbon cycle fluxes, 
    # from parameters, time, and DIC and ALK for ocean and pore space.
    [Ca,EM_H_o,EM_pH_o,EM_co3_o,EM_hco3_o,EM_co2aq_o,EM_ppCO2_o,EM_Ca_o,EM_H_p,EM_pH_p,EM_co3_p,EM_hco3_p,EM_co2aq_p,EM_ppCO2_p,EM_Ca_p,area,T_surface,T_surface_diff,buffer_T,EM_omega_o,EM_omega_p,Relative_outgassing,weatherability,F_outg,carb_weath,Precip_ocean,Precip_pore,F_diss,F_silic] = \
        Cretaceous_cc([y[0],y[1],y[2],y[3]],t0,cl_sens,deep_grad,PG,change_out,F_outgass,W_plus_1,F_carbw,CWF,CO2_dep,Te,n,coef_for_diss,beta,Ebas)
       
    DICo=y[0]-EM_ppCO2_o*s # Correct dissolved inorganic carbon content of ocean by subtracting atmospheric pCO2 from atmosphere-ocean reservoir (this ensures mass balance)
   
    # Calculate time derivatives from current state (equation 6).
    dy0_dt=F_outg/Mo+carb_weath/Mo-W*(DICo-y[2])/Mo-Precip_ocean/Mo # Time derivative of atmosphere-ocean carbon abundance
    dy1_dt=2*carb_weath/Mo-W*(y[1]-y[3])/Mo+2*F_silic/Mo-2*Precip_ocean/Mo #Time derivative of ocean alkalinity
    dy2_dt=W*(DICo-y[2])/Mp-Precip_pore/Mp # Time derivative of pore space carbon abundance
    dy3_dt=W*(y[1]-y[3])/Mp+2*F_diss/Mp-2*Precip_pore/Mp #Time derivative of pore space alkalinity
    
    return [dy0_dt,dy1_dt,dy2_dt,dy3_dt]                   

# This function calculates weatherability as a function of time
def weatherability_func (t0,W_plus_1):
    weather_coef=1-(1-W_plus_1)*abs(t0)/1e8 ## effectively 1+W*t0/1e8 - see equation 4 main text
    return weather_coef
   
## This function computes equilibrium chemistry and carbon cycle fluxes from DIC and ALK, for both the ocean and the pore space.  
## Cretaceous_cc first calculates equilibrium chemistry for the ocean and the pore space using equations 20-25 in the main text, 
## then calculates surface and pore space temperatures from pCO2 using our climate model (equation 25) and surface-deep ocean relationship (equation 11).
## Finally, carbon cycle fluxes are calculated from this information and prescribed [Ca] and shelf area through time.
## Both carbon cycle fluxes and equilibrium chemistry are returned to system_of_equations or forward_model.
def Cretaceous_cc (y,t0,cl_sens,deep_grad,PG,change_out,F_outgass,W_plus_1,F_carbw,CWF,CO2_dep,Te,n,coef_for_diss,beta,Ebas): #inputs include DIC and ALK for ocean and pore space, time, and key parameters
    
    global pK1,pK2,H_CO2,Pore_mod,k_d,k_w,k_c,k_o,Mp,ppCO2_o,Ea,Mo,ko1,ko2
    
    # Perscribed evolution of ocean and pore-space calcium molality (Fig. S6): 
    Ca=7.00658e-27*abs(t0)**3-1.9847e-18*(t0)**2+2.4016e-10*(t0)+0.0100278
       
    #ocean variables
    s=1.8e20/Mo #Constant for ensuring mass conversation
    [c,d]=pylab.roots([y[1]/(10**-pK2*10**-pK1)*(1+s/H_CO2),(y[1]-y[0])/(10**-pK2),(y[1]-2*y[0])])
    EM_H_o=numpy.max([c ,d]) ## Evolving model ocean hydrogen molality
    EM_pH_o=-pylab.log10(EM_H_o) ## Evolving model ocean pH
    EM_co3_o=y[1]/(2+EM_H_o/(10**-pK2)) ## Evolving model ocean carbonate molality
    EM_hco3_o=y[1]-2*EM_co3_o ## Evolving model ocean bicarbonate molality
    EM_co2aq_o=( EM_hco3_o*EM_H_o/(10**-pK1) ) ## Evolving model ocean aqueous CO2 molality
    EM_ppCO2_o = EM_co2aq_o /H_CO2 ## Evolving model atmospheric pCO2
    EM_Ca_o = Ca #Perscribed ocean Ca molality
    
    #pore space variables
    [c,d]=pylab.roots([y[3]/(10**-pK2*10**-pK1),(y[3]-y[2])/(10**-pK2),(y[3]-2*y[2])])
    EM_H_p=numpy.max([c ,d]) ## Evolving model pore space hydrogen molality
    EM_pH_p=-pylab.log10(EM_H_p) ## Evolving model pore space pH
    EM_co3_p=y[3]/(2+EM_H_p/(10**-pK2))  ## Evolving model pore space carbonate molality
    EM_hco3_p=y[3]-2*EM_co3_p ## Evolving model pore space bicarbonate molality
    EM_co2aq_p=( EM_hco3_p*EM_H_p/(10**-pK1) ) ## Evolving model pore space aqueous CO2 molality
    EM_ppCO2_p = EM_co2aq_p /H_CO2 ## Evolving model pore space pCO2 (unused)
    EM_Ca_p = Ca #Percribed pore space Ca molality

    # Perscribed evolution of continental shelf area (Fig. S5):
    area=(.540e-23*abs(t0)**3-3.072e-15*abs(t0)**2+2.9997e-07*abs(t0)+16.089)/16.089    
 
    # Apply climate model to calculate surface temperature (equation 5):
    T_surface_diff=(cl_sens/numpy.log(2))*numpy.log(EM_ppCO2_o/ppCO2_o)-cl_sens*abs(t0)/(0.181e9*1.26) + PG*abs(t0)/1e8
    T_surface=T_surface_diff+285
    # Calculate pore space temperature (equation 11):
    interc=274.037-deep_grad*285
    buffer_T=deep_grad*T_surface+interc

    F_outg=F_outgass+change_out*F_outgass*(abs(t0)/1e8) #Outgassing as a function of time, equation 9
    Relative_outgassing= F_outg/F_outgass #Relative outgassing
    weatherability=weatherability_func(abs(t0),W_plus_1) #Weatherability as a functino of time
    
    #Calculate saturation state of ocean and pore space (equation 19):
    EM_omega_o=EM_Ca_o *EM_co3_o/(Sol_prod(T_surface))
    EM_omega_p=EM_Ca_p *EM_co3_p/(Sol_prod(buffer_T+Pore_mod))    

    # Carbonate weathering flux (equartion 7)
    carb_weath=F_carbw*(1+CWF*abs(t0)/1e8)*(EM_ppCO2_o/ppCO2_o)**(CO2_dep)*numpy.exp((T_surface_diff)/Te)*weatherability   
    # Ocean carbonate precipitation flux (equation 17)
    Precip_ocean=ko1*area*(EM_omega_o - 1)**n+ko2*(EM_omega_o**2.84)  #shelf + pelagic                   
    # Pore space carbonate precipitation flux (equation 18)
    Precip_pore = k_c*(EM_omega_p -1)**n
    #seafloor dissolution flux (equation 12)
    F_diss = k_d*2.88*10**-14*10**(-coef_for_diss*EM_pH_p)*Relative_outgassing**beta*numpy.exp(-Ebas/(8.314*(buffer_T+Pore_mod)))
    #continental silicate weathering flux (equation 2)
    F_silic=k_w*weatherability*(EM_ppCO2_o/ppCO2_o)**CO2_dep*numpy.exp((T_surface_diff)/Te)

    return [Ca,EM_H_o,EM_pH_o,EM_co3_o,EM_hco3_o,EM_co2aq_o,EM_ppCO2_o,EM_Ca_o,EM_H_p,EM_pH_p,EM_co3_p,EM_hco3_p,EM_co2aq_p,EM_ppCO2_p,EM_Ca_p,area,T_surface,T_surface_diff,buffer_T,EM_omega_o,EM_omega_p,Relative_outgassing,weatherability,F_outg,carb_weath,Precip_ocean,Precip_pore,F_diss,F_silic]
