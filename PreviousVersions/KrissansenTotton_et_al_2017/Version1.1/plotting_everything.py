#! /usr/local/bin/python
import numpy
import pylab

import scipy.stats
from scipy.interpolate import interp1d

def mc_plotter(all_output,decision,legend_counter):
    
    F_outgass=6e12 #dodgy, used in error calculations
    
    
    #pH from Pearson Palmer
    tpH1=numpy.array([.085,.98,1.49,3.0,3.31,3.87,6,6.2,9.02,10.39,11.4,11.81,13.06,14.73,14.96,16.23,16.7,18.38,19.85,21.7,23.,23.51])*10**6#,34.84,36.10,39.51])*10**6
    pH1=numpy.array([8.1,8.12,8.13,8.21,8.17,8.14,8.15,8.12,8.20,8.18,8.19,8.16,8.20,8.31,8.26,8.14,8.18,8.20,8.20,8.19,8.12,8.04])
    tpH2=numpy.array([40.12,42.52,44.26,45.69,46.07,46.97,50.33,51.02,52.22,53.24,55.84,57.12,59.88])*10**6
    pH2=numpy.array([7.8,8.07,7.95,7.79,7.54,7.99,7.84,7.92,7.42,7.62,7.48,7.54,7.42])
    #pylab.figure()
    if decision=="n":
        pylab.figure(figsize=(30,15))
    pylab.subplot(3, 4, 1)
    pylab.plot(all_output[4][:],all_output[5][:],'r',label='ocean')
    pylab.plot(all_output[4][:],all_output[7][:],'b',label='pore space')
    if legend_counter==0:
        pylab.plot(tpH1,pH1,'ro',linestyle="-")
        pylab.plot(tpH2,pH2,'ro',linestyle="-")
        pylab.xlabel('Time (yr)')
        pylab.ylabel('pH')
        pylab.legend()
    
    # CO2 
    #Modern CO2 for reference
    ppCO2=10**-6#0.000420318058799 #model modern
    preinudsmod=1.0#280.0
    #Cretaceous CO2: from Hong Lee 2012
    tCO2=numpy.array([65,65.5,66,66.5,67,68,68,75,76.5,80,83,91,95,98,100.5,102,103.5,107.5,108,113.5,115,115,120,122,125,129,143])*10**6
    CO2v=numpy.array([406,782,495,437,340,171,456,1412,656,917,1522,1437,1626,1520,1368,1428,1060,1219,907,449,1117,1325,798,1024,701,309,788])/preinudsmod
    CO2er=numpy.array([5,95,83,96,91,126,201,310,180,218,173,366,700,228,68,128,76,431,424,140,97,333,157,153,511,78,114])/preinudsmod
    #Cenozoic CO2 from Beerling Royer
    CO2_temp=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Cenozoic_CO2_Beerling_Royer.txt',delimiter=',')
    #pylab.figure()
    pylab.subplot(3, 4, 2)
    pylab.plot(all_output[4][:],all_output[6][:]/ppCO2,'r',label='RCO2')
    if legend_counter==0:
        pylab.errorbar(tCO2,CO2v,yerr=CO2er,color='r',marker='o',linestyle="None")
        pylab.plot(CO2_temp[:,0],CO2_temp[:,1]/preinudsmod,color='r',marker='o',linestyle="None")
        pylab.xlabel('Time (yr)')
        pylab.ylabel('CO2 relative to modern')
        pylab.legend(loc=2)
    
        
    
      
    #pylab.figure()
    pylab.subplot(3, 4, 3)
    pylab.plot(all_output[4][:],all_output[9][:],'r',label='Ca ocean')
    pylab.plot(all_output[4][:],all_output[10][:],'b',label='Ca pore')
    #pylab.plot(all_output[4][:],Ca_fun(all_output[4][:]),'rx',label="Ca proxie Horita") #optional curve
    # Alternatively use data points from Horita table 2, and Cretaceous 94 Ma value from Timofeeff 2006
    tCa=numpy.array([5,14,35,37,94])*10**6
    Ca_prox_low=numpy.array([7,9,12,11,20])*10**-3
    Ca_prox_high=numpy.array([15,18,21,20,28])*10**-3
    Ca_prox=numpy.array([12,14,17,16,26])*10**-3
    if legend_counter==0:
        pylab.errorbar(tCa,Ca_prox,yerr=[Ca_prox-Ca_prox_low, Ca_prox_high-Ca_prox],color='r',marker='o',linestyle="None")
            
    pylab.plot(all_output[4][:],all_output[13][:],'k',label='HCO3 ocean')
    pylab.plot(all_output[4][:],all_output[14][:],'g',label='HCO3 pore')
    if legend_counter==0:
        pylab.xlabel('Molality (mol/kg)')
        pylab.xlabel('Time (yr)')
        pylab.legend(loc=2)
    
    #pylab.figure()
    pylab.subplot(3, 4, 4)
    # CO3 from Tyrrel and Zeebe (reconstructed from Ca and omega proxies)
    pylab.plot(all_output[4][:],all_output[11][:],'r',label='CO3 ocean')
    pylab.plot(all_output[4][:],all_output[12][:],'b',label='CO3 pore')
    CO3_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Tyrrell_Zeebe_CO3.txt',delimiter=',')
    CO3_f=interp1d(numpy.flipud(CO3_data[:,0]),numpy.flipud(CO3_data[:,1]))
    if legend_counter==0:
        pylab.plot(all_output[4][:],CO3_f(all_output[4][:]),'ro')
        #pylab.plot(CO3_data[:,0],CO3_data[:,1],'ko')
        pylab.xlabel('Time (yr)')
        pylab.ylabel('Molality (mol/kg)')
        pylab.legend()    
            
    #omega and CCD (from Tyrrel and Zeebe)
    #pylab.figure()
    pylab.subplot(3, 4, 5)
    pylab.plot(all_output[4][:],all_output[15][:],'r',label='ocean')
    pylab.plot(all_output[4][:],all_output[16][:],'b',label='pore space')
    #do proxy stuff
    CCD_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Tyrrell_Zeebe_CCD.txt',delimiter=',')
    CCD_f=interp1d(numpy.flipud(CCD_data[:,0]),numpy.flipud(CCD_data[:,1]))
    # convert CCD to omega using equation 4 in Jansen et al. 2002, find coefficient first
    #k_CCD=all_output[15][0]/numpy.exp(0.176*(CCD_f(0)-3.06))
    k_CCD=all_output[15][0]/numpy.exp(0.189*(CCD_f(0)-3.82))
    CCD=CCD_f(all_output[4][:])
    #omega_proxy=k_CCD*numpy.exp(0.176*(CCD-3.06))
    omega_proxy=k_CCD*numpy.exp(0.189*(CCD-3.82))
    if legend_counter==0:
        pylab.plot(all_output[4][:],omega_proxy,'ro')
        pylab.legend(loc=2)
        pylab.ylabel('Saturation state')
        pylab.xlabel('Time (yr)')
    
    # Temperature Cretaceous plus pre-PETM
    t1=numpy.array([55.,67.5,75.,90.,100.])*10**6
    T1=numpy.array([23.029,19.03,21.56,26.68,25.35])+273.15 #surface
    er_T1=numpy.array([2.0,5.35,3.99,4.06,3.53])
    t1b=numpy.array([67.5,75.,90.,100.])*10**6
    T2=numpy.array([11,10,19,16])+273.15 #shallow
    er_T2=numpy.array([1,1,1,4])
    #pylab.figure()
    pylab.subplot(3, 4, 6)
    pylab.plot(all_output[4][:],all_output[17][:],'r',label='Surface')
    pylab.plot(all_output[4][:],all_output[18][:],'b',label='Deep')
    if legend_counter==0:
        pylab.errorbar(t1,T1,yerr=er_T1,color='r',marker='o',linestyle="None")
        pylab.errorbar(t1b,T2,yerr=er_T2,color='b',marker='o',linestyle="None")
    Tdeep_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Cenozoic_tdeep_Beerling_Royer.txt',delimiter=',') #ultimaetly from Hansen as cited in Beerling and Royer
    pylab.plot(Tdeep_data[:,0],Tdeep_data[:,1]+273.15,color='b',marker='o',linestyle="None")
    ## suerface temperature from Hansen 2013
    hansen_surf=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/hansen2013_surf.txt',delimiter=',')
    #IMPORTANT: be aware Hansen surface is just multipliaction by deep sea change, with coefficient of 1. NOT TRUE RECORD
    if legend_counter==0:
        pylab.plot(hansen_surf[:,0],hansen_surf[:,1]+273.15,color='r',marker='o',linestyle="None")
        pylab.ylabel('Temperature (K)')
        pylab.xlabel('Time (yr)')
        pylab.legend(loc=2)
    
    
    #pylab.figure()
    pylab.subplot(3, 4, 7)
    pylab.plot(all_output[4][:],all_output[20][:],'b',label='Continental weathering')
    pylab.plot(all_output[4][:],all_output[21][:],'g',label='ocean precip.')
    if legend_counter==0:
        pylab.ylabel('Fluxes (mol C/yr)')
        pylab.xlabel('Time (yr)')     
        pylab.legend(loc=2)     
    
    pylab.subplot(3, 4, 8)
    pylab.plot(all_output[4][:],all_output[19][:],'r',label='Seafloor dissolution')
    pylab.plot(all_output[4][:],all_output[22][:],'k',label='pore precip.')
    #Gillis and Coogan estimates for precip are found in Notes2
    t_prec=numpy.array([98.*10**6])
    #get spread function
    F_outg=F_outgass+F_outgass*(t_prec/1e8) #cretacious based on Berner book Fig. 5.14         ###ITERATE THIS ### maybe
    spread= F_outg/F_outgass #bit of a fudge because not taking into account uncertainty in spreading rate
    spread=1
    #prec=numpy.array([spread*5.1])*confidence_Prec_p[1][0]
    #prec_er=2.7*confidence_Prec_p[1][0]
    # alternative  not using relative error
    prec=numpy.array([spread*2.35e12]) #changed from 2.3
    prec_er=0.75e12*spread # % changed from .64
    if legend_counter==0:
        pylab.errorbar(t_prec,prec,yerr=prec_er,color='k',marker='o',linestyle="None")
    
    #Estimates for dissolution change come from Coogan2015_error_Sr_diss.py in main file - basically using their Sr-fitted Ea with uncertainty, and uncertainty temperatures
    t_diss=numpy.array([98.*10**6])
    diss=numpy.array([spread*5.241])*all_output[19][0]
    diss_er=1.911*all_output[19][0]
    #pylab.errorbar(t_diss,diss,yerr=diss_er,color='r',marker='o',linestyle="None") #histogram better
    if legend_counter==0:
        pylab.legend(loc=2)
        pylab.ylabel('Fluxes (mol C/yr)')
        pylab.xlabel('Time (yr)')
    ### NOT CORRECT WAY OF DOING IT: really should be looking at distribution of relative change
    ### AND COMPARING TO CHANGE FROM SR iSOTOPES
    
    #pylab.figure()
    pylab.subplot(3, 4, 9)
    #pylab.hist(Dissolution_change_array)
    #pylab.xlabel('relative change in dissolution seafloor')
    #pylab.ylabel('Number')
    #pylab.errorbar(spread*5.241,10,xerr=spread*1.911,color='r',marker='o',linestyle="None")
    #pylab.title('Seafloor weathering change Sr isotopes')
    
    #pylab.figure()
    pylab.subplot(3, 4, 10)
    pylab.plot(all_output[4][:],all_output[0][:],'r--')
    pylab.plot(all_output[4][:],all_output[0][:],'r',label='DIC_o')
    pylab.plot(all_output[4][:],all_output[1][:],'b',label='ALK_o')
    pylab.plot(all_output[4][:],all_output[2][:],'g',label='DIC_p')
    pylab.plot(all_output[4][:],all_output[3][:],'k',label='ALK_p')
    if legend_counter==0:
        pylab.xlabel('Time (yr)')
        pylab.ylabel('Molality (mol/kg)')
        pylab.legend()
    
    pylab.tight_layout()


    #pylab.show()





def mc_plotter_spread(all_output,decision,legend_counter,spread_best): #identical except spreading change
    
    F_outgass=6e12 #dodgy, used in error calculations
    
    
    #pH from Pearson Palmer
    tpH1=numpy.array([.085,.98,1.49,3.0,3.31,3.87,6,6.2,9.02,10.39,11.4,11.81,13.06,14.73,14.96,16.23,16.7,18.38,19.85,21.7,23.,23.51])*10**6#,34.84,36.10,39.51])*10**6
    pH1=numpy.array([8.1,8.12,8.13,8.21,8.17,8.14,8.15,8.12,8.20,8.18,8.19,8.16,8.20,8.31,8.26,8.14,8.18,8.20,8.20,8.19,8.12,8.04])
    tpH2=numpy.array([40.12,42.52,44.26,45.69,46.07,46.97,50.33,51.02,52.22,53.24,55.84,57.12,59.88])*10**6
    pH2=numpy.array([7.8,8.07,7.95,7.79,7.54,7.99,7.84,7.92,7.42,7.62,7.48,7.54,7.42])
    #pylab.figure()
    if decision=="n":
        pylab.figure(figsize=(30,15))
    pylab.subplot(3, 4, 1)
    pylab.plot(all_output[4][:],all_output[5][:],'r',label='ocean')
    pylab.plot(all_output[4][:],all_output[7][:],'b',label='pore space')
    if legend_counter==0:
        pylab.plot(tpH1,pH1,'ro',linestyle="-")
        pylab.plot(tpH2,pH2,'ro',linestyle="-")
        pylab.xlabel('Time (yr)')
        pylab.ylabel('pH')
        pylab.legend()
    
    # CO2 
    #Modern CO2 for reference
    ppCO2=10**-6#0.000420318058799 #model modern
    preinudsmod=1.0#280.0
    #Cretaceous CO2: from Hong Lee 2012
    tCO2=numpy.array([65,65.5,66,66.5,67,68,68,75,76.5,80,83,91,95,98,100.5,102,103.5,107.5,108,113.5,115,115,120,122,125,129,143])*10**6
    CO2v=numpy.array([406,782,495,437,340,171,456,1412,656,917,1522,1437,1626,1520,1368,1428,1060,1219,907,449,1117,1325,798,1024,701,309,788])/preinudsmod
    CO2er=numpy.array([5,95,83,96,91,126,201,310,180,218,173,366,700,228,68,128,76,431,424,140,97,333,157,153,511,78,114])/preinudsmod
    #Cenozoic CO2 from Beerling Royer
    CO2_temp = numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Cenozoic_CO2_Beerling_Royer.txt',delimiter=',')
    #pylab.figure()
    pylab.subplot(3, 4, 2)
    pylab.plot(all_output[4][:],all_output[6][:]/ppCO2,'r',label='RCO2')
    if legend_counter==0:
        pylab.errorbar(tCO2,CO2v,yerr=CO2er,color='r',marker='o',linestyle="None")
        pylab.plot(CO2_temp[:,0],CO2_temp[:,1]/preinudsmod,color='r',marker='o',linestyle="None")
        pylab.xlabel('Time (yr)')
        pylab.ylabel('CO2 relative to modern')
        pylab.legend(loc=2)
    
        
    
      
    #pylab.figure()
    pylab.subplot(3, 4, 3)
    pylab.plot(all_output[4][:],all_output[9][:],'r',label='Ca ocean')
    pylab.plot(all_output[4][:],all_output[10][:],'b',label='Ca pore')
    #pylab.plot(all_output[4][:],Ca_fun(all_output[4][:]),'rx',label="Ca proxie Horita") #optional curve
    # Alternatively use data points from Horita table 2, and Cretaceous 94 Ma value from Timofeeff 2006
    tCa=numpy.array([5,14,35,37,94])*10**6
    Ca_prox_low=numpy.array([7,9,12,11,20])*10**-3
    Ca_prox_high=numpy.array([15,18,21,20,28])*10**-3
    Ca_prox=numpy.array([12,14,17,16,26])*10**-3
    if legend_counter==0:
        pylab.errorbar(tCa,Ca_prox,yerr=[Ca_prox-Ca_prox_low, Ca_prox_high-Ca_prox],color='r',marker='o',linestyle="None")
            
    pylab.plot(all_output[4][:],all_output[13][:],'k',label='HCO3 ocean')
    pylab.plot(all_output[4][:],all_output[14][:],'g',label='HCO3 pore')
    if legend_counter==0:
        pylab.xlabel('Molality (mol/kg)')
        pylab.xlabel('Time (yr)')
        pylab.legend(loc=2)
    
    #pylab.figure()
    pylab.subplot(3, 4, 4)
    # CO3 from Tyrrel and Zeebe (reconstructed from Ca and omega proxies)
    pylab.plot(all_output[4][:],all_output[11][:],'r',label='CO3 ocean')
    pylab.plot(all_output[4][:],all_output[12][:],'b',label='CO3 pore')
    CO3_data = numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Tyrrell_Zeebe_CO3.txt',delimiter=',')
    CO3_f = interp1d(numpy.flipud(CO3_data[:,0]),numpy.flipud(CO3_data[:,1]))
    if legend_counter==0:
        pylab.plot(all_output[4][:],CO3_f(all_output[4][:]),'ro')
        #pylab.plot(CO3_data[:,0],CO3_data[:,1],'ko')
        pylab.xlabel('Time (yr)')
        pylab.ylabel('Molality (mol/kg)')
        pylab.legend()    
            
    #omega and CCD (from Tyrrel and Zeebe)
    #pylab.figure()
    pylab.subplot(3, 4, 5)
    pylab.plot(all_output[4][:],all_output[15][:],'r',label='ocean')
    pylab.plot(all_output[4][:],all_output[16][:],'b',label='pore space')
    #do proxy stuff
    CCD_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Tyrrell_Zeebe_CCD.txt',delimiter=',')
    CCD_f=interp1d(numpy.flipud(CCD_data[:,0]),numpy.flipud(CCD_data[:,1]))
    # convert CCD to omega using equation 4 in Jansen et al. 2002, find coefficient first
    #k_CCD=all_output[15][0]/numpy.exp(0.176*(CCD_f(0)-3.06))
    k_CCD=all_output[15][0]/numpy.exp(0.189*(CCD_f(0)-3.82))
    CCD=CCD_f(all_output[4][:])
    #omega_proxy=k_CCD*numpy.exp(0.176*(CCD-3.06))
    omega_proxy = k_CCD*numpy.exp(0.189*(CCD-3.82))
    if legend_counter==0:
        pylab.plot(all_output[4][:],omega_proxy,'ro')
        pylab.legend(loc=2)
        pylab.ylabel('Saturation state')
        pylab.xlabel('Time (yr)')
    
    # Temperature Cretaceous plus pre-PETM
    t1=numpy.array([55.,67.5,75.,90.,100.])*10**6
    T1=numpy.array([23.029,19.03,21.56,26.68,25.35])+273.15 #surface
    er_T1=numpy.array([2.0,5.35,3.99,4.06,3.53])
    t1b=numpy.array([67.5,75.,90.,100.])*10**6
    T2=numpy.array([11,10,19,16])+273.15 #shallow
    er_T2=numpy.array([1,1,1,4])
    #pylab.figure()
    pylab.subplot(3, 4, 6)
    pylab.plot(all_output[4][:],all_output[17][:],'r',label='Surface')
    pylab.plot(all_output[4][:],all_output[18][:],'b',label='Deep')
    if legend_counter==0:
        pylab.errorbar(t1,T1,yerr=er_T1,color='r',marker='o',linestyle="None")
        pylab.errorbar(t1b,T2,yerr=er_T2,color='b',marker='o',linestyle="None")
    Tdeep_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Cenozoic_tdeep_Beerling_Royer.txt',delimiter=',') #ultimaetly from Hansen as cited in Beerling and Royer
    pylab.plot(Tdeep_data[:,0],Tdeep_data[:,1]+273.15,color='b',marker='o',linestyle="None")
    ## suerface temperature from Hansen 2013
    hansen_surf=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/hansen2013_surf.txt',delimiter=',')
    #IMPORTANT: be aware Hansen surface is just multipliaction by deep sea change, with coefficient of 1. NOT TRUE RECORD
    if legend_counter==0:
        pylab.plot(hansen_surf[:,0],hansen_surf[:,1]+273.15,color='r',marker='o',linestyle="None")
        pylab.ylabel('Temperature (K)')
        pylab.xlabel('Time (yr)')
        pylab.legend(loc=2)
    
    
    #pylab.figure()
    pylab.subplot(3, 4, 7)
    pylab.plot(all_output[4][:],all_output[20][:],'b',label='Continental weathering')
    pylab.plot(all_output[4][:],all_output[21][:],'g',label='ocean precip.')
    if legend_counter==0:
        pylab.ylabel('Fluxes (mol C/yr)')
        pylab.xlabel('Time (yr)')     
        pylab.legend(loc=2)     
    
    pylab.subplot(3, 4, 8)
    pylab.plot(all_output[4][:],all_output[19][:],'r',label='Seafloor dissolution')
    pylab.plot(all_output[4][:],all_output[22][:],'k',label='pore precip.')
    #Gillis and Coogan estimates for precip are found in Notes2
    t_prec=numpy.array([98.*10**6])
    #get spread function
    F_outg=F_outgass+F_outgass*(t_prec/1e8) #cretacious based on Berner book Fig. 5.14         ###ITERATE THIS ### maybe
    spread= F_outg/F_outgass #bit of a fudge because not taking into account uncertainty in spreading rate
    spread=spread_best
    #spread=1
    #prec=numpy.array([spread*5.1])*confidence_Prec_p[1][0]
    #prec_er=2.7*confidence_Prec_p[1][0]
    # alternative  not using relative error
    prec=numpy.array([spread*2.35e12]) #changed from 2.3
    prec_er=0.75e12*spread # % changed from 0.64
    if legend_counter==0:
        pylab.errorbar(t_prec,prec,yerr=prec_er,color='k',marker='o',linestyle="None")
    
    #Estimates for dissolution change come from Coogan2015_error_Sr_diss.py in main file - basically using their Sr-fitted Ea with uncertainty, and uncertainty temperatures
    t_diss=numpy.array([98.*10**6])
    diss=numpy.array([spread*5.241])*all_output[19][0]
    diss_er=1.911*all_output[19][0]
    #pylab.errorbar(t_diss,diss,yerr=diss_er,color='r',marker='o',linestyle="None") #histogram better
    if legend_counter==0:
        pylab.legend(loc=2)
        pylab.ylabel('Fluxes (mol C/yr)')
        pylab.xlabel('Time (yr)')
    ### NOT CORRECT WAY OF DOING IT: really should be looking at distribution of relative change
    ### AND COMPARING TO CHANGE FROM SR iSOTOPES
    
    #pylab.figure()
    pylab.subplot(3, 4, 9)
    #pylab.hist(Dissolution_change_array)
    #pylab.xlabel('relative change in dissolution seafloor')
    #pylab.ylabel('Number')
    #pylab.errorbar(spread*5.241,10,xerr=spread*1.911,color='r',marker='o',linestyle="None")
    #pylab.title('Seafloor weathering change Sr isotopes')
    
    #pylab.figure()
    pylab.subplot(3, 4, 10)
    pylab.plot(all_output[4][:],all_output[0][:],'r--')
    pylab.plot(all_output[4][:],all_output[0][:],'r',label='DIC_o')
    pylab.plot(all_output[4][:],all_output[1][:],'b',label='ALK_o')
    pylab.plot(all_output[4][:],all_output[2][:],'g',label='DIC_p')
    pylab.plot(all_output[4][:],all_output[3][:],'k',label='ALK_p')
    if legend_counter==0:
        pylab.xlabel('Time (yr)')
        pylab.ylabel('Molality (mol/kg)')
        pylab.legend()
    
    pylab.tight_layout()


    #pylab.show()

def dist_plotter(all_output,spread_output,sd):   
    
    #pH from Pearson Palmer
    tpH1=numpy.array([.085,.98,1.49,3.0,3.31,3.87,6,6.2,9.02,10.39,11.4,11.81,13.06,14.73,14.96,16.23,16.7,18.38,19.85,21.7,23.,23.51])*10**6#,34.84,36.10,39.51])*10**6
    pH1=numpy.array([8.1,8.12,8.13,8.21,8.17,8.14,8.15,8.12,8.20,8.18,8.19,8.16,8.20,8.31,8.26,8.14,8.18,8.20,8.20,8.19,8.12,8.04])
    tpH2=numpy.array([40.12,42.52,44.26,45.69,46.07,46.97,50.33,51.02,52.22,53.24,55.84,57.12,59.88])*10**6
    pH2=numpy.array([7.8,8.07,7.95,7.79,7.54,7.99,7.84,7.92,7.42,7.62,7.48,7.54,7.42])
    confidence_pH_o=scipy.stats.scoreatpercentile(all_output[:,5,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_pH_p=scipy.stats.scoreatpercentile(all_output[:,7,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    #pylab.figure()

    pylab.figure(figsize=(30,15))
    pylab.subplot(3, 4, 1)
    pylab.plot(all_output[0,4,:],confidence_pH_o[1],'r',label='ocean')
    pylab.plot(all_output[0,4,:],confidence_pH_p[1],'b',label='pore space')
    pylab.fill_between(all_output[0,4,:], confidence_pH_o[0], confidence_pH_o[2], color='red',  alpha=0.4)
    pylab.fill_between(all_output[0,4,:], confidence_pH_p[0], confidence_pH_p[2], color='blue', alpha=0.4)
    pylab.plot(tpH1,pH1,'ro',linestyle="-")
    pylab.plot(tpH2,pH2,'ro',linestyle="-")
    pylab.xlabel('Time (yr)')
    pylab.ylabel('pH')
    observ_pH_saved=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_pH.txt',delimiter=',')
    pylab.errorbar(observ_pH_saved[0,:],observ_pH_saved[1,:],observ_pH_saved[2,:],color='r')
    pylab.legend()

    ############################
    
    # CO2 
    #Modern CO2 for reference
    ppCO2=10**-6#0.000420318058799 #model modern
    preinudsmod=1.0#280.0
    #Cretaceous CO2: from Hong Lee 2012
    tCO2=numpy.array([65,65.5,66,66.5,67,68,68,75,76.5,80,83,91,95,98,100.5,102,103.5,107.5,108,113.5,115,115,120,122,125,129,143])*10**6
    CO2v=numpy.array([406,782,495,437,340,171,456,1412,656,917,1522,1437,1626,1520,1368,1428,1060,1219,907,449,1117,1325,798,1024,701,309,788])/preinudsmod
    CO2er=numpy.array([5,95,83,96,91,126,201,310,180,218,173,366,700,228,68,128,76,431,424,140,97,333,157,153,511,78,114])/preinudsmod
    #Cenozoic CO2 from Beerling Royer
    CO2_temp=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Cenozoic_CO2_Beerling_Royer.txt',delimiter=',')
    confidence_CO2o=scipy.stats.scoreatpercentile(all_output[:,6,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    #pylab.figure()
    pylab.subplot(3, 4, 2)
    pylab.plot(all_output[0,4,:],confidence_CO2o[1]/ppCO2,'r',label='RCO2')
    pylab.fill_between(all_output[0,4,:], confidence_CO2o[0]/ppCO2, confidence_CO2o[2]/ppCO2, color='red', alpha=0.4)    
    pylab.errorbar(tCO2,CO2v,yerr=CO2er,color='r',marker='o',linestyle="None")
    pylab.plot(CO2_temp[:,0],CO2_temp[:,1]/preinudsmod,color='r',marker='o',linestyle="None")
    pylab.xlabel('Time (yr)')
    pylab.ylabel('CO2 relative to modern')
    pylab.legend(loc=2)
    
    observ_CO2 = numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_CO2.txt',delimiter=',')
    pylab.errorbar(observ_CO2[0,:],observ_CO2[1,:],observ_CO2[2,:],color='r')
    
    #########################
    confidence_Ca_o=scipy.stats.scoreatpercentile(all_output[:,9,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Ca_p=scipy.stats.scoreatpercentile(all_output[:,10,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_CO3_o=scipy.stats.scoreatpercentile(all_output[:,11,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_CO3_p=scipy.stats.scoreatpercentile(all_output[:,12,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_HCO3_o=scipy.stats.scoreatpercentile(all_output[:,13,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_HCO3_p=scipy.stats.scoreatpercentile(all_output[:,14,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)        
    
      
    #pylab.figure()
    pylab.subplot(3, 4, 3)
    pylab.plot(all_output[0,4,:],confidence_Ca_o[1],'r',label='Ca ocean')
    pylab.plot(all_output[0,4,:],confidence_Ca_p[1],'b',label='Ca pore')
    pylab.fill_between(all_output[0,4,:], confidence_Ca_o[0], confidence_Ca_o[2], color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Ca_p[0], confidence_Ca_p[2], color='red', alpha=0.4)  
    #pylab.plot(all_output[4][:],Ca_fun(all_output[4][:]),'rx',label="Ca proxie Horita") #optional curve
    # Alternatively use data points from Horita table 2, and Cretaceous 94 Ma value from Timofeeff 2006
    tCa=numpy.array([5,14,35,37,94])*10**6
    Ca_prox_low=numpy.array([7,9,12,11,20])*10**-3
    Ca_prox_high=numpy.array([15,18,21,20,28])*10**-3
    Ca_prox=numpy.array([12,14,17,16,26])*10**-3
    pylab.errorbar(tCa,Ca_prox,yerr=[Ca_prox-Ca_prox_low, Ca_prox_high-Ca_prox],color='r',marker='o',linestyle="None")    
    pylab.plot(all_output[0,4,:],confidence_HCO3_o[1],'k',label='HCO3 ocean')
    pylab.plot(all_output[0,4,:],confidence_HCO3_p[1],'g',label='HCO3 pore')
    pylab.fill_between(all_output[0,4,:],confidence_HCO3_o[0],confidence_HCO3_o[2], color='grey', alpha=0.4)
    pylab.fill_between(all_output[0,4,:],confidence_HCO3_p[0],confidence_HCO3_p[2], color='green', alpha=0.4)   
    pylab.xlabel('Molality (mol/kg)')
    pylab.xlabel('Time (yr)')
    pylab.legend(loc=2)
    
    #pylab.figure()
    pylab.subplot(3, 4, 4)
    # CO3 from Tyrrel and Zeebe (reconstructed from Ca and omega proxies)
    pylab.plot(all_output[0,4,:],confidence_CO3_o[1],'r',label='CO3 ocean')
    pylab.plot(all_output[0,4,:],confidence_CO3_p[1],'b',label='CO3 pore')
    pylab.fill_between(all_output[0,4,:], confidence_CO3_o[0], confidence_CO3_o[2], color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_CO3_p[0], confidence_CO3_p[2], color='blue', alpha=0.4)  
    CO3_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Tyrrell_Zeebe_CO3.txt',delimiter=',')
    CO3_f=interp1d(numpy.flipud(CO3_data[:,0]),numpy.flipud(CO3_data[:,1]))
    pylab.plot(all_output[0,4,:],CO3_f(all_output[0,4,:]),'ro')
    #pylab.plot(CO3_data[:,0],CO3_data[:,1],'ko')
    pylab.xlabel('Time (yr)')
    pylab.ylabel('Molality (mol/kg)')
    pylab.legend()    
        
    #omega and CCD (from Tyrrel and Zeebe)
    #pylab.figure()
    pylab.subplot(3, 4, 5)
    confidence_omega_o=scipy.stats.scoreatpercentile(all_output[:,15,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_omega_p=scipy.stats.scoreatpercentile(all_output[:,16,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.plot(all_output[0,4,:],confidence_omega_o[1],'r',label='ocean')
    pylab.plot(all_output[0,4,:],confidence_omega_p[1],'b',label='pore space')
    pylab.fill_between(all_output[0,4,:], confidence_omega_o[0],confidence_omega_o[2], color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_omega_p[0], confidence_omega_p[2], color='blue', alpha=0.4)  
    #do proxy stuff
    CCD_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Tyrrell_Zeebe_CCD.txt',delimiter=',')
    CCD_f=interp1d(numpy.flipud(CCD_data[:,0]),numpy.flipud(CCD_data[:,1]))
    # convert CCD to omega using equation 4 in Jansen et al. 2002, find coefficient first
    #k_CCD=all_output[0,15,0]/numpy.exp(0.176*(CCD_f(0)-3.06))
    k_CCD=all_output[0,15,0]/numpy.exp(0.189*(CCD_f(0)-3.82))
    CCD=CCD_f(all_output[0,4,:])
    #omega_proxy=k_CCD*numpy.exp(0.176*(CCD-3.06))
    omega_proxy=k_CCD*numpy.exp(0.189*(CCD-3.82))
    pylab.plot(all_output[0,4,:],omega_proxy,'ro')
    pylab.legend(loc=2)
    pylab.ylabel('Saturation state')
    pylab.xlabel('Time (yr)')
    
    # Temperature Cretaceous plus pre-PETM
    t1=numpy.array([55.,67.5,75.,90.,100.])*10**6
    T1=numpy.array([23.029,19.03,21.56,26.68,25.35])+273.15 #surface
    er_T1=numpy.array([2.0,5.35,3.99,4.06,3.53])
    t1b=numpy.array([67.5,75.,90.,100.])*10**6
    T2=numpy.array([11,10,19,16])+273.15 #shallow
    er_T2=numpy.array([1,1,1,4])
    confidence_Tsurf=scipy.stats.scoreatpercentile(all_output[:,17,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Tdeep=scipy.stats.scoreatpercentile(all_output[:,18,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    #pylab.figure()
    pylab.subplot(3, 4, 6)
    pylab.fill_between(all_output[0,4,:], confidence_Tsurf[0], confidence_Tsurf[2], color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Tdeep[0], confidence_Tdeep[2], color='blue', alpha=0.4)  
    pylab.plot(all_output[0,4,:],confidence_Tsurf[1],'r',label='Surface')
    pylab.plot(all_output[0,4,:],confidence_Tdeep[1],'b',label='Deep')
    pylab.errorbar(t1,T1,yerr=er_T1,color='r',marker='o',linestyle="None")
    pylab.errorbar(t1b,T2,yerr=er_T2,color='b',marker='o',linestyle="None")
    Tdeep_data=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/Cenozoic_tdeep_Beerling_Royer.txt',delimiter=',') #ultimaetly from Hansen as cited in Beerling and Royer
    pylab.plot(Tdeep_data[:,0],Tdeep_data[:,1]+273.15,color='b',marker='o',linestyle="None")
    ## suerface temperature from Hansen 2013
    hansen_surf=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/hansen2013_surf.txt',delimiter=',')
    #IMPORTANT: be aware Hansen surface is just multipliaction by deep sea change, with coefficient of 1. NOT TRUE RECORD
    pylab.plot(hansen_surf[:,0],hansen_surf[:,1]+273.15,color='r',marker='o',linestyle="None")
    pylab.ylabel('Temperature (K)')
    pylab.xlabel('Time (yr)')
    pylab.legend(loc=2)
    
    
    #pylab.figure()
    confidence_Fd=scipy.stats.scoreatpercentile(all_output[:,19,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Fs=scipy.stats.scoreatpercentile(all_output[:,20,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Prec_o=scipy.stats.scoreatpercentile(all_output[:,21,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_Prec_p=scipy.stats.scoreatpercentile(all_output[:,22,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    pylab.subplot(3, 4, 7)
    pylab.plot(all_output[0,4,:], confidence_Fs[1],'b',label='Continental weathering')
    pylab.plot(all_output[0,4,:],confidence_Prec_o[1],'g',label='ocean precip.')
    pylab.fill_between(all_output[0,4,:], confidence_Fs[0], confidence_Fs[2], color='blue', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Prec_o[0], confidence_Prec_o[2], color='green', alpha=0.4)  
    pylab.ylabel('Fluxes (mol C/yr)')
    pylab.xlabel('Time (yr)')     
    pylab.legend(loc=2)     
    
    pylab.subplot(3, 4, 8)
    pylab.plot(all_output[0,4,:],confidence_Fd[1],'r',label='Seafloor dissolution')
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1],'k',label='pore precip.')
    pylab.fill_between(all_output[0,4,:], confidence_Fd[0], confidence_Fd[2], color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Prec_p[0], confidence_Prec_p[2], color='grey', alpha=0.4)  
    #Gillis and Coogan estimates for precip are found in Notes2
    t_prec=numpy.array([98.*10**6])
    
    #get spread function
    spread=1 #for no spreading dependence
    if sd=="y":
        spread=numpy.median(spread_output[:]) #bit of fudge to use median, but only for visualization
   
    
    #prec=numpy.array([spread*5.1])*confidence_Prec_p[1][0]
    #prec_er=2.7*confidence_Prec_p[1][0]
    # alternative  not using relative error
    prec=numpy.array([spread*2.35e12]) #changed 2.3
    prec_er=0.75e12*spread # % changed 0.64
    ## try this new thing
    ##spread_intervals=scipy.stats.scoreatpercentile(spread_output*2.3e12,[16,50,84], interpolation_method='fraction',axis=0)
    ##prec=spread_intervals[1]
    ##pre_er=
    ##

# alternative spread function
    spread_dist=[]
    for kk in range(0,10000):
        rate=1+numpy.random.uniform(0.2,1.5)
        beta_plot=numpy.random.uniform(0.0,1.0)
        #precip_e=(numpy.random.uniform(-.64,0.64)+2.3)*10**12
        precip_e=(numpy.random.uniform(-.75,0.75)+2.35)*10**12 #this is a range estimate remember
        #precip_e=(0.75*numpy.random.randn()+2.35)*10**12
        spread_dist.append(rate**beta_plot*precip_e) #### WHERE IS BETA?? Uh oh.   
    spread_dist=numpy.array(spread_dist)  
    [Slow,Smed,Shigh]=scipy.stats.scoreatpercentile(spread_dist,[5,50,95], interpolation_method='fraction',axis=0)    
    prec=Smed#0.5*(Slow+Shigh)
    #print "CHECK HERE",prec,Smed
    prec_er=numpy.array([[Smed-Slow],[Shigh-Smed]])#0.5*(Shigh-Slow)
# end alternative spread function    
    
    pylab.errorbar(t_prec,prec,yerr=prec_er,color='k',marker='o',linestyle="None")

    
    #Estimates for dissolution change come from Coogan2015_error_Sr_diss.py in main file - basically using their Sr-fitted Ea with uncertainty, and uncertainty temperatures
    t_diss=numpy.array([98.*10**6])
    diss=numpy.array([spread*5.241])*all_output[0,19,0]
    diss_er=1.911*all_output[0,19,0]
    #pylab.errorbar(t_diss,diss,yerr=diss_er,color='r',marker='o',linestyle="None") #histogram better
    pylab.legend(loc=2)
    pylab.ylabel('Fluxes (mol C/yr)')
    pylab.xlabel('Time (yr)')
    ### NOT CORRECT WAY OF DOING IT: really should be looking at distribution of relative change
    ### AND COMPARING TO CHANGE FROM SR iSOTOPES
    
    #pylab.figure()
    pylab.subplot(3, 4, 9)
    pylab.hist(numpy.mean(all_output[:,19,98:],axis=1)/all_output[:,19,0])
    pylab.xlabel('relative change in dissolution seafloor')
    pylab.ylabel('Number')
    pylab.errorbar(spread*5.241,10,xerr=spread*1.911,color='r',marker='o',linestyle="None")
    pylab.title('Seafloor weathering change Sr isotopes')
    
    #pylab.figure()
    pylab.subplot(3, 4, 10)
    confidence_DICo=scipy.stats.scoreatpercentile(all_output[:,0,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_ALKo=scipy.stats.scoreatpercentile(all_output[:,1,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_DICp=scipy.stats.scoreatpercentile(all_output[:,2,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)
    confidence_ALKp=scipy.stats.scoreatpercentile(all_output[:,3,:],[2.5,50,97.5], interpolation_method='fraction',axis=0)

    pylab.plot(all_output[0,4,:],confidence_DICo[1],'r',label='DIC_o')
    pylab.plot(all_output[0,4,:],confidence_ALKo[1],'b',label='ALK_o')
    pylab.plot(all_output[0,4,:],confidence_DICp[1],'g',label='DIC_p')
    pylab.plot(all_output[0,4,:],confidence_ALKp[1],'k',label='ALK_p')
    
    pylab.fill_between(all_output[0,4,:], confidence_DICo[0], confidence_DICo[2], color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_ALKo[0], confidence_ALKo[2], color='blue', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_DICp[0], confidence_DICp[2], color='green', alpha=0.4)  
    #pylab.fill_between(all_output[0,4,:], confidence_ALKp[0], confidence_ALKp[2], color='grey', alpha=0.4)  
    
    pylab.xlabel('Time (yr)')
    pylab.ylabel('Molality (mol/kg)')
    pylab.legend()
    
    pylab.subplot(3, 4, 11)
    pylab.hist(numpy.mean(all_output[:,20,98:],axis=1)/all_output[:,20,0])
    pylab.xlabel('relative change in silicate weathering')
    pylab.ylabel('Number')
    #pylab.errorbar(spread*5.241,10,xerr=spread*1.911,color='r',marker='o',linestyle="None")
    #pylab.title('')
    pylab.tight_layout()
    
    
    
    ########################################
    ## FIGURE FOR PAPER
    ########################################
    
    all_output[0,4,:]=all_output[0,4,:]/1e6
    
    strt_lim=-.01e2
    fin_lim=1.01e2
    pylab.figure(figsize=(30,15))
    pylab.subplot(3, 3, 1)
    pylab.plot(all_output[0,4,:],confidence_pH_o[1],'k',label='ocean')
    pylab.fill_between(all_output[0,4,:], confidence_pH_o[0], confidence_pH_o[2], color='grey', alpha=0.4)
    #pylab.plot(tpH1/1e6,pH1,'ko',linestyle="--")
    #pylab.plot(tpH2/1e6,pH2,'ko',linestyle="--")
    ##
    observ_pH_saved=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_pH.txt',delimiter=',')
    pylab.errorbar(observ_pH_saved[0,:],observ_pH_saved[1,:],observ_pH_saved[2,:],color='k',marker='o',linestyle="None")
    ##
    pylab.xlabel('Time (Ma)')
    pylab.ylabel('Ocean pH')
    #pylab.legend()
    pylab.xlim([strt_lim,fin_lim])
    pylab.text(-10, 8.35, 'A', #transform=ax.transAxes,
      fontsize=16, fontweight='bold', va='top')

    pylab.subplot(3, 3, 2)
    pylab.plot(all_output[0,4,:],confidence_CO2o[1]/ppCO2,'k',label='Atmospheric CO2')
    pylab.fill_between(all_output[0,4,:], confidence_CO2o[0]/ppCO2, confidence_CO2o[2]/ppCO2, color='grey', alpha=0.4)    
    #pylab.errorbar(tCO2/1e6,CO2v,yerr=CO2er,color='k',marker='o',linestyle="None")
    #pylab.plot(CO2_temp[:,0]/1e6,CO2_temp[:,1]/preinudsmod,color='k',marker='o',linestyle="None")
    ##
    observ_CO2=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_CO2.txt',delimiter=',')
    pylab.errorbar(observ_CO2[0,:],observ_CO2[1,:],observ_CO2[2,:],color='k',marker='o',linestyle="None")
    ##
    pylab.xlabel('Time (Ma)')
    pylab.ylabel('Atmospheric pCO2 (ppm)')
    #pylab.legend(loc=2)
    pylab.xlim([strt_lim,fin_lim])
    pylab.text(-10, 3400, 'B', #transform=ax.transAxes,
      fontsize=16, fontweight='bold', va='top')

    pylab.subplot(3,3 ,3)
    pylab.plot(all_output[0,4,:],confidence_omega_o[1],'k',label='ocean')
    pylab.fill_between(all_output[0,4,:], confidence_omega_o[0],confidence_omega_o[2], color='grey', alpha=0.4)  
    #do proxy stuff
    #pylab.plot(all_output[0,4,:],omega_proxy,'ko')
    ##
    observ_omega=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_omega_calc.txt',delimiter=',') #calcite
    pylab.errorbar(observ_omega[0,:],observ_omega[1,:],observ_omega[2,:],color='k',marker='o',linestyle="None")
    ##
    #pylab.legend(loc=2)
    pylab.ylabel('Saturation state')
    pylab.xlabel('Time (Ma)')
    pylab.xlim([strt_lim,fin_lim])
    pylab.text(-10, 3.58, 'C', #transform=ax.transAxes,
      fontsize=16, fontweight='bold', va='top')
    
    pylab.subplot(3, 3, 4)
    pylab.fill_between(all_output[0,4,:], confidence_Tsurf[0], confidence_Tsurf[2], color='grey', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Tdeep[0], confidence_Tdeep[2], color='red', alpha=0.4)  
    pylab.plot(all_output[0,4,:],confidence_Tsurf[1],'k',label='Surface')
    pylab.plot(all_output[0,4,:],confidence_Tdeep[1],'r',label='Deep ocean')
    #pylab.errorbar(t1/1e6,T1,yerr=er_T1,color='k',marker='o',linestyle="None")
    #pylab.errorbar(t1b/1e6,T2,yerr=er_T2,color='r',marker='o',linestyle="None")
    #pylab.plot(Tdeep_data[:,0]/1e6,Tdeep_data[:,1]+273.15,color='r',marker='o',linestyle="None")
    ## suerface temperature from Hansen 2013
    ##
    observ_Td=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_Td.txt',delimiter=',')
    pylab.errorbar(observ_Td[0,:],observ_Td[1,:]+273.15,observ_Td[2,:],color='r',marker='o',linestyle="None")#time_CO2
    observ_T=numpy.loadtxt('/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/obs_T.txt',delimiter=',')
    pylab.errorbar(observ_T[0,:],observ_T[1,:]+273.15,observ_T[2,:],color='k',marker='o',linestyle="None")#time_CO2    
    ##
    pylab.ylabel('Temperature (K)')
    pylab.xlabel('Time (Ma)')
    pylab.legend(loc=2)
    pylab.xlim([strt_lim,fin_lim])
    pylab.text(-10, 304, 'D', #transform=ax.transAxes,
      fontsize=16, fontweight='bold', va='top')
    
    pylab.subplot(3,3,5)
    pylab.plot(all_output[0,4,:], confidence_Fs[1]/1e12,'r',label='Cont. weathering')
    pylab.plot(all_output[0,4,:],confidence_Prec_o[1]/1e12,'k',label='ocean precip.')
    pylab.fill_between(all_output[0,4,:], confidence_Fs[0]/1e12, confidence_Fs[2]/1e12, color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Prec_o[0]/1e12, confidence_Prec_o[2]/1e12, color='grey', alpha=0.4)  
    pylab.ylabel('Fluxes (Tmol/yr)')
    pylab.xlabel('Time (Ma)')     
    pylab.legend(loc=2)  
    pylab.xlim([strt_lim,fin_lim]) 
    pylab.text(-10, 24, 'E', #transform=ax.transAxes,
      fontsize=16, fontweight='bold', va='top')
    
    pylab.subplot(3, 3, 6)
    pylab.plot(all_output[0,4,:],confidence_Fd[1]/1e12,'r',label='Seafloor dissolution')
    pylab.plot(all_output[0,4,:],confidence_Prec_p[1]/1e12,'k',label='pore precip.')
    pylab.fill_between(all_output[0,4,:], confidence_Fd[0]/1e12, confidence_Fd[2]/1e12, color='red', alpha=0.4)  
    pylab.fill_between(all_output[0,4,:], confidence_Prec_p[0]/1e12, confidence_Prec_p[2]/1e12, color='grey', alpha=0.4)  
    pylab.errorbar(t_prec/1e6,prec/1e12,yerr=prec_er/1e12,color='k',marker='o',linestyle="None")
    pylab.legend(loc=2) 
    pylab.ylabel('Fluxes (Tmol/yr)')
    pylab.xlabel('Time (Ma)')   
    pylab.xlim([strt_lim,fin_lim])
    pylab.text(-10, 5.9, 'F', fontsize=16, fontweight='bold', va='top')
    
    ## old figures for relative change
    #pylab.subplot(3, 3, 7)
    #pylab.hist(numpy.mean(all_output[:,19,98:],axis=1)/all_output[:,19,0], bins=30,color='grey',density=True)
    #pylab.xlabel('relative change dissolution seafloor')
    #pylab.ylabel('Probability density')
    ##pylab.errorbar(spread*5.241,10,xerr=spread*1.911,color='r',marker='o',linestyle="None")
    ##pylab.title('Seafloor weathering change Sr isotopes')
    dist_1=numpy.mean(all_output[:,19,98:],axis=1)/all_output[:,19,0]
    print (numpy.percentile(dist_1,2.5),numpy.percentile(dist_1,50),numpy.percentile(dist_1,97.5))
    ## pylab.text(-0.15, 1.8, 'G', fontsize=16, fontweight='bold', va='top')

    #
    #pylab.subplot(3, 3, 8)
    #pylab.hist(numpy.mean(all_output[:,20,98:],axis=1)/all_output[:,20,0],bins=30,color='grey',density=True)
    #pylab.xlabel('relative change silicate weathering')
    #pylab.ylabel('Probability density')
    dist_2=numpy.mean(all_output[:,20,98:],axis=1)/all_output[:,20,0]
    print (numpy.percentile(dist_2,2.5),numpy.percentile(dist_2,50),numpy.percentile(dist_2,97.5))
    ## pylab.text(0.45, 13.5, 'H', fontsize=16, fontweight='bold', va='top')


    pylab.subplot(3, 3, 8)
    #pylab.hist2d(numpy.mean(all_output[:,20,98:],axis=1)/all_output[:,20,0], numpy.mean(all_output[:,19,98:],axis=1)/all_output[:,19,0], bins=30,density=True,cmap=pylab.cm.jet)
    #pylab.colorbar(label='Probability density')
    pylab.xlabel('Relative change continental silicate weathering')
    pylab.ylabel('Relative change seafloor weathering')
    pylab.show()
    
    pylab.tight_layout()
