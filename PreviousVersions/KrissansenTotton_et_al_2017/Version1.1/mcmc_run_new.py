import emcee
import numpy
import matplotlib.pyplot as pl
from corner import corner
import corner
from model_functions import forward_model
#from foward_model_for_mcmc import forward_model_old

### Data constraints:
###########################
## Surface Temp
observ_T=numpy.loadtxt('obs_T.txt',delimiter=',')
time_T=observ_T[0,:]
time_T=numpy.array(time_T.astype(int))
obs_T=observ_T[1,:]+273.15
er_T=observ_T[2,:]

#######################
## Deep ocean Temp
observ_Td=numpy.loadtxt('obs_Td.txt',delimiter=',')
time_Td=observ_Td[0,:]#time_CO2
time_Td=numpy.array(time_Td.astype(int))
obs_Td=observ_Td[1,:]+273.15
er_Td=observ_Td[2,:]
##############################

## pCO2 constraints
preinudsmod=280.0
observ_CO2=numpy.loadtxt('obs_CO2.txt',delimiter=',')
time_CO2=observ_CO2[0,:]#time_CO2
time_CO2=numpy.array(time_CO2.astype(int))
obs_CO2=observ_CO2[1,:]/preinudsmod
er_CO2=observ_CO2[2,:]/preinudsmod

########################
# ocean carbonate precipitation
obs_prec_zero=2.3e12 #need spread modifier later
obs_prec_zero=2.35e12 #need spread modifier later
time_prec=99
er_prec_zero=0.64e12 # need spread modifier late
er_prec_zero=0.75e12 # need spread modifier late
########################

## ocean saturation state
#observ_omega=numpy.loadtxt('obs_omega.txt',delimiter=',') #aragonite
observ_omega=numpy.loadtxt('obs_omega_calc.txt',delimiter=',') #calcite
time_omega=observ_omega[0,:]#time_CO2
time_omega=time_omega.astype(int)
obs_omega_o=observ_omega[1,:]
er_omega_o=observ_omega[2,:]
#########################

#pH constraints
observ_pH_saved=numpy.loadtxt('obs_pH.txt',delimiter=',')
time_pH=observ_pH_saved[0,:]#numpy.array([5,15,35,45,55])
time_pH=time_pH.astype(int)
obs_pH=observ_pH_saved[1,:]
er_pH=observ_pH_saved[2,:]
##########################

### This is the likelihood function:
def LnLike(x):
  #import pdb
  #pdb.set_trace()
  ### Each range defines the prior for unknown variables that wea are trying to determine
  if (x[0] > 0.5 ) or (x[0]<0.2): # exponent for pCO2-dependence continental weeathering (default)
  #if (x[0] > 1.0 ) or (x[0]<0.01): # alternative exponent for plant modified weathering 
  #if (x[0] > 0.045 ) or (x[0]<0.025): # alternative exponent for runoff parameterization
    return -numpy.inf    ## likelihood is -infinity if parameter is outside prior range
  if (x[1] > 50 ) or (x[1]<5): #Te parameter (K)
    return -numpy.inf    
  if (x[2] > 1.2 ) or (x[2] < 0.2): # This is the weatherability parameter, W+1 (so W: -0.8 to 1.2)
    return -numpy.inf 
  if (x[3] > 8.0) or (x[3] < 1.5): # climate sensitivity 
    return -numpy.inf 
  if (x[4] > 1.5) or (x[4] < 0.2): ### Volcanic outgassing parameter, V (I think this is V not V+1)
    return -numpy.inf 
  if (x[5] > 1.5) or (x[5] < -0.9): ## Carbonate weatherability parameter
    return -numpy.inf 
  if (x[6] > 10e12) or (x[6] < 4e12): #modern outgassing, Tmol C/yr
    return -numpy.inf 
  if (x[7] > 14e12) or (x[7] < 7e12): #modern carbonate weathering, Tmol C/yr
    return -numpy.inf 
  if (x[8] > 1e6) or (x[8] < 2e4): # Mixing time for pore-space, yrs. Mixing flux will be Mass_ocean / x[8]
    return -numpy.inf 
  if (x[9] > 2.5) or (x[9] < 1.0): # Exponent for carbonate precipitation, dependence on saturation state
    return -numpy.inf 
  if (x[10] > 1.5) or (x[10] < 0.5): # fraction of carbonate precipitation seafloor relative to alkalinity release of seafloor? 
    return -numpy.inf 
  if (x[11] > 1.4) or (x[11] < 0.8): # gradient deep oecan temperature to surface temperature
    return -numpy.inf 
  if (x[12] > 0.5) or (x[12] < 0.0): #dissolution exponent for seafloor weathering pH-dependence
    return -numpy.inf
  if (x[13] > 110000) or (x[13] < 40000): #activation energy seafloor weathering
    return -numpy.inf 
  if (x[14] > 0.6) or (x[14] < 0.4): #pelagic carbonate fraction
    return -numpy.inf 
  if (x[15] > 1.0) or (x[15] < 0.0): #" beta" parameter for outgassing dependence seafloor weathering
    return -numpy.inf 
  if (x[16] > 5.0) or (x[16] < 0.0): # Paleogeography fudge factor
    return -numpy.inf 
    
  # Now try and run the parameters for the parameters above.
  try: #climp,tdep_weath,lfrac,cl_sens,change_out,CWF,F_outgass,F_carbw,tau_mix,n,alt_frac,deep_grad,coef_for_diss,Eact,fpel,beta,PG)
         #x0     1          2     3       4         5     6        7      8      9  10       11        12            13  14 15     16                                                                             
    [[out0,out1,out2,out3,time,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array],imbalance]=forward_model(x[8],x[6],x[9],x[0],x[1],0.45e12,x[10],0.01,x[2],x[3],x[4],x[7],x[14],x[5],x[11],x[12],x[15],x[13],x[16])
    #[out0,out1,out2,out3,time,pH_array_o,CO2_array_o,pH_array_p,CO2_array_p,Ca_array_o,Ca_array_p,CO3_array_o,CO3_array_p,HCO3_array_o,HCO3_array_p,omega_o,omega_p,Tsurf_array,Tdeep_array,Fd_array,Fs_array,Precip_ocean_array,Precip_pore_array]=forward_model_old(x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16])
    #print(out0[0],out00[0])
    #print(pH_array_p[0],pH_array_p0[0])
    #print(Fs_array[0],Fs_array0[0])
    #print(Fd_array[0],Fd_array0[0])
    #import pdb
    #pdb.set_trace()
    #  W,F_outgass,    n, CO2_dep,Te,mod_sea,alt_frac,Mp_frac,W_plus_1,cl_sens,change_out,F_carbw,frac_pel,CWF,deep_grad,coef_for_diss,beta,Ebas,PG
    #  x[8], x[6], x[9],x[0],x[1],0.45e12,x[10], 0.01,   x[2],     x[3],  x[4],        x[7],   x[14], x[5], x[11],    x[12],       x[15], x[13], x[16]
  except:
    return -numpy.inf ## if it fails, throw it away

  # dont allow unbounded or stupid solutions 
  if numpy.isnan(pH_array_p[98]):
    return -numpy.inf 
  if HCO3_array_p[98]<0.0:
    return -numpy.inf 
  if out2[98]<0.0:
    return -numpy.inf   
  
  spread= ( (x[6]+x[4]*x[6])/x[6] )**x[15] ## Spreading rate at 100 Myr for this particular run
  ## will be used to calculate model carbonate precipitation in seafloor at 100 Myr below
  
  ### Now get outputs from model at the times at which we have data to compare
  model_T=Tsurf_array[time_T]
  model_Td=Tdeep_array[time_Td]
  model_CO2=CO2_array_o[time_CO2]/CO2_array_o[0] 
  model_prec=numpy.mean(Precip_pore_array[time_prec]) 
  model_omega_o=omega_o[time_omega] 
  model_pH=pH_array_o[time_pH]

  # Also need to use spreading rates from 100 Myr to get "observed" precipitation data
  obs_prec=obs_prec_zero*spread
  er_prec=er_prec_zero*spread
  
  ## Now calculate the log-likelihood i.e. compare data the model outputs
  ## The larger log_like, the better the fit
  log_terms=numpy.sum(numpy.log(2*numpy.pi*er_T**2))+numpy.sum(numpy.log(2*numpy.pi*er_Td**2)) +numpy.sum(numpy.log(2*numpy.pi*er_omega_o**2)) +numpy.log(2*numpy.pi*er_prec**2) +numpy.sum(numpy.log(2*numpy.pi*er_pH**2))   + numpy.sum(numpy.log(2*numpy.pi*er_CO2**2))  
  log_like =  -0.5 * (   numpy.sum((obs_T-model_T)**2/er_T**2) + numpy.sum((obs_Td-model_Td)**2/er_Td**2) + (obs_prec-model_prec)**2/er_prec**2 +  numpy.sum((obs_CO2-model_CO2)**2/er_CO2**2)  + 
                         numpy.sum((obs_omega_o-model_omega_o)**2/er_omega_o**2) + log_terms +  numpy.sum((obs_pH-model_pH)**2/er_pH**2) )               
  return log_like # return the likelihood value to the program running the inverse analysis


# Set up the sampler
ndim = 17 ### This is the number of parameters we are solving for (x[0], x[1], ... x[16] = 17 parameters)
### Next two parameters control the number of walks and the number of steps to take when solving inverse problem
### nwalk * nsteps = total number of model runs
### To do a full model run, nwalk = 1000 and nsteps = 10000 is appropriate 
### For a quick run, nwalk = 100 and nsteps = 1000 will tell you if the code is working.
# nwalk  = 500 
# nsteps = 2000
# nwalk = 200
# nsteps = 500
### --- quick test run values --- ###
nwalk  = 50
nsteps = 500 
### Important: if you make nsteps < 1000, you will need to modify some of the plotting stuff below
### This is because, ideally, you want to throw out the first 1000 steps, as it takes a while for the Markov
### chains to converge on the posteriors. But if you are just doing a quick test run to see if the code works
### you might have <1000 steps, and so you definitely don't want to throw away your whole sample!


### Define the initial random valules for the 17 variables. This is where the walkers start.
#p0 = numpy.vstack([[0.01+0.99*numpy.random.random() for i in range(nwalk)], #for plant pCO2 dependence 0-1 WEATH1
#p0 = numpy.vstack([[0.025+0.02*numpy.random.random() for i in range(nwalk)], #for runoff no pCo2 dependence WEATH2
p0 = numpy.vstack([[0.2+0.3*numpy.random.random() for i in range(nwalk)], #Default
                 [5+45*numpy.random.random() for i in range(nwalk)],
                [0.2+1.0*numpy.random.random() for i in range(nwalk)],
                [1.5+6.5*numpy.random.random() for i in range(nwalk)],
                [0.2+1.3*numpy.random.random() for i in range(nwalk)],
                [-0.9+2.4*numpy.random.random() for i in range(nwalk)],
                [4e12+6e12*numpy.random.random() for i in range(nwalk)],
                [7e12+7e12*numpy.random.random() for i in range(nwalk)],
                [20000+980000.*numpy.random.random() for i in range(nwalk)],
                [1.0+1.5*numpy.random.random() for i in range(nwalk)],
                [0.5+1.0*numpy.random.random() for i in range(nwalk)],
                [0.8+0.6*numpy.random.random() for i in range(nwalk)],
                [0.+0.5*numpy.random.random() for i in range(nwalk)],
                [40000.+70000.*numpy.random.random() for i in range(nwalk)],
                [0.4+0.2*numpy.random.random() for i in range(nwalk)],
                [1.0*numpy.random.random() for i in range(nwalk)],
                [5.0*numpy.random.random() for i in range(nwalk)]]).T
                
## Actually do the inverse analysis:
## Note that the "threads" variable below is very important for code parallelization.
## This needs to be chosen carefully to reflect the number of cores on your computer.
## If you can max out all your cores, that would be ideal. 
## But if threads > number of cores then it'll run slow.

#from multiprocessing import Pool
#with Pool() as pool:
sampler = emcee.EnsembleSampler(nwalk, ndim, LnLike,threads=14)
                                #args = (t, y, error), 
                                #pool = emcee.interruptible_pool.InterruptiblePool()) #for parallel
pos, lnprob, rstate=sampler.run_mcmc(p0, nsteps)

## Save the ouputs just in case the plotting code crashes (don't want to have to redo the inverse analysis!)
numpy.save('newchain2',sampler.chain[:,:,:]) 
numpy.save('newln2',sampler.lnprobability) 

## autocorrelation parameter (useful for checking if you did enough model runs to have valid posteriors)
try:
    print ("ESS",nsteps * nwalk / numpy.nanmax(sampler.acor))
except:
    print ("couldn't calculate autocorrelation")

chain=numpy.load('newchain2.npy') 
lnprob=numpy.load('newln2.npy') 

##### The rest is all plotting

# Plot a few individual chains
fig, ax = pl.subplots(4) 
for n in range(nwalk):
  ax[0].plot(chain[n,:,0])
  ax[1].plot(chain[n,:,1])
  ax[2].plot(chain[n,:,2]) 
  ax[3].plot(chain[n,:,3]) 

#find highest likelihood run
logprob=numpy.array(lnprob)
values=chain
ii,jj = numpy.unravel_index(logprob.argmax(), logprob.shape)
print ("indeces for best",ii,jj)
print ("loglikelihood and values",logprob[ii,jj],values[ii,jj,:])

# Plot the corner plot, discarding the first 1000 steps as burn-in
production = chain[:,1000:,:]
#production = chain[:,0:,:] ## Use this if you have <1000 steps
s = production.shape
flatchain = production.reshape(s[0] * s[1], s[2])

flatchain2=numpy.copy(flatchain)
## convert some variables:
flatchain2[:,4]=flatchain2[:,4]+1 #make outgassing relative Cretaceous outgassing, V+1
## Weatherability is already Cretaceous weatherability, W+1, assuming w=1+W with plus sign
flatchain2[:,5]=flatchain2[:,5]+1 #Carbonate weathering modifier
flatchain2[:,8]=flatchain2[:,8]/1000.0 #Convert to ky
flatchain2[:,13]=flatchain2[:,13]/1000.0 #Convert to kJ/mol
flatchain2[:,6]=flatchain2[:,6]/1e12 #Convert to Tmol
flatchain2[:,7]=flatchain2[:,7]/1e12 #Convert to Tmol

from matplotlib import rc
## Plot posteriors as corner plots (compare to Fig. 6 in the paper)
corner.corner(flatchain2[:,[0,1,2,3,4]], quantiles=[0.16, 0.5, 0.84],labels=[r"CO$_2$-dependence, $\alpha$", "Temp. dep. cont.\nweathering, $T_e$ (K)", "Relative Cretaceous\nweatherability, 1+$W$","Climate sensitivity,\n${\Delta}T_{2X}$ (K)","Relative Cretaceous\noutgassing, 1+$V$"])#,truths=values[ii,jj,:])
corner.corner(flatchain2[:,[5,8,11,13,6,16]], quantiles=[0.16, 0.5, 0.84],labels=["Carbonate weath.\nmodifier, 1+$C_{WF}$",r"Circulation time, $\tau$ (kyr)","Surface-deep\ntemp. gradient, $a_{grad}$","Temp. dependence\nseafloor, $E_{bas}$ (kJ/mol)","Modern outgassing,\n$F_{out}^{mod}$ (Tmol C/yr)","Paleogeography parameter,\n${\Delta}P$ (K)"])#,truths=values[ii,jj,:])
corner.corner(flatchain2[:,[7,9,10,12,14,15]], quantiles=[0.16, 0.5, 0.84],labels=["Modern carb.\nweathering, $F_{carb}^{mod}$ (Tmol C/yr)","Carb. precip.\ncoefficient, $n$","Modern seafloor diss.\nrelative precip.","pH dependence\nseafloor, $\gamma$","Modern pelagic\nfraction",r"Spreading rate dep., $\beta$"])#,truths=values[ii,jj,:])


## Confidence intervals for unknown parameter:
ab, bc, cd,de,ef ,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr= map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),zip(*numpy.percentile(flatchain, [16, 50, 84],axis=0)))
print ("median values with errors", numpy.array([ab, bc, cd,de,ef,fg,gh,hi,ij,jk,kl,lm,mn,no,op,pq,qr]))
print ("confidence intervals")
map(lambda v: (v[0], v[1], v[2]),zip(*numpy.percentile(flatchain, [5, 50, 95],axis=0)))

from plotting_everything import mc_plotter_spread,dist_plotter

## Can't remember what this does - probably not important
import pylab
pylab.figure(figsize=(30,15))
legend_counter=0
for x_ex in flatchain[numpy.random.randint(len(flatchain), size=100)]:
    #print (x_ex)
    #outputs=forward_model_old(x_ex[0],x_ex[1],x_ex[2],x_ex[3],x_ex[4],x_ex[5],x_ex[6],x_ex[7],x_ex[8],x_ex[9],x_ex[10],x_ex[11],x_ex[12],x_ex[13],x_ex[14],x_ex[15],x_ex[16])
    [outputs,imbalance] = forward_model(x_ex[8],x_ex[6],x_ex[9],x_ex[0],x_ex[1],0.45e12,x_ex[10],0.01,x_ex[2],x_ex[3],x_ex[4],x_ex[7],x_ex[14],x_ex[5],x_ex[11],x_ex[12],x_ex[15],x_ex[13],x_ex[16])
    
    sp=((x_ex[6]+x_ex[4]*x_ex[6])/x_ex[6])**x_ex[15]
    mc_plotter_spread(outputs,"y",legend_counter,sp)
    legend_counter=legend_counter+1

### This is important. This takes 1000 sets of parameter values from your posterior
### and re-runs the forward model 1000 times to get distributions for the time-evolution
### of different model parameters e.g. Fig. 5a-f
### Another script is called to actually do the plotting
mega_output=[]
spread_output=[]
carbw_factor=[]
sil_change=[]
seafloor_change=[]
for x_ex in flatchain[numpy.random.randint(len(flatchain), size=1000)]:
    #print (x_ex)
    #outputs=forward_model_old(x_ex[0],x_ex[1],x_ex[2],x_ex[3],x_ex[4],x_ex[5],x_ex[6],x_ex[7],x_ex[8],x_ex[9],x_ex[10],x_ex[11],x_ex[12],x_ex[13],x_ex[14],x_ex[15],x_ex[16])
    [outputs,imbalance]= forward_model(x_ex[8],x_ex[6],x_ex[9],x_ex[0],x_ex[1],0.45e12,x_ex[10],0.01,x_ex[2],x_ex[3],x_ex[4],x_ex[7],x_ex[14],x_ex[5],x_ex[11],x_ex[12],x_ex[15],x_ex[13],x_ex[16])
    spread_output.append( ((x_ex[6]+x_ex[4]*x_ex[6])/x_ex[6])**x_ex[15])
    mega_output.append(outputs)
    carbw_factor.append((1+x_ex[5])*outputs[20][99]/outputs[20][0])
    sil_change.append(outputs[20][99]-outputs[20][0])
    seafloor_change.append(outputs[19][99]-outputs[19][0])
mega_output=numpy.array(mega_output)
spread_output=numpy.array(spread_output)
dist_plotter(mega_output,spread_output,"y")
sil_change=numpy.array(sil_change)/1e12
seafloor_change=numpy.array(seafloor_change)/1e12


change_precip_array=mega_output[:,21,99]/mega_output[:,21,0]
print (numpy.percentile(numpy.array(change_precip_array),2.5),numpy.percentile(numpy.array(change_precip_array),50),numpy.percentile(numpy.array(change_precip_array),97.5))
print (numpy.percentile(numpy.array(change_precip_array),16),numpy.percentile(numpy.array(change_precip_array),50),numpy.percentile(numpy.array(change_precip_array),84))

## The rest has something to do with Fig. 5g-i,

#final subplot  
pylab.subplot(3, 3, 9)
pylab.hist2d(sil_change, seafloor_change, range=[[-1, 6], [0, 6]],bins=30,normed=True,cmap=pylab.cm.jet)
pylab.colorbar(label='Probability density')
pylab.xlabel('Decrease in continental weathering\nsince mid Cretaceous (Tmol/yr)')
pylab.ylabel('Decrease in seafloor weathering\nsince mid Cretaceous (Tmol/yr)')

# carbonate plot:
pylab.subplot(3, 3, 7)
pylab.hist(numpy.array(carbw_factor),bins=30,color='grey',normed=True)
pylab.xlabel('Relative change carbonate weathering')
pylab.ylabel('Probability density')
print (numpy.percentile(numpy.array(carbw_factor),2.5),numpy.percentile(numpy.array(carbw_factor),50),numpy.percentile(numpy.array(carbw_factor),97.5))

pylab.show()


pylab.figure()
pylab.hist(sil_change,bins=30,color='grey',normed=True)
pylab.xlabel('Absolute change silicate weathering')
pylab.ylabel('Probability density')

pylab.figure()
pylab.hist(seafloor_change,bins=30,color='grey',normed=True)
pylab.xlabel('Absolute change seafloor weathering')
pylab.ylabel('Probability density')

pylab.figure()
pylab.hist2d(sil_change, seafloor_change, range=[[-1, 6], [-1, 6]],bins=30,normed=True,cmap=pylab.cm.jet)
pylab.colorbar(label='Probability density')
pylab.xlabel('Decrease in continental silicate weathering\nsince mid Cretaceous (Tmol/yr)')
pylab.ylabel('Decrease in seafloor weathering\nsince mid Cretaceous (Tmol/yr)')
pylab.show()
