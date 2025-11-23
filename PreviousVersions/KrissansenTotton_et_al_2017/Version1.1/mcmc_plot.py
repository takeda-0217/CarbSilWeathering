import emcee
import numpy
import matplotlib.pyplot as pl
from corner import corner
import corner
from model_functions import forward_model
#from foward_model_for_mcmc import forward_model_old
from multiprocessing import Pool
print("BACKEND:", pl.get_backend())

nwalk  = 1000
nsteps = 10000

chain  = numpy.load(r'/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/my_output/inverse_model/full_run/newchain2.npy') 
lnprob = numpy.load(r'/home/taketomo0217/my_study/CarbSilWeathering/PreviousVersions/KrissansenTotton_et_al_2017/Version1.1/my_output/inverse_model/full_run/newln2.npy') 

##### The rest is all plotting

# Plot a few individual chains
fig, ax = pl.subplots(4) 
for n in range(nwalk):
  ax[0].plot(chain[n,:,0])  # CO2 dependence alpha
  ax[1].plot(chain[n,:,1])  # e-folding temperature of cont. weathering Te
  ax[2].plot(chain[n,:,2])  # Weatherability factor W
  ax[3].plot(chain[n,:,3])  # Climate sensitivity parameter delata_T2x

# find highest likelihood run
logprob=numpy.array(lnprob)
values=chain
ii,jj = numpy.unravel_index(logprob.argmax(), logprob.shape)
print ("indeces for best",ii,jj)
print ("loglikelihood and values",logprob[ii,jj],values[ii,jj,:])

# Plot the corner plot, discarding the first 1000 steps as burn-in
production = chain[:,0:,:]
# production = chain[:,0:,:] ## Use this if you have <1000 steps

s         = production.shape                      # 3D array (nwalkers, nsteps_after_burn_in, ndim)
flatchain = production.reshape(s[0] * s[1], s[2]) # convert to 2D array for corner plot

flatchain2=numpy.copy(flatchain)
## convert some variables:
flatchain2[:,4]=flatchain2[:,4]+1           # make outgassing relative Cretaceous outgassing, V+1
## Weatherability is already Cretaceous weatherability, W+1, assuming w=1+W with plus sign
flatchain2[:,5]=flatchain2[:,5]+1           # Carbonate weathering modifier
flatchain2[:,8]=flatchain2[:,8]/1000.0      # Convert to ky: Mixing time for pore-space
flatchain2[:,13]=flatchain2[:,13]/1000.0    # Convert to kJ/mol: activation energy seafloor weathering
flatchain2[:,6]=flatchain2[:,6]/1e12        # Convert to Tmol: modern outgassing
flatchain2[:,7]=flatchain2[:,7]/1e12        # Convert to Tmol: modern carbonate weathering

from matplotlib import rc
## Plot posteriors as corner plots (compare to Fig. 6 in the paper)
# corner.corner(flatchain2[:,[0,1,2,3,4]],     quantiles=[0.16, 0.5, 0.84],labels = [r"CO$_2$-dependence, $\alpha$", "Temp. dep. cont.\nweathering, $T_e$ (K)", "Relative Cretaceous\nweatherability, 1+$W$","Climate sensitivity,\n${\Delta}T_{2X}$ (K)","Relative Cretaceous\noutgassing, 1+$V$"])#,truths=values[ii,jj,:])
# corner.corner(flatchain2[:,[5,8,11,13,6,16]], quantiles=[0.16, 0.5, 0.84],labels = ["Carbonate weath.\nmodifier, 1+$C_{WF}$",r"Circulation time, $\tau$ (kyr)","Surface-deep\ntemp. gradient, $a_{grad}$","Temp. dependence\nseafloor, $E_{bas}$ (kJ/mol)","Modern outgassing,\n$F_{out}^{mod}$ (Tmol C/yr)","Paleogeography parameter,\n${\Delta}P$ (K)"])#,truths=values[ii,jj,:])
# corner.corner(flatchain2[:,[7,9,10,12,14,15]], quantiles=[0.16, 0.5, 0.84],labels = ["Modern carb.\nweathering, $F_{carb}^{mod}$ (Tmol C/yr)","Carb. precip.\ncoefficient, $n$","Modern seafloor diss.\nrelative precip.","pH dependence\nseafloor, $\gamma$","Modern pelagic\nfraction",r"Spreading rate dep., $\beta$"])#,truths=values[ii,jj,:])

# pl.figure()
# pl.plot([0, 1], [0, 1])
# pl.show()

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
    if numpy.any(~numpy.isfinite(outputs)):
        continue
    sp=((x_ex[6]+x_ex[4]*x_ex[6])/x_ex[6])**x_ex[15]  # sp = (1+V)**beta: spreading rate relative to modern
    mc_plotter_spread(outputs,"y",legend_counter,sp)
    legend_counter=legend_counter+1
    print("RUN", legend_counter)

### This is important. This takes 1000 sets of parameter values from your posterior
### and re-runs the forward model 1000 times to get distributions for the time-evolution
### of different model parameters e.g. Fig. 5a-f
### Another script is called to actually do the plotting
mega_output=[]
spread_output=[]
carbw_factor=[]
sil_change=[]
seafloor_change=[]
for x_ex in flatchain[numpy.random.randint(len(flatchain), size=100)]: # default size=1000
    #print (x_ex)
    #outputs=forward_model_old(x_ex[0],x_ex[1],x_ex[2],x_ex[3],x_ex[4],x_ex[5],x_ex[6],x_ex[7],x_ex[8],x_ex[9],x_ex[10],x_ex[11],x_ex[12],x_ex[13],x_ex[14],x_ex[15],x_ex[16])
    [outputs,imbalance]= forward_model(x_ex[8],x_ex[6],x_ex[9],x_ex[0],x_ex[1],0.45e12,x_ex[10],0.01,x_ex[2],x_ex[3],x_ex[4],x_ex[7],x_ex[14],x_ex[5],x_ex[11],x_ex[12],x_ex[15],x_ex[13],x_ex[16])
    if numpy.any(~numpy.isfinite(outputs)):
        continue
    spread_output.append( ((x_ex[6]+x_ex[4]*x_ex[6])/x_ex[6])**x_ex[15])
    mega_output.append(outputs)
    carbw_factor.append((1+x_ex[5])*outputs[20][99]/outputs[20][0])
    sil_change.append(outputs[20][99]-outputs[20][0])
    seafloor_change.append(outputs[19][99]-outputs[19][0])
    print("RUN", legend_counter)
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
pylab.hist2d(sil_change, seafloor_change, range=[[-1, 6], [0, 6]],bins=30,density=True,cmap=pylab.cm.jet)
pylab.colorbar(label='Probability density')
pylab.xlabel('Decrease in continental weathering\nsince mid Cretaceous (Tmol/yr)')
pylab.ylabel('Decrease in seafloor weathering\nsince mid Cretaceous (Tmol/yr)')

# carbonate plot:
pylab.subplot(3, 3, 7)
pylab.hist(numpy.array(carbw_factor),bins=30,color='grey',density=True)
pylab.xlabel('Relative change carbonate weathering')
pylab.ylabel('Probability density')
print (numpy.percentile(numpy.array(carbw_factor),2.5),numpy.percentile(numpy.array(carbw_factor),50),numpy.percentile(numpy.array(carbw_factor),97.5))

pylab.show()


pylab.figure()
pylab.hist(sil_change,bins=30,color='grey',density=True)
pylab.xlabel('Absolute change silicate weathering')
pylab.ylabel('Probability density')

pylab.figure()
pylab.hist(seafloor_change,bins=30,color='grey',density=True)
pylab.xlabel('Absolute change seafloor weathering')
pylab.ylabel('Probability density')

pylab.figure()
pylab.hist2d(sil_change, seafloor_change, range=[[-1, 6], [-1, 6]],bins=30,density=True,cmap=pylab.cm.jet)
pylab.colorbar(label='Probability density')
pylab.xlabel('Decrease in continental silicate weathering\nsince mid Cretaceous (Tmol/yr)')
pylab.ylabel('Decrease in seafloor weathering\nsince mid Cretaceous (Tmol/yr)')
pylab.show()