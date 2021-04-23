import numpy as np
from PyEcoLib.PopSimulator import PopSimulator
meanbirthsize = 1 # micron '  8      
doubling_time = 18 #min 
gr = np.log(2)/doubling_time
tmax = 7*doubling_time #min 
sample_time = 0.1*doubling_time #min 
div_steps = 10
ncells = 20
CV2sz=0.01
v0=meanbirthsize*np.random.gamma(shape=1/CV2sz,scale=CV2sz,size=ncells)
sim = PopSimulator(ncells=ncells,gr = gr, sb=meanbirthsize, steps = div_steps,nu=1,V0array=v0,CV2div=0.005,CV2gr=0.01) #Initializing the simulator
sim.szdyn(tmax = tmax, sample_time = 0.1*doubling_time, FileName=  "./data2Pop.csv", DivEventsFile="./DivEvents2.csv")
