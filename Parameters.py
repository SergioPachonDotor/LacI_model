import numpy as np 

#_____Time_____
tmax = 1440 #minutes
dt = 0.01
steps = int(tmax/dt)

#_____Cell_Number_____
cells = 1

#_____TMG_____
tmg = 40
tmg_range = 1

#_____STATE_____
STATE = 'off'
P_on = 1
P_off = 0

#_____Algorithm_____
algorithm = 'gillespie'

#_____Parameters_____
gamma = 0.0231
gamma_r = 0.1
x_0 = 0.12
alpha = 1000
k_R = 15.4
R_0 = 0.04
k_r = 1
k_p = 0.00231
n = 2
p = 167
landa = 1.5
#landa = 0.015
