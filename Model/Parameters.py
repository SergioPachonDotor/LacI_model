import numpy as np

#_____Time_____
tmax = 75 #minutes
dt = 0.1
steps = int(tmax/dt)

#_____Cell_Number_____
cells = 1000

#_____TMG_____
tmg = 10
tmg_range = 9

#_____STATE_____
# Choose Algorithm Initial State
STATE = 'off'
P_on = 1
P_off = 0

#_____Algorithm_____
#   You can choose between: 
#     - Euler: euler
#     - Poisson: poisson, poisson_delta_tmg
#     - Gillespie: gillespie, gillespie_full_noise, gillespie_delta_tmg, gillespie_on_off, gillespie_on_off_double, simulate_adder

algorithm = 'simulate_adder'

#_____Parameters_____
gamma = 0.0231
gamma_r = 0.2
x_0 = 0.12
alpha = 1000
k_R = 15.4
R_0 = 0.04
k_r = 2
k_p = 0.00231
n = 2
q = 0
p = 167
landa = 1.5

#_____Extra_____

sampling_time = 0.1
c_assoc = 0.00001               # Association Constant. 0.00001 Associations/s
c_diss = 0.005                   # Complex Degradation Constant 0.005/s
