from Model.Parameters import *
import numpy as np

def beta_funct(tmg):
    return 0.00123 * pow(tmg,0.6)

def extmg(beta, y):
    return  beta * y

def R_T(R_mono):
    return R_mono/100

def R_(x, R_T):
    result = 1/(1 + (x/x_0)**n) * R_T
    return result

def dy(r_y, y, dt):
    result = (alpha * k_p * r_y - gamma * y)*dt 
    return result

def dy_tau(r_y, y, tau):
    yss = ((alpha * k_p * r_y)/gamma)
    return  yss + (y - yss) * np.exp(-gamma*tau)

def dr_y():
    return k_r*dt

def dr_y_euler(R,r_y, dt):
    return (k_r/(1 + (R/R_0))- gamma_r*r_y) * dt

def dR_mono(R_mono, dt):
    return (k_R - gamma * R_mono)*dt

def dR_mono_euler():
    return k_R/ gamma

def promoter_on(R_0,R_T):
    return landa*(R_0/R_T)

def promoter_off(R,R_T):
    return landa*(R/R_T)
