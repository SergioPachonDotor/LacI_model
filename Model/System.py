from Model.Parameters import *

def beta_funct(tmg):
    return 0.00123 * pow(tmg,0.6)

def extmg(beta, y):
    return  beta * y

def R_T(R_mono):
    return R_mono/100

def R_(x, R_T):
    return 1/(1 + pow((x/x_0),n)) * R_T

def dy(r_y, y, dt):
    return (alpha * k_p * r_y - gamma * y)*dt 

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
