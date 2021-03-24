from Parameters import *

def beta():
    return 0.00123 * pow(4,0.6)

def beta():
    return 0.00123 * pow(4,0.6)

def extmg(beta, y):
    """TMG Concentration"""
    return  beta * y

def R_T(R_mono):
    return R_mono/100

def R_(x, R_T):
    return 1/(1 + pow((x/x_0),n)) * R_T

def dy(r_y, y):
    return (alpha * k_p * r_y - gamma * y)*dt 

def dr_y():
    return k_r*dt

def dR_mono(R_mono):
    return (k_R - gamma * R_mono)*dt

def promoter_on(R_0,R_T):
    return landa*(R_0/R_T)

def promoter_off(R,R_T):
    return landa*(R/R_T)
