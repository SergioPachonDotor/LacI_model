# from Model.Parameters import *
from os import register_at_fork
import numpy as np

class Parameters:
    def __init__(self) -> None:
        self.gamma = 0.0231
        self.gamma_r = 0.2
        self.x_0 = 0.12
        self.alpha = 1000
        self.k_R = 15.4
        self.R_0 = 0.04
        self.k_r = 2
        self.k_p = 0.00231
        self.n = 2
        self.q = 0
        self.p = 167
        self.landa = 1.5

class System(Parameters):

    def __init__(
                    self, 
                    beta=0, 
                    Extracellular_TMG=20,
                    LacI_Tetramer=0,
                    Active_LacI=0, 
                    Permease=0, 
                    mRNA=0, 
                    LacI_monomer=0,
                    Intracellular_TMG = 0,
                    Promoter_State=0,

                    ) -> None:
        super().__init__()
        self.beta = beta
        self.permease = Permease
        self.extracellular_tmg = Extracellular_TMG
        self.lacI_tetramer = LacI_Tetramer  
        self.active_lacI = Active_LacI
        self.mrna = mRNA
        self.lacI_monomers = LacI_monomer   
        self.intracellular_tmg = Intracellular_TMG  
        self.promoter_state = Promoter_State


    def set_beta(self, tmg):
        self.beta =  0.00123 * pow(tmg,0.6)

    def set_extracellular_tmg(self):
        self.extracellular_tmg = self.beta * self.permease

    def set_laci_tetramers(self):
        self.lacI_tetramer = self.lacI_monomer/100

    def set_laci_monomers(self):
        self.lacI_monomers = 1/(1 + (self.intracellular_tmg/self.x_0)**self.n) * self.lacI_tetramer
        
    def set_permease(self, r_y, y, dt):
        self.permease = (self.alpha * self.k_p * r_y - self.gamma * y)*dt 

    def set_permease_tau(self, tau):
        yss = ((self.alpha * self.k_p * self.mrna)/self.gamma)
        self.permease = yss + (self.permease - yss) * np.exp( -self.gamma*tau)

    def set_d_mrna(self, dt):
        self.mrna = self.k_r*dt

    def set_mrna_euler(self, dt):
        self.mrna (self.k_r/(1 + (self.active_lacI/self.R_0))- self.gamma_r * self.mrna) * dt

    def set_d_laci_monomers(self, dt):
        self.laci_monomers = (self.k_R - self.gamma * self.lacI_monomers)*dt

    def set_laci_monomers_euler(self):
        self.laci_monomers = self.k_R/ self.gamma

    def promoter_on(self):
        return self.landa*(self.R_0/self.lacI_tetramer)

    def promoter_off(self):
        return self.landa*(self.active_lacI/self.lacI_tetramer)

    def set_promoter_state(self):
        if self.promoter_state == 'on':
            pass
        elif self.promoter_state == 'off':
            pass
