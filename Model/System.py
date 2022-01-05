# from Model.Parameters import *
from DataTools.ModelData import Data
import warnings
import numpy as np
import copy

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
                    Active_LacI=0.04, 
                    Permease=0, 
                    mRNA=0, 
                    LacI_monomer=0,
                    Intracellular_TMG = 0,
                    Promoter = 0,
                    On_Time = 0,
                    Off_Time = 0,
                    Size = 0,
                    DNA = 1,
                    propensities = 4
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
        self.dna = DNA
        self.promoter = Promoter
        self.on_time = On_Time
        self.off_time = Off_Time
        self.size = Size
        self.aj = [0 for i in range(propensities)]
        self.on_counter = 0
        self.off_counter = 0
        
    def set_beta(self, tmg):
        self.beta =  0.00123 * pow(tmg,0.6)

    def set_intracellular_tmg(self):
        self.intracellular_tmg = self.beta * self.permease

    def set_laci_tetramers(self):
        self.lacI_tetramer = self.lacI_monomers/100

    def set_active_laci(self):
        self.active_lacI = 1/(1 + (self.intracellular_tmg/self.x_0)**self.n) * self.lacI_tetramer
        
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
        self.lacI_monomers += (self.k_R - self.gamma * self.lacI_monomers)*dt

    def set_laci_monomers_euler(self):
        self.laci_monomers = self.k_R/ self.gamma

    def promoter_on(self):
        return self.landa*(self.R_0/self.lacI_tetramer)

    def promoter_off(self):
        return self.landa*(self.active_lacI/self.lacI_tetramer)

    def set_promoter_state(self):

        if self.promoter == 'on':
            self.permease = 1000
            self.mrna = 10

        elif self.promoter == 'off':
            self.permease = 0
            self.mrna = 0

    def model_to_save(self, tautime, cell, tmg):
        data = Data(
                    time = round(tautime,4), 
                    cell = round(cell,4), 
                    beta = round(self.beta,4), 
                    Extracellular_TMG= tmg,
                    LacI_Tetramer = round(self.lacI_tetramer,4),
                    Active_LacI = round(self.active_lacI,4), 
                    Permease = round(self.permease,4), 
                    mRNA = round(self.mrna,4), 
                    LacI_monomer = round(self.lacI_monomers,4),
                    Intracellular_TMG = round(self.intracellular_tmg,4),
                    Promoter_State = self.promoter,
                    On_Time = self.on_time,
                    Off_Time = self.off_time
                    )
        return data.data_model()

    def solve_eqs(self, tau):
        self.set_permease_tau(tau=tau)
        self.set_d_laci_monomers(dt=tau)
        self.set_intracellular_tmg()
        self.set_laci_tetramers()
        self.set_active_laci()

    def set_aj(self, reaction_number=4):
        self.aj = [0 for i in range(reaction_number)]

    def set_reactions_classic_gillespie(self):
        if self.lacI_tetramer > 0:

            self.aj[0] = self.k_r * self.dna                  # Transcription  
            self.aj[1] = self.gamma_r * self.mrna              # mRNA Degradation
            self.aj[2] = self.landa * (self.active_lacI/self.lacI_tetramer) 
            self.aj[3] = self.landa * (self.R_0/self.lacI_tetramer)
        
        elif self.lacI_tetramer <= 0:

            self.aj[0] = self.k_r * self.dna                  # Transcription  
            self.aj[1] = self.gamma_r * self.mrna              # mRNA Degradation
            self.aj[2] = 0
            self.aj[3] = 0

    def set_q_tau(self):
        a_total = np.sum(self.aj)
        tau = (1/a_total)*np.log(1/np.random.rand())
        dart = a_total * np.random.rand()
        sum_a = 0
        q = 0
        for i in range(len(self.aj)):
            sum_a += self.aj[i]
            if sum_a > dart:
                q = i + 1
                break
        return q, tau

    def react(self, q):
        if q == 1:
            # Transcription
            self.mrna += 1

        elif q == 2:
            # mRNA degradation
            self.mrna -= 1
        
        elif q == 3:
            if self.dna > 0:
                self.dna -= 1
            elif self.dna <= 0:
                self.dna = 0
            self.promoter = 'off'

        elif q == 4:
            if self.dna == 0:
                self.dna += 1
                
            elif self.dna > 0:
                self.dna = 1
            self.promoter = 'on'

    def switch_counter(self, tautime):
            if self.promoter == 'off':
                if self.permease >= 250 and self.on_counter == 0:
                    self.on_time = round(tautime,4)
                    self.on_counter = 1

            if self.promoter == 'on':
                if self.permease < 250 and self.off_counter == 0:
                    self.off_time = round(tautime,4)
                    self.off_counter = 1

    def get_single_switch(self, tautime=0, cell=1, tmg=0):
        """An efficient version of the switch counter"""

        if self.promoter == 'off':
            if tmg == 0:
                print(f'Unable to simulate\nTMG: {tmg}')
                exit()

            if self.permease >= 100 and self.on_counter == 0:

                self.on_time = round(tautime,4)
                data = Data(
                            cell = cell,
                            On_Time = self.on_time,
                            Off_Time= self.off_time,
                            Promoter_State= self.promoter
                            )  
                model = data.data_time_distribution_model()
                self.on_counter = 1
                return model

        if self.promoter == 'on':

            if self.permease < 100 and self.off_counter == 0:

                self.off_time = round(tautime,4)
                data = Data(
                            cell = cell,
                            On_Time = self.on_time,
                            Off_Time= self.off_time,
                            Promoter_State= self.promoter
                            )  
                model = data.data_time_distribution_model()
                self.off_counter = 1
                return model

    def get_double_switch(self, tautime=0, cell=1, tmg=0):

        if tmg < 5:
            print(f'Unable to simulate\nTMG: {tmg}')
            exit()

        if self.permease >= 100 and self.on_counter == 0:

            self.on_time = round(tautime,4)
            data = Data(
                        cell = cell,
                        On_Time = self.on_time,
                        Off_Time= self.off_time,
                        Promoter_State= self.promoter
                        )  
            model = data.data_time_distribution_model()

            self.on_counter = 1
            self.promoter = 'on'
            return model

        if self.permease < 100 and self.off_counter == 0:

            self.off_time = round(tautime,4)
            data = Data(
                        cell = cell,
                        On_Time = self.on_time,
                        Off_Time= self.off_time,
                        Promoter_State= self.promoter
                        )  
            model = data.data_time_distribution_model()

            self.off_counter = 1
            self.promoter = 'off'
            return model


class DivisionSystem(System, Parameters):

    def __init__(
                    self, 
                    birth_size: float, 
                    division_time: float,
                    Extracellular_TMG: float,
                    Promoter: str):
        super().__init__()

        self.size = birth_size
        self.division_size = division_time
        self.mu = np.log(2)/division_time
        self.reference_division_time = (1/self.mu) * np.log(((birth_size + np.random.normal(loc=1.0, scale=0.05))/birth_size))
        self.extracellular_tmg = Extracellular_TMG
        self.promoter = Promoter
               
    def calculate_tau(self, propensity):
        try:
            return -(1/propensity) * np.log(np.random.rand())
        except ZeroDivisionError: 
            return 100

    def calculate_division_tau(self, propensity):

        try:

            warnings.simplefilter("ignore")
            p = (1/self.mu) * np.log(1 - (self.mu * np.log(np.random.rand())/(propensity * self.size)))
            return p

        except ZeroDivisionError:
            return 100

    def set_calculated_sorted_tau(self):
        calculated_tau = [self.calculate_tau(propensity=i) for i in self.aj]
        calculated_tau[0] = self.calculate_division_tau(propensity=self.aj[0])
        calculated_tau.append(self.reference_division_time)
        if self.mrna <= 0:
            calculated_tau[1] = np.inf
        tau = np.min(calculated_tau)
        q = np.argmin(calculated_tau)

        return q + 1, tau

    def division_react(self, q):
        if q == 1:
            # Transcription
            self.mrna += 1

        elif q == 2:
            # mRNA degradation
            if self.mrna <= 0:
                self.mrna = 0

            elif self.mrna > 0:
                self.mrna -= 1
        
        elif q == 3:
            if self.dna > 0:
                self.dna -= 1

            elif self.dna <= 0:
                self.dna = 0

            self.promoter = 'off'

        elif q == 4:
            if self.dna == 0:
                self.dna += 1
                
            elif self.dna > 0:
                self.dna = 1

            self.promoter = 'on'

        elif q == 5:
            # Set all division conditions like dilution of species
            self.size /= 2
            self.reference_division_time = (1/self.mu) * np.log((self.size + np.random.normal(loc=1.0, scale=0.05))/self.size)   
            
            # dilute_species
            self.mrna = np.random.binomial(self.mrna, 0.5)
            self.permease = np.random.binomial(self.permease, 0.5)
            self.lacI_monomers =np.random.binomial(self.lacI_monomers, 0.5)
            self.active_lacI =np.random.binomial(self.active_lacI, 0.5)
            self.lacI_tetramer =np.random.binomial(self.lacI_tetramer, 0.5)

    def update_tau(self, tau):
            self.size *= np.exp(self.mu * tau)

    def switch_counter(self, tautime):
        ## Change it to detect size over permease concentration
        if self.promoter == 'off':
            if self.permease >= 250/2 and self.on_counter == 0:
                self.on_time = round(tautime,4)
                self.on_counter = 1

        if self.promoter == 'on':
            if self.permease < 250/2 and self.off_counter == 0:
                self.off_time = round(tautime,4)
                self.off_counter = 1

    def division_model_to_save(self, tautime, cell, tmg):
        data = Data(
                    time = round(tautime,4), 
                    cell = round(cell,4), 
                    beta = round(self.beta,4), 
                    Extracellular_TMG= tmg,
                    LacI_Tetramer = round(self.lacI_tetramer,4),
                    Active_LacI = round(self.active_lacI,4), 
                    Permease = round(self.permease,4), 
                    mRNA = round(self.mrna,4), 
                    LacI_monomer = round(self.lacI_monomers,4),
                    Intracellular_TMG = round(self.intracellular_tmg,4),
                    Promoter_State = self.promoter,
                    On_Time = self.on_time,
                    Off_Time = self.off_time,
                    size = round(self.size, 4)
                    )
        return data.data_division_model()
