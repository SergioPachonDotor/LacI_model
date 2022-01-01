from DataTools.ModelData import *
from DataTools.dbTools import *
from Model.System import System as sys
from DataTools.View import *
from tqdm import tqdm
from numpy.random import rand

import warnings
import numpy as np
import csv

class Utilities(sys):
    def __init__(self) -> None:
        pass
    
    def define_state(self):
        if self.state == 'on':
            self.state = 'on'
        elif self.state == 'on':
            self.state = 'off'

class System:
    def __init__(
                    self, 
                    TMG: float, 
                    cells: int, 
                    file: str, 
                    tmax: int, 
                    steps: int,
                    state: str,
                    ) -> None:

        self.tmg: float = TMG
        self.cells: int = cells
        self.file: str = file
        self.dt: float = tmax/steps
        self.steps: int = steps
        self.tmax: int = tmax
        self.state: str = state

    def gillespie(self):
        SCHEMA = Schema.get_classic_schema()
        with open(self.file, mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter=',')
        
            for cell in tqdm(range(self.cells)):

                beta = sys.beta_funct(tmg=self.tmg)

                on_counter = 0
                on_time = 0
                off_counter = 0
                off_time = 0
                

                if self.state == 'off':
                    y = 0
                    r_y = 0
                    dna = 0
                    promoter = 'off'

                elif self.state == 'on':
                    y = 1000
                    r_y = 10
                    dna = 1
                    promoter = 'on'

                R_mono = 0
                R = 0.04
                R_T = 0
                dna = 1

                tautime = 0
                reference_time = 0

                a_j = [0 for i in range(4)]

                while tautime < tmax:

                    if R_T != 0:

                        a_j[0] = k_r * dna                  # Transcription  
                        a_j[1] = gamma_r * r_y              # mRNA Degradation
                        a_j[2] = landa * (R/R_T) 
                        a_j[3] = landa * (R_0/R_T)
                    
                    elif R_T == 0:

                        a_j[0] = k_r * dna                  # Transcription  
                        a_j[1] = gamma_r * r_y              # mRNA Degradation
                        a_j[2] = 0
                        a_j[3] = 0
                    
                    a_total = np.sum(a_j)

                    tau = (1/a_total)*np.log(1/rand())

                    dart = a_total * rand()
                    sum_a = 0
                    q = 0

                    if tautime + tau > reference_time: 
                        tau = reference_time - tautime
                        y += sys.dy(r_y, y, dt=tau)
                        R_mono += sys.dR_mono(R_mono, dt=tau)
                        x = sys.extmg(beta,y)
                        R_T = R_mono/100
                        R = sys.R_(x,R_T)

                        tautime = reference_time
                        data = Data(
                                time = round(tautime,4), 
                                cell = round(cell,4), 
                                beta = round(beta,4), 
                                Extracellular_TMG= self.tmg,
                                LacI_Tetramer = round(R_T,4),
                                Active_LacI = round(R,4), 
                                Permease = round(y,4), 
                                mRNA = round(r_y,4), 
                                LacI_monomer = round(R_mono,4),
                                Intracellular_TMG = round(x,4),
                                Promoter_State = dna,
                                On_Time = on_time,
                                Off_Time = off_time
                                )

                        model = data.data_model()
                        writer.writerow(model)
                        data = None
                        
                        reference_time += sampling_time
                    
                    else:
                        y += sys.dy(r_y, y, dt=tau)
                        R_mono += sys.dR_mono(R_mono, dt=tau)
                        x = sys.extmg(beta,y)
                        R_T = R_mono/100
                        R = sys.R_(x,R_T)

                        for i in range(len(a_j)):
                            sum_a += a_j[i]
                            if sum_a > dart:
                                q = i + 1
                                break
                        if q == 1:
                            # Transcription
                            r_y += 1

                        elif q == 2:
                            # mRNA degradation
                            r_y -= 1
                        
                        elif q == 3:
                            if dna > 0:
                                dna -= 1
                            elif dna <= 0:
                                dna = 0
                            promoter = 'off'

                        elif q == 4:
                            if dna == 0:
                                dna += 1
                                
                            elif dna > 0:
                                dna = 1
                            promoter == 'on'

                        tautime += tau

                    if self.state == 'off':
                        if y >= 250 and on_counter == 0:
                            on_time = round(tautime,4)
                            on_counter = 1

                    if self.state == 'on':
                        if y < 250 and off_counter == 0:
                            off_time = round(tautime,4)
                            off_counter = 1
        
        f.close()

