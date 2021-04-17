from DataTools.ModelData import *
from DataTools.dbTools import *
from Model.System import *
from Model.Parameters import *
from DataTools.View import *
from tqdm import tqdm
from numpy.random import rand
import numpy as np
import csv


init_state = STATE

def gillespie(TMG,Cells):
    
    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        
        for cell in tqdm(range(Cells)):

            beta = beta_funct(tmg=TMG)
            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            dna = 1

            if init_state == 'off':
                y = 0
                r_y = 0

            elif init_state == 'on':
                y = 1000
                r_y = 10

            R_mono = 0
            R = 0.04
            tautime = 0
            reference_time = 0

            a_j = [i for i in range(2)]

            while tautime < tmax:

                a_j[0] = k_r * dna          
                a_j[1] = gamma_r * r_y

                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                if tautime + tau > reference_time: 
                    tau = reference_time - tautime

                    y += dy(r_y, y, dt=tau)
                    R_mono += dR_mono(R_mono, dt=tau)
                    x = extmg(beta,y)
                    R_T = R_mono/100
                    R = R_(x,R_T)

                    tautime = reference_time
                    data = Data(
                            time=round(tautime,4), 
                            cell=round(cell,4), 
                            beta=round(beta,4), 
                            Extracellular_TMG= TMG,
                            LacI_Tetramer=round(R_T,4),
                            Active_LacI=round(R,4), 
                            Permease=round(y,4), 
                            mRNA=round(r_y,4), 
                            LacI_monomer=round(R_mono,4),
                            Intracellular_TMG = round(x,4),
                            Promoter_State='None',
                            On_Time=on_time,
                            Off_Time=off_time
                            )

                    model = data.data_model()
                    writer.writerow(model)
                    data = None

                    reference_time += sampling_time
                
                else:
                    y += dy(r_y, y, dt=tau)
                    R_mono += dR_mono(R_mono, dt=tau)
                    x = extmg(beta,y)
                    R_T = R_mono/100
                    R = R_(x,R_T)
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
                    tautime += tau

                if init_state == 'off':
                    if y >= 500 and on_counter == 0:
                        on_time = round(tautime,4)
                        on_counter = 1

                if init_state == 'on':
                    if y < 500 and off_counter == 0:
                        off_time = round(tautime,4)
                        off_counter = 1


    f.close()


def gillespie_full_noise(TMG,Cells):
    
    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        
        for cell in tqdm(range(Cells)):

            beta = beta_funct(tmg=TMG)

            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            

            if init_state == 'off':
                y = 0
                r_y = 0
                dna = 0
                promoter = 'off'

            elif init_state == 'on':
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

            a_j = [0 for i in range(6)]

            while tautime < tmax:

                if R_T != 0:

                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = alpha * k_p * r_y          # Protein Synthesis
                    a_j[3] = gamma * y                  # Protein Degradation
                    a_j[4] = landa * (R/R_T) 
                    a_j[5] = landa * (R_0/R_T)
                
                elif R_T == 0:

                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = alpha * k_p * r_y          # Protein Synthesis
                    a_j[3] = gamma * y                  # Protein Degradation
                    a_j[4] = 0
                    a_j[5] = 0
                
                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                if tautime + tau > reference_time: 
                    tau = reference_time - tautime

                    R_mono += dR_mono(R_mono, dt=tau)
                    x = extmg(beta,y)
                    R_T = R_mono/100
                    R = R_(x,R_T)

                    tautime = reference_time
                    data = Data(
                            time = round(tautime,4), 
                            cell = round(cell,4), 
                            beta = round(beta,4), 
                            Extracellular_TMG= TMG,
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

                    R_mono += dR_mono(R_mono, dt=tau)
                    x = extmg(beta,y)
                    R_T = R_mono/100
                    R = R_(x,R_T)
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
                        y += 1
                    
                    elif q == 4:
                        y -= 1
                    
                    elif q == 5:
                        if dna > 0:
                            dna -= 1
                        elif dna <= 0:
                            dna = 0
                        promoter = 'off'

                    elif q == 6:
                        if dna == 0:
                            dna += 1
                            
                        elif dna > 0:
                            dna = 1
                        promoter == 'on'

                    tautime += tau

                if init_state == 'off':
                    if y >= 500 and on_counter == 0:
                        on_time = round(tautime,4)
                        on_counter = 1

                if init_state == 'on':
                    if y < 500 and off_counter == 0:
                        off_time = round(tautime,4)
                        off_counter = 1
    f.close()


if __name__ == '__main__':

    if algorithm == 'gillespie':

        file = f'./simulation_data/gillespie/gillespie_{init_state}_{tmg}_cells_{cells}.csv'
        Tools(file=file).delete()
        Tools(name=file, schema=SCHEMA).create()
        gillespie(TMG=tmg,Cells=cells)
        main_view()
    
    elif algorithm == 'gillespie_full_noise':

        file = f'./simulation_data/gillespie_full_noise/gillespie_full_noise_{init_state}_{tmg}_cells_{cells}.csv'
        Tools(file=file).delete()
        Tools(name=file, schema=SCHEMA).create()
        gillespie_full_noise(TMG=tmg,Cells=cells)
        main_view()

    elif algorithm == 'gillespie_delta_tmg':
        
        for tmg_concentration in tqdm(range(tmg_range)):
            file = f'./simulation_data/gillespie_delta_tmg/gillespie_delta_tmg_{init_state}_{tmg}_cells_{cells}.csv'
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()
            gillespie(TMG=tmg,Cells=cells)
        main_view()