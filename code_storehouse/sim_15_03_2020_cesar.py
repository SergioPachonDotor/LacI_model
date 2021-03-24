from Parameters import *
from ModelData import *
from dbTools import *
from Model import *
from View import *
from numpy.random import rand
import numpy as np
import csv


def laci_poisson (TMG=5) -> None:
    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        for cell in range(cells): 
            state = STATE
            beta = 0.00123 * pow(TMG, 0.6)

            if state == 'on':
                y = 10000
                r_y = 10
            elif state == 'off':
                y = 0
                r_y = 0         
            R_mono = k_R/gamma
            R = 0.04                
            time_p = 0 
            promoter = 'off'   

            for step in range(steps):

                time_p += dt
                y += dy(r_y,y)
                R_mono += dR_mono(R_mono)
                x = extmg(beta, y)
                R_T = R_mono/100
                R = R_(x,R_T)
                rann = rand()
                rann_2 = rand()

                if promoter == 'off':
                    if rann < landa*(R_0/R_T)*dt:
                        if rand() < dr_y(): #modificar
                            r_y += 1
                        else:
                            r_y += 0
                        promoter = 'on'
                    else:
                        r_y += 0
                        promoter = 'off'

                elif promoter == 'on':
                    if rann_2 < landa*(R/R_T)*dt:
                        r_y += 0
                        promoter = 'off'
                    else:
                        if rand() < dr_y():
                            r_y += 1
                        else:
                            r_y += 0
                        promoter = 'on'

                if rand() < (k_R/10)*dt:
                    R_mono += 10

                R_mono -= gamma*R*dt
                
                if rand() < gamma_r*r_y*dt:
                            r_y -= 1   

                if y > -1:

                    data = Data(time=round(time_p,4), 
                                cell=round(cell,4), 
                                beta=round(beta,4), 
                                Extracellular_TMG= TMG,
                                LacI_Tetramer=round(R_T,4),
                                Active_LacI=round(R,4), 
                                Permease=round(y,4), 
                                mRNA=round(r_y,4), 
                                LacI_monomer=round(R_mono,4),
                                Intracellular_TMG = round(x,4),
                                Promoter_State=promoter)
                    model = data.data_model()
                    writer.writerow(model)
                    data = None

    f.close()

def laci_gillespie_1(TMG=5) -> None:

    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')

        for cell in range(cells): 

            #Initial Values
            tau_time = 0
            sampling_time= 0.1
            tref=0
            #tau_time_p = 0
            #_Species_Quantity_
            r_y = 0
            y = 0
            R = 0.04 

            beta = 0.00123 * pow(TMG, 0.6)
            R_mono = k_R/gamma
            promoter = 'off'  
               
            counter = 0
            reaction_number = 4
            a_j_p = [i for i in range(2)]
            a_j = [i for i in range(reaction_number)]

            while tau_time < tmax:
                x = extmg(beta, y)
                R_mono += dR_mono(R_mono)
                R_T = np.floor(R_mono/100)
                R = R_(x,R_T)

                #_Propensities_
                if promoter == 'off':
                    kry=0
                else:
                    kry=1
                a_j[0] = (k_r * kry)/2                    # Transcription rate (0.3/min)
                a_j[1] = gamma_r * r_y                   # mRNA Degradation Rate (2.5 min) 
                if promoter == 'off':
                    a_j[2] = landa * (R_0/R_T)
                else:
                    a_j[2] = landa * (R/R_T) 

                #a_total
                a_total = sum(a_j)
                
                #tau
                tau = (1/a_total)*np.log(1/rand())
                if tau_time+tau<tref:
                    tau_time += tau
                    dart = a_total * rand()
                    cms = np.cumsum(a_j)
                    nn=0
                    for m in cms:
                        if dart<m:
                            if nn == 0:
                                r_y+=1
                            elif nn == 1:
                                r_y-=1
                            elif nn == 2:
                                if promoter == 'off':
                                    promoter = 'on'
                                else:
                                    promoter = 'off' 
                            break
                        else: 
                            nn+=1
                        
                    tau_time+=tau
                
                    y+=(alpha * k_p*2 * r_y - gamma * y)*tau
                
                    if rand() < (k_R/50)*tau:
                        R_mono += 50

                    R_mono -= gamma*R_mono*tau
                else:
                    
                
                    y+=(alpha * k_p*2 * r_y - gamma * y)*(tref-tau_time)
                
                    if rand() < (k_R/50)*(tref-tau_time):
                        R_mono += 50

                    R_mono -= gamma*R_mono*(tref-tau_time)
                    
                    tau_time=tref
                    tref+=sampling_time
                    
                    data = Data(time=round(tau_time,4), 
                                cell=round(cell,4), 
                                beta=round(beta,4), 
                                Extracellular_TMG= TMG,
                                LacI_Tetramer=round(R_T,4),
                                Active_LacI=round(R,4), 
                                Permease=round(y,4), 
                                mRNA=round(r_y,4), 
                                LacI_monomer=round(R_mono,4),
                                Intracellular_TMG = round(x,4),
                                Promoter_State=promoter)

                    model = data.data_model()
                    writer.writerow(model)
                    data = None
                    

    f.close()


if __name__ == '__main__':

    if algorithm == 'poisson':
        for concentration in range(tmg_range):
            file = f'./data/simulations/cells_{cells}_tmg_{tmg}_state_{STATE}_poisson.csv'
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()
            p_sim = laci_poisson(TMG=tmg)
            main_view()

    elif algorithm == 'gillespie':
        for concentration in range(tmg_range):
            file = f'./data/simulations/cells_{cells}_tmg_{tmg}_state_{STATE}_gillespie.csv'
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()
            p_sim = laci_gillespie_1(TMG=tmg)
            main_view()