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
            tau_time_p = 0
            #___Species_Quantity___
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

            while counter < steps:
                t_counter = 0
                x = extmg(beta, y)
                R_mono += dR_mono(R_mono)
                R_T = np.floor(R_mono/100)
                R = R_(x,R_T)

                #___Propensities___
                a_j[0] = 0.3                            # Transcription rate (0.3/min)
                a_j[1] = 10.02 * r_y                    # Translation rate (10.02/min)
                a_j[2] = 0.120 * r_y                    # mRNA Degradation Rate (2.5 min)
                a_j[3] = gamma * y                      # Protein Degradation Rate (60 min)
                
                #___Promotor_Propensities___
                a_j_p[0] = landa * (R_0/R_T)            # Promoter On            
                a_j_p[1] = landa * (R/R_T)              # Promoter Off

                #a_total
                a_total = sum(a_j)
                a_total_p = sum(a_j_p)
                
                #tau

                tau_p = (1/a_total_p)*np.log(1/rand())
                tau_time_p += tau_p

                #q
                dart = a_total * rand()
                sum_a = 0
                q = 0

                #q promoter
                dart_p = a_total_p * rand()
                sum_a_p = 0
                q_p = 0

                for i in range(2):
                    sum_a_p = sum_a_p + a_j_p[i]
                    if sum_a_p > dart_p:
                        q_p = i
                        break
                if q_p == 0:
                    promoter = 'on'
                    #tau
                    tau = (1/a_total)*np.log(1/rand())
                    tau_time += tau
                    #q
                    dart = a_total * rand()
                    sum_a = 0
                    q = 0

                    for i in range(reaction_number):
                        sum_a = sum_a + a_j[i]
                        if sum_a > dart:
                            q = i
                            break
                    if q == 0 :
                        r_y += 1

                    elif q == 1:
                        y += 1

                    elif q == 2:
                        if r_y != 0:
                            r_y -= 1

                    elif q == 3:
                        if y != 0:
                            y -= 1
                
                elif q_p == 1:
                    promoter = 'off'
                    #tau
                    tau = (1/a_total)*np.log(1/rand())
                    tau_time += tau
                    #q
                    dart = a_total * rand()
                    sum_a = 0
                    q = 0
                    for i in range(reaction_number):
                        sum_a = sum_a + a_j[i]
                        if sum_a > dart:
                            q = i
                            break
                    if q == 0 :
                        r_y += 0

                    elif q == 1:
                        if r_y != 0:
                            y += 1

                    elif q == 2:
                        if r_y != 0:
                            r_y -= 1

                    elif q == 3:
                        if y != 0:
                            y -= 1

                if rand() < (k_R/10)*dt:
                    R_mono += 10

                R_mono -= gamma*R*dt

                counter += 1 
                t_counter += tau

                if t_counter > -0.01:

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
        

