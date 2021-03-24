from Parameters import *
from ModelData import *
from dbTools import *
from Model import *
from View import *
from numpy.random import rand
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

if __name__ == '__main__':

    if algorithm == 'poisson':
        for concentration in range(tmg_range):
            file = f'./data/simulations/cells_{cells}_tmg_{tmg}_state_{STATE}_poisson.csv'
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()
            p_sim = laci_poisson(TMG=tmg)
            main_view()

        