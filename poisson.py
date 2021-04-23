from DataTools.ModelData import *
from DataTools.dbTools import *
from Model.System import *
from Model.Parameters import *
from DataTools.View import *
from numpy.random import rand
from tqdm import tqdm
import csv



def laci_poisson (TMG=5, Cells=cells) -> None:
    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        for cell in tqdm(range(Cells)): 
            state = STATE
            beta = beta_funct(tmg=TMG)

            if state == 'on':
                y = 1000
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
                y += dy(r_y,y, dt=dt)
                R_mono += dR_mono(R_mono,dt=dt)
                x = extmg(beta, y)
                R_T = R_mono/100
                R = R_(x,R_T)

                if promoter == 'off':
                    if rand() < landa*(R_0/R_T)*dt:
                        if rand() < dr_y(): #modificar
                            r_y += 1
                        else:
                            r_y += 0
                        promoter = 'on'
                    else:
                        r_y += 0
                        promoter = 'off'

                elif promoter == 'on':
                    if rand() < landa*(R/R_T)*dt:
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

                    data = Data(
                                time = round(time_p,4), 
                                cell = round(cell,4), 
                                beta = round(beta,4), 
                                Extracellular_TMG = TMG,
                                LacI_Tetramer = round(R_T,4),
                                Active_LacI = round(R,4), 
                                Permease = round(y,4), 
                                mRNA = round(r_y,4), 
                                LacI_monomer = round(R_mono,4),
                                Intracellular_TMG = round(x,4),
                                Promoter_State = promoter
                                )
                    model = data.data_model()
                    writer.writerow(model)
                    data = None
    f.close()

if __name__ == '__main__':

    if algorithm == 'poisson':

        file = f'./simulation_data/poisson/cells_{cells}_tmg_{tmg}_state_{STATE}_poisson.csv'
        Tools(file=file).delete()
        Tools(name=file, schema=SCHEMA).create()
        p_sim = laci_poisson(TMG=tmg, Cells=cells)
        main_view()

    elif algorithm == 'poisson_delta_tmg':

        for tmg_concentration in tqdm(range(tmg_range)):
            file = f'./simulation_data/poisson_delta_tmg/cells_{cells}_tmg_{tmg}_state_{STATE}_poisson_delta_tmg.csv'
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()
            p_sim = laci_poisson(TMG=tmg, Cells=cells)
            main_view()

        