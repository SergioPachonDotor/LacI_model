from DataTools.ModelData import *
from DataTools.dbTools import *
from Model.System import *
from Model.Parameters import *
from DataTools.View import *
import csv

init_state = STATE

def euler(TMG):
    
    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')

        for cell in range(cells):
            beta = beta_funct(tmg=TMG)
            time_ = 0

            if init_state == 'off':
                y = 0
                r_y = 0

            elif init_state == 'on':
                y = 1000
                r_y = 10

            R_mono = 0
            R = 0.04
            counter = 0

            for step in range(steps):

                time_ += dt
                counter += dt
                y += dy(r_y, y)
                r_y += dr_y_euler(R,r_y)
                R_mono = dR_mono_euler()
                x = extmg(beta,y)
                R_T = R_mono/100
                R = R_(x,R_T)

                if step == steps - 1:   
                    data = Data(time=round(time_,4), 
                            cell=round(cell,4), 
                            beta=round(beta,4), 
                            Extracellular_TMG= TMG,
                            LacI_Tetramer=round(R_T,4),
                            Active_LacI=round(R,4), 
                            Permease=round(y,4), 
                            mRNA=round(r_y,4), 
                            LacI_monomer=round(R_mono,4),
                            Intracellular_TMG = round(x,4),
                            Promoter_State='None')
                    model = data.data_model()
                    writer.writerow(model)
                    data = None
    f.close()


if __name__ == '__main__':
    file = f'./simulation_data/euler/deterministic_steady_state_euler_{init_state}.csv'
    Tools(file=file).delete()
    Tools(name=file, schema=SCHEMA).create()
    for tmg in range(120):
        euler(TMG=tmg/2)