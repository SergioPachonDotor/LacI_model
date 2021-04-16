
import csv
import numpy as np
import os

init_state = 'on'
file = f'deterministic_steady_state_euler_{init_state}.csv'

SCHEMA = ['Time', 'Cell', 'beta', 'x', 'y', 'r_y', 'R_mono', 'R_T', 'R', 'T_']

class Data:
    
    def __init__(self, 
                 time=float, 
                 cell=int, 
                 beta=float, 
                 x=float,
                 R_T=float,
                 R=float, 
                 y=float, 
                 r_y=float, 
                 R_mono=float,
                 T = float):

        #Time (dt)| Sample (trayectory id) | Variables
        self.time = time
        self.cell = cell
        self.beta = beta
        self.x = x
        self.R_T = R_T
        self.R = R
        self.y = y
        self.r_y = r_y
        self.R_mono = R_mono
        self.T = T
    
    def data_model(self):
        model = {'Time':f'{self.time}', 
                'Cell':f'{self.cell}', 
                'beta':f'{self.beta}', 
                'x':f'{self.x}',
                'y':f'{self.y}',
                'r_y':f'{self.r_y}',
                'R_mono':f'{self.R_mono}',
                'R_T':f'{self.R_T}',
                'R':f'{self.R}',
                'T_':f'{self.T}'
                }
        return model

class CRUD:

    def __init__(self, data=dict, file=str, name='simulation_data.csv', schema=list):
        self.data = data
        self.file = file
        self.name = name
        self.schema = schema
    
    def create(self):
        try:
            with open(self.name, mode='x') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.schema, delimiter='|')
                writer.writeheader()
                csvfile.close()
        except:
            print(f'File "{self.name}" already exist')
            pass

    def update(self):      
        with open(self.file, mode='a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
            writer.writerow(self.data)
            f.close()

    def delete(self):
        try:
            os.remove(self.file)
        except:
            pass

class Simulation: 

    def __init__(self, steps, cells):
        self.tmax = 11520 #minutes
        self.steps = steps
        self.dt = self.tmax/self.steps
        self.cells = cells
        #_____Parameters_____
        self.gamma = 0.0231
        self.gamma_r = 0.1
        self.x_0 = 0.12
        self.alpha = 1000
        self.k_R = 15.4
        self.R_0 = 0.04
        self.k_r = 1
        self.k_p = 0.00231
        self.q = 0
        self.n = 2
        self.p = 167

    def dW(self, q=1):
        dW = np.random.normal(0.0, 1.0 * 1) * np.sqrt(self.dt) # instert the initial condition
        return q * dW 

    def x(beta, y):
        return  beta * y
    
    def R_T(R_mono):
        return R_mono/100

    def R(self, x, R_T):
        return 1/(1 + pow((x/self.x_0),self.n)) * R_T

    def dy(self, r_y, y):
        return (self.alpha * self.k_p * r_y - self.gamma * y)*self.dt 
    
    def dr_y(self, R,r_y):
        return (self.k_r/(1 + (R/self.R_0))-self.gamma_r*r_y)*self.dt

    def dR_mono(self, R_mono):
        return self.k_R/self.gamma


    def simulation(self,TMG):
        
        with open(file, mode='a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
            for cell in range(self.cells):
                beta = 0.00123 * pow(TMG,0.6)
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
                for step in range(self.steps):
                    time_ += self.dt
                    counter += self.dt
                    y += self.dy(r_y, y)
                    r_y += self.dr_y(R, r_y) #+self.dW()
                    R_mono = self.k_R /self.gamma
                    x = beta * y
                    R_T = R_mono/100
                    R = self.R(x,R_T)
                    if step == self.steps - 1:   
                        data = Data(time=round(time_,4), 
                                    cell=round(cell,4), 
                                    beta= round(beta,4), 
                                    x=round(x,4), 
                                    y=round(y,4), 
                                    r_y=round(r_y,4), 
                                    R_mono=round(R_mono,4), 
                                    R_T=round(R_T,4), 
                                    R=round(R,4),
                                    T=TMG)
                        model = data.data_model()
                        writer.writerow(model)
                        data = None
        f.close()


if __name__ == '__main__':

    CRUD(file=file).delete()
    CRUD(name=file, schema=SCHEMA).create()
    sim = Simulation(steps=100000,cells=1)
    [sim.simulation(TMG=i/2) for i in range(120)]
