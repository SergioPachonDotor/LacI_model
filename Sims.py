from DataTools.dbTools import Schema
from DataTools.dbTools import Tools
from Model.System import System as sys
from Model.System import DivisionSystem as divsys
from tqdm import tqdm
import numpy as np

import csv

class Simulation:
    
    def __init__(
                    self, 
                    TMG: float, 
                    cells: int, 
                    file: str, 
                    tmax: int = 720, 
                    steps: int = 1000,
                    state: str = 'on',
                    tmg_range: list = [0, 5, 10, 15, 20, 25, 30, 35, 40],
                    ) -> None:

        self.tmg: float = TMG
        self.cells: int = cells
        self.file: str = file
        # self.dt: float = tmax/steps
        self.dt: float = 0.1
        self.steps: int = steps
        self.tmax: int = tmax
        self.state: str = state
        self.delta_tmg: list = tmg_range

    def gillespie(self):
        SCHEMA = Schema().get_classic_schema()

        Tools(file=self.file).delete()
        Tools(file=self.file, schema=SCHEMA).create()

        with open(f'{self.file}.csv', mode='a', newline='') as f:

            writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter=',')
            
            for cell in tqdm(range(1, self.cells + 1)):
                
                system = sys(Extracellular_TMG=self.tmg, Promoter=self.state)

                system.set_beta(tmg=self.tmg)
                system.set_promoter_state()
                tautime = 0
                reference_time = 0

                while tautime < self.tmax:
                    system.set_reactions_classic_gillespie()
                    q, tau = system.set_q_tau()

                    if tautime + tau < reference_time:
                        system.solve_eqs(tau)
                        system.react(q=q)
                        tautime += tau
                        system.switch_counter(tautime=tautime)
        
                    elif tautime + tau > reference_time: 
                        tau = reference_time - tautime
                        system.solve_eqs(tau)
                        tautime = reference_time
                        
                        model = system.model_to_save(tautime=tautime, cell=cell, tmg=self.tmg)
                        writer.writerow(model)
                        
                        reference_time += self.dt

    def gillespie_delta_tmg(self):
        for tmg in tqdm(self.delta_tmg):
            self.tmg = tmg
            self.file = f'{self.file}_{tmg}'
            self.gillespie()

    def simulate_adder(self, birth_size, division_time):
        # np.random.seed(212121)
        # np.random.seed(123456)

        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reaction and division simulations.
        """

        SCHEMA_DIV = Schema().get_division_schema()

        Tools(file=self.file).delete()
        Tools(file=self.file, schema=SCHEMA_DIV).create()

        with open(f'{self.file}.csv', mode='a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames= SCHEMA_DIV, delimiter=',')

            for cell in tqdm(range(1, self.cells + 1)):
                div_sys = divsys(Extracellular_TMG=self.tmg, Promoter=self.state, birth_size=birth_size, division_time=division_time)                
                div_sys.set_beta(tmg=self.tmg)
                div_sys.set_promoter_state()
                
                tautime = 0
                reference_time = 0

                while tautime < self.tmax:
                    div_sys.set_reactions_classic_gillespie()
                    q, tau = div_sys.set_calculated_sorted_tau()

                    if tautime + tau < reference_time:          
                        div_sys.solve_eqs(tau)
                        div_sys.division_react(q=q)
                        div_sys.update_tau(tau=tau)

                        tautime += tau   

                        div_sys.switch_counter(tautime=tautime)
                            
                    elif tautime + tau > reference_time:
                        
                        div_sys.solve_eqs(tau)

                        div_sys.reference_division_time -= (reference_time - tautime)
                        div_sys.update_tau(tau=(reference_time - tautime))
                        tautime = reference_time  

                        # ___Save Data___ # 
                        model = div_sys.division_model_to_save(tautime=tautime, cell=cell, tmg=self.tmg)

                        writer.writerow(model)
                        reference_time += self.dt

    def get_time_distribution(self, mode='single'):

        SCHEMA = Schema().get_classic_time_distribution_schema()
        Tools(file=self.file).delete()
        Tools(file=self.file, schema=SCHEMA).create()

        with open(f'{self.file}.csv', mode='a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter=',')
            
            for cell in tqdm(range(1, self.cells + 1)):

                system = sys(Extracellular_TMG=self.tmg, Promoter=self.state)

                system.set_beta(tmg=self.tmg)
                system.set_promoter_state()
                tautime = 0
                flag = True

                while flag:
                    system.set_reactions_classic_gillespie()
                    q, tau = system.set_q_tau()
                    system.solve_eqs(tau)
                    system.react(q)

                    tautime += tau

                    if mode == 'single':
                        model = system.get_single_switch(tautime, cell, self.tmg)
                        if model is None:
                            pass
                        elif model is not None:
                            writer.writerow(model)
                            flag = False

                    elif mode == 'double':
                        model = system.get_double_switch(tautime, cell, self.tmg)
                        if model is None:
                            pass
                        elif model is not None:
                            writer.writerow(model)

                        if system.on_counter == 1 and system.off_counter == 1:
                            flag = False

    def get_time_division_distribution(self, mode='single'):
            
            SCHEMA = Schema().get_division_time_distribution_schema()
            Tools(file=self.file).delete()
            Tools(file=self.file, schema=SCHEMA).create()
    
            with open(f'{self.file}.csv', mode='a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter=',')
                
                for cell in tqdm(range(1, self.cells + 1)):
    
                    div_sys = divsys(Extracellular_TMG=self.tmg, Promoter=self.state)
    
                    div_sys.set_beta(tmg=self.tmg)
                    div_sys.set_promoter_state()
                    tautime = 0
                    flag = True
    
                    while flag:
                        div_sys.set_reactions_classic_gillespie()
                        q, tau = div_sys.set_calculated_sorted_tau()
                        div_sys.solve_eqs(tau)
                        div_sys.division_react(q)
                        div_sys.update_tau(tau=tau)
    
                        tautime += tau
    
                        if mode == 'single':
                            model = div_sys.get_single_switch(tautime, cell, self.tmg)
                            if model is None:
                                pass
                            elif model is not None:
                                writer.writerow(model)
                                flag = False
    
                        elif mode == 'double':
                            model = div_sys.get_double_switch(tautime, cell, self.tmg)
                            if model is None:
                                pass
                            elif model is not None:
                                writer.writerow(model)
    
                            if div_sys.on_counter == 1 and div_sys.off_counter == 1:
                                flag = False

if __name__ == '__main__':
    sim = Simulation(
                TMG=0.2, 
                cells=100, 
                file='.tests/test_OOP', 
                tmax=100, 
                steps=1000,
                state='on'
                )
    sim.gillespie()