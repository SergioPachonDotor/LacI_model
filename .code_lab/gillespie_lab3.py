import numpy as np 
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd 

from tqdm import tqdm
from numpy.random import rand

class Gillespie:
    def __init__(self, tmax=10, species = dict, propensities = dict, reactions= dict):
        self.sampling_time = 0.001          # Sampling Time
        self.time = 0                       # Iitial Time
        self.reference_time = 0             # Reference Time
        self.tmax = tmax                    # Total Time
        self.propensities = propensities    # {Reaction: Propensity}
        self.species = species              # {Species: Quantity}
        self.reactions = reactions          # {Reaction: Species}

    def simulate(self):
        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """
        while self.time < self.tmax:

            reactions_list = list(self.propensities.keys())
            propenisties_list = list(self.propensities.values())

            reaction_times = lambda propensity: -(1/(propensity)) * np.log(rand())
            calculated_propensities = map(reaction_times, propenisties_list)

            tau = np.min(propenisties_list)     # tau time
            q = np.argmin(propenisties_list)    # Reaction that occurs

            if self.time + tau < self.reference_time:
                pass

                self.time += tau
                
            else:
                self.time = self.reference_time  

                # Save Data
                # tarray.append(t)
                # protarr.append(proteins)

                self.reference_time += self.sampling_time   

            reactions_activity = dict(zip(reactions_list, calculated_propensities))

        print(reactions_activity)

if __name__ == '__main__':
    
    
    species = {
                'dna': 1, 
                'rna': 0, 
                'protein': 0, 
                'complex': 0
              }

    reactions = {
                'create': 1,
                'destroy': 1
                }

    propensities = {
                'transcription':       0.005, 
                'translation':         0.167, 
                'mrna_degradation':    np.log(2)/150, 
                'protein_degradation': np.log(2)/3600,
                'complex_creation':    0.00001,
                'complex_degradation': 0.005
                }

    print(list(reactions.values()))
    # sim_3 = Gillespie(tmax=10, species = species, propensities= propensities, reactions= reactions)
    # sim_3.simulate()


class Reaction:

    def __init__(self, name=str, species=list, propensities=list, action=str):
        self.name = name
        self.species = species
        self.propensities = propensities
        self.action = action

    def model(self):
        reaction_model = True
        return reaction_model

    def birth(self):
        return
    
    def death(self):
        return