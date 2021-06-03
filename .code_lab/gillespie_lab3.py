import numpy as np 
import matplotlib.pyplot as plt 
import copy
import warnings

class Gillespie:

    def __init__(self, tmax=10, sampling_time = 0.01,reaction_model = object):
        self.sampling_time = sampling_time      # Sampling Time
        self.time = 0                           # Iitial Time
        self.reference_time = 0                 # Reference Time
        self.tmax = tmax                        # Total Time
        self._reaction_model_copy = copy.copy(reaction_model)
        self.reaction_model = reaction_model    # Object of reaction model
        self.tau = 0
        self.tarr = [0]
        self.protarr = [0]

    def simulate(self):
        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reactions simulations.
        """

        while self.time < self.tmax:

            propensities = self.reaction_model.calculate_propensities()
            propenisties_list = list(propensities.values())
            possible_reactions = None

            warnings.simplefilter("ignore")
            reaction_times = lambda propensity, random_num: -(1/propensity) * np.log(random_num)
            random_array = [np.random.rand() for k in range(len(propenisties_list))]

            calculated_propensities = list(map(reaction_times, propenisties_list, random_array))
            reactions_dict = dict(enumerate(list(self.reaction_model.show_reactions().keys())))
            
            self.tau = np.min(calculated_propensities)     # tau time
            q = np.argmin(calculated_propensities)    # Reaction that occurs

            if self.time + self.tau < self.reference_time:
                possible_reactions = self.reaction_model.show_q()[reactions_dict[q]]

                create_species = possible_reactions['create']
                [self.reaction_model.create(name=species) for species in create_species if species != None]

                destroy_species = possible_reactions['destroy']
                [self.reaction_model.destroy(name=species) for species in destroy_species if species != None]


                self.time += self.tau
                
            else:
                self.time = self.reference_time  

                self.tarr.append(self.time)
                self.protarr.append(self.reaction_model.show_species()['protein'])

                self.reference_time += self.sampling_time
        
        plt.plot(self.tarr, self.protarr)
        plt.savefig('test.jpg')
        self.reaction_model = self._reaction_model_copy

class ReactionModel:
    
    def __init__(self, reactions=dict,species=dict, propensities=dict, q=dict):
        self.reactions = reactions
        self.species = species 
        self.propensities = propensities
        self.q = q
    
    def show_species(self):
        return self.species

    def show_reactions(self):
        return self.reactions

    def show_propensities(self):
        return self.propensities
    
    def show_q(self):
        return self.q
    
    def calculate_propensities(self):
        
        complete_propensity = list(self.reactions.values())
        reaction_names = self.reactions.keys()
        species_copy = copy.copy(self.species)
        species_copy.update(self.propensities)

        species_calculus = []

        for i in range(len(complete_propensity)):
            propensities_array = []
            species_calculus.append(propensities_array)

            for k in complete_propensity[i]:
                propensities_array.append(species_copy[k])
                
        propensities_result = list(map(np.prod, species_calculus))
        
        return dict(zip(reaction_names, propensities_result))
    
    def create(self, name=str):
        self.species[name] += 1
    
    def destroy(self, name=str):
        self.species[name] -= 1
    
    def set_q(self):
        pass

if __name__ == '__main__':

    species = {
                'dna': 1, 
                'rna': 0, 
                'protein': 0, 
                'complex': 0
                }

    propensities = {
                    'trc_c':    0.3, 
                    'trl_c':    10.02, 
                    'pdeg_d':   np.log(2)/2.5, 
                    'cmplx_c':  np.log(2)/30,
                    'cmplx_d':  0.3
                }

    reactions = {
                'transcription': ['trc_c', 'dna'], 
                'translation':   ['trl_c', 'rna'], 
                'degradation':   ['pdeg_d', 'protein'], 
                'complex_c':     ['cmplx_c', 'dna', 'protein'],
                'complex_d':     ['cmplx_d', 'complex']
                }

    q = {
            'transcription': {'create': ['rna'],            'destroy': [None]}, 
            'translation':   {'create': ['protein'],        'destroy': [None]}, 
            'degradation':   {'create': [None],             'destroy': ['protein']}, 
            'complex_c':     {'create': ['complex'],        'destroy': ['dna', 'protein']},
            'complex_d':     {'create': ['dna', 'protein'], 'destroy': ['complex']}
        }  

    model = ReactionModel(reactions=reactions, species=species, propensities=propensities, q=q)
    Gillespie(tmax=720, sampling_time=0.1,reaction_model = model).simulate(); 

    species_2 = {
                    'protein': 0
                    }

    propensities_2 = {
                        'kr': 100, 
                        'gamma':10
                        }
    reactions_2 = {
                    'translation':['kr'], 
                    'degradation':['gamma', 'protein']
                    }
    q_2 = {
        'translation': {'create': ['protein'],  'destroy': [None]}, 
        'degradation': {'create': [None],       'destroy': ['protein']}
        }

    # model_2 = ReactionModel(reactions=reactions_2, species=species_2, propensities=propensities_2, q=q_2)
    # Gillespie(tmax=10, sampling_time=0.01, reaction_model=model_2).simulate()