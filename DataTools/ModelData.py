#_____SCHEMA_____

SCHEMA = [
            'Time',
            'Cell',
            'beta',
            'Intracellular_TMG',
            'Permease',
            'mRNA',
            'LacI_monomer',
            'LacI_tetramer',
            'Active_LacI',
            'Extracellular_TMG',
            'Promoter_State'
         ]

class Data:
    
    def __init__(self, 
                 time=float, 
                 cell=int, 
                 beta=float, 
                 Extracellular_TMG=float,
                 LacI_Tetramer=float,
                 Active_LacI=float, 
                 Permease=float, 
                 mRNA=float, 
                 LacI_monomer=float,
                 Intracellular_TMG = float,
                 Promoter_State=str):

        #Time (dt)| Sample (trayectory id) | Variables
        self.time = time
        self.cell = cell
        self.beta = beta
        self.x = Intracellular_TMG
        self.R_T = LacI_Tetramer
        self.R = Active_LacI
        self.y = Permease
        self.r_y = mRNA
        self.R_mono = LacI_monomer
        self.T = Extracellular_TMG
        self.promoter = Promoter_State

    def data_model(self):
        model = {
                 SCHEMA[0]:f'{self.time}', 
                 SCHEMA[1]:f'{self.cell}', 
                 SCHEMA[2]:f'{self.beta}', 
                 SCHEMA[3]:f'{self.x}',
                 SCHEMA[4]:f'{self.y}',
                 SCHEMA[5]:f'{self.r_y}',
                 SCHEMA[6]:f'{self.R_mono}',
                 SCHEMA[7]:f'{self.R_T}',
                 SCHEMA[8]:f'{self.R}',
                 SCHEMA[9]:f'{self.T}',
                 SCHEMA[10]:f'{self.promoter}'
                }
        return model