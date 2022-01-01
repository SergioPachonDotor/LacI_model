#_____SCHEMA_____

class Schema:

    def __init__(self) -> None:

        self.classic:list[str] = [
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
                        'Promoter_State',
                        'On_Time',
                        'Off_Time'
                    ]

        self.division: list[str] = [
                            'Time',
                            'Cell',
                            'Size',
                            'beta',
                            'Intracellular_TMG',
                            'Permease',
                            'mRNA',
                            'LacI_monomer',
                            'LacI_tetramer',
                            'Active_LacI',
                            'Extracellular_TMG',
                            'Promoter_State',
                            'On_Time',
                            'Off_Time'
                        ]

        self.new_schema: list[str] = []

    def get_classic_schema(self):
        return self.classic

    def get_division_schema(self):
        return self.division

    def set_new_schema(self, new_schema: list[str]):
        self.new_schema = new_schema

class Data(Schema):
    
    def __init__(self, 
                 time=None, 
                 cell=None, 
                 beta=None, 
                 Extracellular_TMG=None,
                 LacI_Tetramer=None,
                 Active_LacI=None, 
                 Permease=None, 
                 mRNA=None, 
                 LacI_monomer=None,
                 Intracellular_TMG = None,
                 Promoter_State=None,
                 On_Time=None,
                 Off_Time=None,
                 size=None
                 ):

        super().__init__()
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
        self.On_Time = On_Time
        self.Off_Time = Off_Time
        self.size = size

    def data_model(self):
        model = {
                 self.classic[0]:f'{self.time}', 
                 self.classic[1]:f'{self.cell}', 
                 self.classic[2]:f'{self.beta}', 
                 self.classic[3]:f'{self.x}',
                 self.classic[4]:f'{self.y}',
                 self.classic[5]:f'{self.r_y}',
                 self.classic[6]:f'{self.R_mono}',
                 self.classic[7]:f'{self.R_T}',
                 self.classic[8]:f'{self.R}',
                 self.classic[9]:f'{self.T}',
                 self.classic[10]:f'{self.promoter}',
                 self.classic[11]:f'{self.On_Time}',
                 self.classic[12]:f'{self.Off_Time}'
                }
        return model

    def data_division_model(self):
        model = {
                self.division[0]:f'{self.time}', 
                self.division[1]:f'{self.cell}', 
                self.division[2]:f'{self.size}',
                self.division[3]:f'{self.beta}',
                self.division[4]:f'{self.x}',
                self.division[5]:f'{self.y}',
                self.division[6]:f'{self.r_y}',
                self.division[7]:f'{self.R_mono}',
                self.division[8]:f'{self.R_T}',
                self.division[9]:f'{self.R}',
                self.division[10]:f'{self.T}',
                self.division[11]:f'{self.promoter}',
                self.division[12]:f'{self.On_Time}',
                self.division[13]:f'{self.Off_Time}'
        }

        return model