from Model.Parameters import *
from DataTools.ModelData import *
import csv
import os

class Tools:

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
