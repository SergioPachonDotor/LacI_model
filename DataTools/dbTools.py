from re import S
from Model.Parameters import *
from DataTools.ModelData import *
import csv
import os

class Tools:

    def __init__(self, data=dict, file='simulation_data.csv', schema=list):
        self.data = data
        self.file = f'{file}.csv'
        self.schema = schema

    def create(self):
        try:
            with open(self.file, mode='x') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self.schema, delimiter=',')
                writer.writeheader()
                csvfile.close()
        except:
            print(f'File "{self.file}" already exist')
            pass

    def update(self):      
        with open(self.file, mode='a', newline='') as f:
            SCHEMA = Schema().get_classic_schema()
            writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter=',')
            writer.writerow(self.data)
            f.close()

    def delete(self):
        try:
            os.remove(self.file)
        except:
            pass
