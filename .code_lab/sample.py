import numpy as np
import os

def create_file():
    try:
        with open('cell_1.py', mode='x') as f:
            f.write('Hello World')
            f.close()
    except FileExistsError as file_exist:
        print(file_exist)

x = 2
y = x * 4
m = 'sst'


