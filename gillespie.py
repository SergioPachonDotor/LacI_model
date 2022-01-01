from DataTools.ModelData import Data as data
from DataTools.ModelData import Schema as schema
from DataTools.dbTools import *
from Model.System import System as sys
from Model.Parameters import *
from DataTools.View import *
from tqdm import tqdm
from numpy.random import rand

import warnings
import numpy as np
import csv


init_state = STATE
SCHEMA = schema.get_classic_schema()
SCHEMA_DIV = Schema.get_division_schema()

def gillespie(TMG,Cells):
    
    SCHEMA = schema.get_classic_schema()

    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        
        for cell in tqdm(range(Cells)):

            beta = sys.beta_funct(tmg=TMG)
            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            dna = 1

            if init_state == 'off':
                y = 0
                r_y = 0

            elif init_state == 'on':
                y = 1000
                r_y = 10

            R_mono = 0
            R = 0.04
            tautime = 0
            reference_time = 0

            a_j = [i for i in range(2)]

            while tautime < tmax:

                a_j[0] = k_r * dna          
                a_j[1] = gamma_r * r_y

                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                if tautime + tau > reference_time: 
                    tau = reference_time - tautime

                    y += sys.dy(r_y, y, dt=tau)
                    R_mono += sys.dR_mono(R_mono, dt=tau)
                    x = sys.extmg(beta,y)
                    R_T = R_mono/100
                    R = sys.R_(x,R_T)

                    tautime = reference_time
                    data = Data(
                            time=round(tautime,4), 
                            cell=round(cell,4), 
                            beta=round(beta,4), 
                            Extracellular_TMG= TMG,
                            LacI_Tetramer=round(R_T,4),
                            Active_LacI=round(R,4), 
                            Permease=round(y,4), 
                            mRNA=round(r_y,4), 
                            LacI_monomer=round(R_mono,4),
                            Intracellular_TMG = round(x,4),
                            Promoter_State='None',
                            On_Time=on_time,
                            Off_Time=off_time
                            )

                    model = data.data_model()
                    writer.writerow(model)
                    data = None

                    reference_time += sampling_time
                
                else:
                    y += sys.dy(r_y, y, dt=tau)
                    R_mono += sys.dR_mono(R_mono, dt=tau)
                    x = sys.extmg(beta,y)
                    R_T = R_mono/100
                    R = sys.R_(x,R_T)
                    for i in range(len(a_j)):
                        sum_a += a_j[i]
                        if sum_a > dart:
                            q = i + 1
                            break
                    if q == 1:
                        # Transcription
                        r_y += 1

                    elif q == 2:
                        # mRNA degradation
                        r_y -= 1
                    tautime += tau

                if init_state == 'off':
                    if y >= 500 and on_counter == 0:
                        on_time = round(tautime,4)
                        on_counter = 1

                if init_state == 'on':
                    if y < 500 and off_counter == 0:
                        off_time = round(tautime,4)
                        off_counter = 1


    f.close()


def gillespie_full_noise(TMG,Cells):
    
    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        
        for cell in tqdm(range(Cells)):

            beta = sys.beta_funct(tmg=TMG)

            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            

            if init_state == 'off':
                y = 0
                r_y = 0
                dna = 0
                promoter = 'off'

            elif init_state == 'on':
                y = 1000
                r_y = 10
                dna = 1
                promoter = 'on'

            R_mono = 0
            R = 0.04
            R_T = 0
            dna = 1

            tautime = 0
            reference_time = 0

            a_j = [0 for i in range(6)]

            while tautime < tmax:

                if R_T != 0:

                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = alpha * k_p * r_y          # Protein Synthesis
                    a_j[3] = gamma * y                  # Protein Degradation
                    a_j[4] = landa * (R/R_T) 
                    a_j[5] = landa * (R_0/R_T)
                
                elif R_T == 0:

                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = alpha * k_p * r_y          # Protein Synthesis
                    a_j[3] = gamma * y                  # Protein Degradation
                    a_j[4] = 0
                    a_j[5] = 0
                
                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                if tautime + tau > reference_time: 
                    tau = reference_time - tautime

                    R_mono += sys.dR_mono(R_mono, dt=tau)
                    x = sys.extmg(beta,y)
                    R_T = R_mono/100
                    R = sys.R_(x,R_T)

                    tautime = reference_time
                    data = Data(
                            time = round(tautime,4), 
                            cell = round(cell,4), 
                            beta = round(beta,4), 
                            Extracellular_TMG= TMG,
                            LacI_Tetramer = round(R_T,4),
                            Active_LacI = round(R,4), 
                            Permease = round(y,4), 
                            mRNA = round(r_y,4), 
                            LacI_monomer = round(R_mono,4),
                            Intracellular_TMG = round(x,4),
                            Promoter_State = dna,
                            On_Time = on_time,
                            Off_Time = off_time
                            )

                    model = data.data_model()
                    writer.writerow(model)
                    data = None

                    reference_time += sampling_time
                
                else:

                    R_mono += sys.dR_mono(R_mono, dt=tau)
                    x = sys.extmg(beta,y)
                    R_T = R_mono/100
                    R = sys.R_(x,R_T)
                    for i in range(len(a_j)):
                        sum_a += a_j[i]
                        if sum_a > dart:
                            q = i + 1
                            break
                    if q == 1:
                        # Transcription
                        r_y += 1

                    elif q == 2:
                        # mRNA degradation
                        r_y -= 1

                    elif q == 3:
                        y += 1
                    
                    elif q == 4:
                        y -= 1
                    
                    elif q == 5:
                        if dna > 0:
                            dna -= 1
                        elif dna <= 0:
                            dna = 0
                        promoter = 'off'

                    elif q == 6:
                        if dna == 0:
                            dna += 1
                            
                        elif dna > 0:
                            dna = 1
                        promoter == 'on'

                    tautime += tau

                if init_state == 'off':
                    if y >= 500 and on_counter == 0:
                        on_time = round(tautime,4)
                        on_counter = 1

                if init_state == 'on':
                    if y < 500 and off_counter == 0:
                        off_time = round(tautime,4)
                        off_counter = 1
    f.close()


def gillespie_delta_tmg(TMG,Cells):
    
    with open(file, mode='a', newline='') as f:

        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter=',')
        
        for cell in tqdm(range(Cells)):

            beta = sys.beta_funct(tmg=TMG)

            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            

            if init_state == 'off':
                y = 0
                r_y = 0
                dna = 0
                promoter = 'off'

            elif init_state == 'on':
                y = 1000
                r_y = 10
                dna = 1
                promoter = 'on'

            R_mono = 0
            R = 0.04
            R_T = 0
            dna = 1

            tautime = 0
            reference_time = 0

            a_j = [0 for i in range(4)]

            while tautime < tmax:

                if R_T != 0:

                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = landa * (R/R_T) 
                    a_j[3] = landa * (R_0/R_T)
                
                elif R_T == 0:

                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = 0
                    a_j[3] = 0
                
                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                if tautime + tau > reference_time: 
                    tau = reference_time - tautime
                    y += sys.dy(r_y, y, dt=tau)
                    R_mono += sys.dR_mono(R_mono, dt=tau)
                    x = sys.extmg(beta,y)
                    R_T = R_mono/100
                    R = sys.R_(x,R_T)

                    tautime = reference_time
                    data = Data(
                            time = round(tautime,4), 
                            cell = round(cell,4), 
                            beta = round(beta,4), 
                            Extracellular_TMG= TMG,
                            LacI_Tetramer = round(R_T,4),
                            Active_LacI = round(R,4), 
                            Permease = round(y,4), 
                            mRNA = round(r_y,4), 
                            LacI_monomer = round(R_mono,4),
                            Intracellular_TMG = round(x,4),
                            Promoter_State = dna,
                            On_Time = on_time,
                            Off_Time = off_time
                            )

                    model = data.data_model()
                    writer.writerow(model)
                    data = None
                    
                    reference_time += sampling_time
                
                else:
                    y += sys.dy(r_y, y, dt=tau)
                    R_mono += sys.dR_mono(R_mono, dt=tau)
                    x = sys.extmg(beta,y)
                    R_T = R_mono/100
                    R = sys.R_(x,R_T)

                    for i in range(len(a_j)):
                        sum_a += a_j[i]
                        if sum_a > dart:
                            q = i + 1
                            break
                    if q == 1:
                        # Transcription
                        r_y += 1

                    elif q == 2:
                        # mRNA degradation
                        r_y -= 1
                    
                    elif q == 3:
                        if dna > 0:
                            dna -= 1
                        elif dna <= 0:
                            dna = 0
                        promoter = 'off'

                    elif q == 4:
                        if dna == 0:
                            dna += 1
                            
                        elif dna > 0:
                            dna = 1
                        promoter == 'on'

                    tautime += tau

                if init_state == 'off':
                    if y >= 250 and on_counter == 0:
                        on_time = round(tautime,4)
                        on_counter = 1

                if init_state == 'on':
                    if y < 250 and off_counter == 0:
                        off_time = round(tautime,4)
                        off_counter = 1
    
    f.close()


def gillespie_on_off(TMG,Cells, mode='single_switch'):

    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        
        for cell in tqdm(range(Cells)):

            beta = sys.beta_funct(tmg=TMG)

            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            for_dswitch = init_state
            flag = True
            
            if init_state == 'off':
                y = 0
                r_y = 0
                dna = 0
                promoter = 'off'

            elif init_state == 'on':
                y = 1000
                r_y = 10
                dna = 1
                promoter = 'on'

            R_mono = 0
            R = 0.04
            R_T = 0
            dna = 1

            tautime = 0

            a_j = [0 for i in range(4)]

            while flag:
                if R_T != 0:
                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = landa * (R/R_T) 
                    a_j[3] = landa * (R_0/R_T)
                
                elif R_T == 0:
                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = 0
                    a_j[3] = 0
                
                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                y = sys.dy_tau(r_y, y, tau)
                R_mono += sys.dR_mono(R_mono, dt=tau)
                x = sys.extmg(beta,y)
                R_T = R_mono/100
                R = sys.R_(x,R_T)

                for i in range(len(a_j)):
                    sum_a += a_j[i]
                    if sum_a > dart:
                        q = i + 1
                        break
                if q == 1:
                    # Transcription
                    r_y += 1

                elif q == 2:
                    # mRNA degradation
                    r_y -= 1
                
                elif q == 3:
                    if dna > 0:
                        dna -= 1
                    elif dna <= 0:
                        dna = 0
                    promoter = 'off'

                elif q == 4:
                    if dna == 0:
                        dna += 1
                        
                    elif dna > 0:
                        dna = 1
                    promoter == 'on'

                tautime += tau

                if mode == 'single_switch':
                    if init_state == 'off':

                        if TMG == 0:
                            break

                        if y >= 250 and on_counter == 0:

                            on_time = round(tautime,4)
                            data = Data(cell = cell,On_Time = on_time)                               
                            model = data.data_model()
                            writer.writerow(model)
                            data = None
                            on_counter = 1
                            flag = False

                    if init_state == 'on':

                        if y < 250 and off_counter == 0:

                            off_time = round(tautime,4)
                            data = Data(cell = cell,Off_Time = off_time)                               
                            model = data.data_model()
                            writer.writerow(model)
                            data = None
                            off_counter = 1
                            flag = False
                        
                elif mode == 'double_switch':

                    if for_dswitch == 'off':

                        if TMG < 5:
                            break

                        if y > 250 and on_counter == 0:

                            on_time = round(tautime,4)
                            data = Data( 
                                    cell = round(cell,4), 
                                    beta = round(beta,4), 
                                    Extracellular_TMG= TMG,
                                    LacI_Tetramer = round(R_T,4),
                                    Active_LacI = round(R,4), 
                                    Permease = round(y,4), 
                                    mRNA = round(r_y,4), 
                                    LacI_monomer = round(R_mono,4),
                                    Intracellular_TMG = round(x,4),
                                    Promoter_State = dna,
                                    On_Time = on_time,
                                    )                            
                            model = data.data_model()
                            writer.writerow(model)

                            data = None
                            on_counter = 1
                            for_dswitch = 'on'

                    if for_dswitch == 'on':

                        if y < 250 and off_counter == 0:

                            off_time = round(tautime,4)
                            data = Data( 
                                    cell = round(cell,4), 
                                    beta = round(beta,4), 
                                    Extracellular_TMG= TMG,
                                    LacI_Tetramer = round(R_T,4),
                                    Active_LacI = round(R,4), 
                                    Permease = round(y,4), 
                                    mRNA = round(r_y,4), 
                                    LacI_monomer = round(R_mono,4),
                                    Intracellular_TMG = round(x,4),
                                    Promoter_State = dna,
                                    Off_Time = off_time
                                    )                               
                            model = data.data_model()
                            writer.writerow(model)

                            data = None
                            off_counter = 1
                            for_dswitch = 'off'
                        
                    if on_counter == 1 and off_counter == 1:
                        flag = False
 
    f.close()


def calculate_propensity(propensity):
    try:
        warnings.simplefilter("ignore")
        p = -(1/propensity) * np.log(rand())
        return p
    except ZeroDivisionError:
        return 0
        

def calculate_division_tau(propensity, mu, size):
    # try:
        warnings.simplefilter("ignore")
        p = (1/mu) * np.log(1 - (mu * np.log(np.random.rand())/(propensity * size)))
        return p

    # except ZeroDivisionError:
    #     return 0
    

def simulate_adder_on_off(Cells, TMG, birth_size, division_time):
        # np.random.seed(212121)
        # np.random.seed(123456)

        """ 
            Performs Cesar Nieto et-al Gillespie algorithm
            for chemical reaction and division simulations.
        """
        with open(file, mode='a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames= SCHEMA_DIV, delimiter='|')

            for cell in tqdm(range(1, Cells+1)):
                size = birth_size
                mu = np.log(2)/division_time
                reference_division_time = (1/mu) * np.log(((birth_size + np.random.normal(loc=1.0, scale=0.05))/birth_size))

                beta = sys.beta_funct(tmg=TMG)

                on_counter = 0
                on_time = 0
                off_counter = 0
                off_time = 0
                

                if init_state == 'off':
                    y = 0
                    r_y = 0
                    dna = 0
                    promoter = 'off'

                elif init_state == 'on':
                    y = 1000
                    r_y = 10
                    dna = 1
                    promoter = 'on'

                R_mono = 0
                R = 6
                R_T = 6
                dna = 1

                tautime = 0
                reference_time = 0

                a_j = [0 for i in range(4)]

                while tautime < tmax:
                    if R_T > 0:

                        a_j[0] = k_r * dna                  # Transcription  
                        a_j[1] = gamma_r * r_y              # mRNA Degradation
                        a_j[2] = landa * (R/R_T) 
                        a_j[3] = landa * (R_0/R_T)
                    
                    elif R_T <= 0:

                        a_j[0] = k_r * dna                  # Transcription  
                        a_j[1] = gamma_r * r_y              # mRNA Degradation
                        a_j[2] = 0
                        a_j[3] = 0
                    

                    tau_arr = [calculate_propensity(i) for i in a_j]
                    tau_arr[0] = calculate_division_tau(propensity=a_j[0], mu=mu, size=size)

                    if r_y <= 0:
                        tau_arr[1] = 100                   

                    
                    tau_arr.append(reference_division_time)

                    tau = np.min(tau_arr)

                    if tautime + tau < reference_time:          
                        # print(tautime)
                        y = sys.dy_tau(r_y, y, tau)
                        R_mono += sys.dR_mono(R_mono, dt=tau)
                        x = sys.extmg(beta,y)
                        R_T = R_mono/100
                        R = sys.R_(x,R_T)

                        q = np.argmin(tau_arr) + 1

                        if q == 1:
                            # Transcription
                            r_y += 1

                        elif q == 2:
                            # mRNA degradation
                            if r_y >= 1:
                                r_y -= 1
                        
                        elif q == 3:
                            if dna > 0:
                                dna -= 1
                            elif dna <= 0:
                                dna = 0
                            promoter = 'off'

                        elif q == 4:
                            if dna == 0:
                                dna += 1
                                
                            elif dna > 0:
                                dna = 1
                            promoter == 'on'
                        
                        elif q == 5:
                            # Set all division conditions like dilution of species
                            size /= 2
                            reference_division_time = (1/mu) * np.log((size + np.random.normal(loc=1.0, scale=0.05))/size)   
                            # dilute_species
                            r_y /= 2
                            y /= 2
                            R_mono /= 2
                            R /= 2
                            R_T /= 2
                            # print('Divide')

                        size *= np.exp(mu*tau)

                        tautime += tau   

                        if init_state == 'off':
                            if y >= 250 and on_counter == 0:
                                on_time = round(tautime,4)
                                on_counter = 1

                        if init_state == 'on':
                            if y < 250 and off_counter == 0:
                                off_time = round(tautime,4)
                                off_counter = 1
                            
                    elif tautime + tau > reference_time:
                        y = sys.dy_tau(r_y, y, tau)
                        R_mono += sys.dR_mono(R_mono, dt=tau)
                        x = sys.extmg(beta,y)
                        R_T = R_mono/100
                        R = sys.R_(x,R_T)

                        reference_division_time -= (reference_time - tautime)
                        size *= np.exp(mu * (reference_time - tautime))
                        tautime = reference_time  

                        # ___Save Data___ # 
                        data = Data(
                                time = round(tautime,4), 
                                cell = round(cell,4), 
                                beta = round(beta,4), 
                                Extracellular_TMG= TMG,
                                LacI_Tetramer = round(R_T,4),
                                Active_LacI = round(R,4), 
                                Permease = round(y,4), 
                                mRNA = round(r_y,4), 
                                LacI_monomer = round(R_mono,4),
                                Intracellular_TMG = round(x,4),
                                Promoter_State = dna,
                                On_Time = on_time,
                                Off_Time = off_time,
                                size = round(size, 4)
                                )
                        model = data.data_division_model()
                        writer.writerow(model)
                        data = None
                        reference_time += sampling_time

        f.close()


def division_gillespie(TMG, Cells, mode='single_switch'):

    with open(file, mode='a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames= SCHEMA, delimiter='|')
        
        for cell in tqdm(range(Cells)):

            beta = sys.beta_funct(tmg=TMG)

            on_counter = 0
            on_time = 0
            off_counter = 0
            off_time = 0
            for_dswitch = init_state
            flag = True
            
            if init_state == 'off':
                y = 0
                r_y = 0
                dna = 0
                promoter = 'off'

            elif init_state == 'on':
                y = 1000
                r_y = 10
                dna = 1
                promoter = 'on'

            R_mono = 0
            R = 0.04
            R_T = 0
            dna = 1

            tautime = 0

            a_j = [0 for i in range(4)]

            while flag:
                if R_T != 0:
                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = landa * (R/R_T) 
                    a_j[3] = landa * (R_0/R_T)
                
                elif R_T == 0:
                    a_j[0] = k_r * dna                  # Transcription  
                    a_j[1] = gamma_r * r_y              # mRNA Degradation
                    a_j[2] = 0
                    a_j[3] = 0
                
                a_total = np.sum(a_j)

                tau = (1/a_total)*np.log(1/rand())

                dart = a_total * rand()
                sum_a = 0
                q = 0

                y = sys.dy_tau(r_y, y, tau)
                R_mono += sys.dR_mono(R_mono, dt=tau)
                x = sys.extmg(beta,y)
                R_T = R_mono/100
                R = sys.R_(x,R_T)

                for i in range(len(a_j)):
                    sum_a += a_j[i]
                    if sum_a > dart:
                        q = i + 1
                        break
                if q == 1:
                    # Transcription
                    r_y += 1

                elif q == 2:
                    # mRNA degradation
                    r_y -= 1
                
                elif q == 3:
                    if dna > 0:
                        dna -= 1
                    elif dna <= 0:
                        dna = 0
                    promoter = 'off'

                elif q == 4:
                    if dna == 0:
                        dna += 1
                        
                    elif dna > 0:
                        dna = 1
                    promoter == 'on'

                tautime += tau

                if mode == 'single_switch':
                    if init_state == 'off':

                        if TMG == 0:
                            break

                        if y >= 250 and on_counter == 0:

                            on_time = round(tautime,4)
                            data = Data(cell = cell,On_Time = on_time)                               
                            model = data.data_model()
                            writer.writerow(model)
                            data = None
                            on_counter = 1
                            flag = False

                    if init_state == 'on':

                        if y < 250 and off_counter == 0:

                            off_time = round(tautime,4)
                            data = Data(cell = cell,Off_Time = off_time)                               
                            model = data.data_model()
                            writer.writerow(model)
                            data = None
                            off_counter = 1
                            flag = False
                        
                elif mode == 'double_switch':

                    if for_dswitch == 'off':

                        if TMG < 5:
                            break

                        if y > 250 and on_counter == 0:

                            on_time = round(tautime,4)
                            data = Data( 
                                    cell = round(cell,4), 
                                    beta = round(beta,4), 
                                    Extracellular_TMG= TMG,
                                    LacI_Tetramer = round(R_T,4),
                                    Active_LacI = round(R,4), 
                                    Permease = round(y,4), 
                                    mRNA = round(r_y,4), 
                                    LacI_monomer = round(R_mono,4),
                                    Intracellular_TMG = round(x,4),
                                    Promoter_State = dna,
                                    On_Time = on_time,
                                    )                            
                            model = data.data_model()
                            writer.writerow(model)

                            data = None
                            on_counter = 1
                            for_dswitch = 'on'

                    if for_dswitch == 'on':

                        if y < 250 and off_counter == 0:

                            off_time = round(tautime,4)
                            data = Data( 
                                    cell = round(cell,4), 
                                    beta = round(beta,4), 
                                    Extracellular_TMG= TMG,
                                    LacI_Tetramer = round(R_T,4),
                                    Active_LacI = round(R,4), 
                                    Permease = round(y,4), 
                                    mRNA = round(r_y,4), 
                                    LacI_monomer = round(R_mono,4),
                                    Intracellular_TMG = round(x,4),
                                    Promoter_State = dna,
                                    Off_Time = off_time
                                    )                               
                            model = data.data_model()
                            writer.writerow(model)

                            data = None
                            off_counter = 1
                            for_dswitch = 'off'
                        
                    if on_counter == 1 and off_counter == 1:
                        flag = False
 
    f.close()


if __name__ == '__main__':

    if algorithm == 'gillespie':

        file = f'./simulation_data/gillespie/gillespie_{init_state}_{tmg}_cells_{cells}.csv'
        Tools(file=file).delete()
        Tools(name=file, schema=SCHEMA).create()
        gillespie(TMG=tmg,Cells=cells)
        main_view()
    
    elif algorithm == 'gillespie_full_noise':

        file = f'./simulation_data/gillespie_full_noise/gillespie_full_noise_{init_state}_{tmg}_cells_{cells}.csv'
        Tools(file=file).delete()
        Tools(name=file, schema=SCHEMA).create()
        gillespie_full_noise(TMG=tmg,Cells=cells)
        main_view()

    elif algorithm == 'gillespie_delta_tmg':
        factor = 5
        for tmg_concentration in tqdm(range(tmg_range)):

            file = f'./simulation_data/gillespie_delta_tmg/gillespie_delta_tmg_{init_state}_{tmg_concentration*factor}_cells_{cells}.csv'
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()
            gillespie_delta_tmg(TMG=(tmg_concentration*factor),Cells=cells)

        main_view()

    elif algorithm == 'gillespie_on_off':
        factor = 5
        for tmg_concentration in tqdm(range(1)):

            file = f'./simulation_data/gillespie_on_off/gillespie_on_off_state_{init_state}_tmg_{tmg_concentration * factor}_cells_{cells}.csv'
            
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()

            gillespie_on_off(TMG=tmg_concentration * factor, Cells=cells, mode='single_switch')

        main_view()

    elif algorithm == 'gillespie_on_off_double':

        factor = 5

        for tmg_concentration in tqdm(range(1,tmg_range)):
            file = f'./simulation_data/gillespie_on_off_double/gillespie_on_off_double_state_{init_state}_tmg_{tmg_concentration * factor}_cells_{cells}.csv'
            
            Tools(file=file).delete()
            Tools(name=file, schema=SCHEMA).create()

            gillespie_on_off(TMG=tmg_concentration * factor, Cells=cells, mode='double_switch')

        main_view()

    elif algorithm == 'simulate_adder':
        birth_size = 1
        # division_time = np.random.normal(loc=2/2, scale=0.1*(2/2))
        division_time = 18
        file = f'./simulation_data/division_adder/adder_tmg_{tmg}_bs_{birth_size}_divtime_{division_time}_{init_state}_{cells}.csv'
        Tools(file=file).delete()
        Tools(name=file, schema=SCHEMA_DIV).create()
        simulate_adder_on_off(TMG=tmg, Cells=cells, birth_size=birth_size, division_time=division_time)
        main_view()