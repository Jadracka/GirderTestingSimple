# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:54:43 2021

@author: jbarker
"""

import math as m
"""

   _____             __ _                       _   _                __ _ _      
  / ____|           / _(_)                     | | (_)              / _(_) |     
 | |     ___  _ __ | |_ _  __ _ _   _ _ __ __ _| |_ _  ___  _ __   | |_ _| | ___ 
 | |    / _ \| '_ \|  _| |/ _` | | | | '__/ _` | __| |/ _ \| '_ \  |  _| | |/ _ \
 | |___| (_) | | | | | | | (_| | |_| | | | (_| | |_| | (_) | | | | | | | | |  __/
  \_____\___/|_| |_|_| |_|\__, |\__,_|_|  \__,_|\__|_|\___/|_| |_| |_| |_|_|\___|
                           __/ |                                                 
                          |___/                                                  
"""
"""Data analysis tools"""
"""CAN BE CHANGED"""
Which_epochs = (0,2) #the comma must stay, otherwise the variable will be int
Names_of_magnets = ['PQK36','PQL6','PQK62']
Instruments_6DoF = False

Pico_StDev_basic = 0.00000001 #mm, PicoScale Standard Deviation
IFM_StDev_basic = (0.0004,0.00015) #mm, LT IFM, based on Leica's white paper
ADM_StDev_basic = 0.010 #mm, from 2-5m based on Leica's typical errors
Hz_StDev_basic = 0.15 #mgon same like Leica TM50
V_StDev_basic = 0.15 #mgon same like Leica TM50
Constraint_StDev_basic = 0.000001 #mm


Max_diff_from_line = 1.7 #mm - maximum distance from line, report the excess
LSM_Threshold = 1e-6
LSM_Max_iterations = 10

Sigma_0 = 0.000001



# Epoch specific factors
Epoch_factors = {0:{}, 1:{}, 2:{}}
Epoch_factors[0]['Pico'] = 1
Epoch_factors[0]['IFM'] = 5.75
Epoch_factors[0]['ADM'] = 0.90
Epoch_factors[0]['Hz'] = 3.95
Epoch_factors[0]['V'] = 12.66
Epoch_factors[0]['Con'] = 1

Epoch_factors[1]['Pico'] = 1
Epoch_factors[1]['IFM'] = 6.32
Epoch_factors[1]['ADM'] = 0.90
Epoch_factors[1]['Hz'] = 8.90
Epoch_factors[1]['V'] = 48.8
Epoch_factors[1]['Con'] = 1

Epoch_factors[2]['Pico'] = 1
Epoch_factors[2]['IFM'] = 5.52
Epoch_factors[2]['ADM'] = 1.33
Epoch_factors[2]['Hz'] = 8.9
Epoch_factors[2]['V'] = 37
Epoch_factors[2]['Con'] = 1

Print_FIDs = False


# Print troubleshooting messages:
Line_differences_checking = False
Print_all_troubleshooting = False
Using_nominal_compare = False
Print_typos =  False #Printing error messages
Print_2F_checks = False
Print_real2nominal_checks = False
Print_epoch_checks = False


"""Measured Lines
Standard SA format with spaces as delimiters and no comments:
Group aka Line name 'space' Point name 'space' Hz [gon] 'space' V [gon] 'space' Sd [mm]
"""
Epochs_dictionary = {'LoS':{},'Pol':{},'Coord':{}}
Epochs_dictionary['LoS'][-1] = "LoS_measurements_01Sep21.txt"
Epochs_dictionary['Pol'][-1] = "Polar_measurements_01Sep21.txt"
Epochs_dictionary['Coord'][-1] = "Point_List_20Sep21.txt"
Epochs_dictionary['LoS'][0] = "LoS_measurements_20Sep21.txt"
Epochs_dictionary['Pol'][0] = "Polar_measurements_20Sep21.txt"
Epochs_dictionary['Coord'][0] = "Point_List_20Sep21.txt"
Epochs_dictionary['LoS'][1] = "LoS_measurements_22Sep21.txt"
Epochs_dictionary['Pol'][1] = "Polar_measurements_22Sep21.txt"
Epochs_dictionary['Coord'][1] = "Point_List_22Sep21.txt"
Epochs_dictionary['LoS'][2] = "LoS_measurements_29Sep21.txt"
Epochs_dictionary['Pol'][2] = "Polar_measurements_29Sep21.txt"
Epochs_dictionary['Coord'][2] = "Point_List_29Sep21.txt"


"""Which Epochs gonna be used and adding the data into the code"""
if len(Which_epochs) == 1:
    LoS_Measurements_file_name = Epochs_dictionary['LoS'][Which_epochs[0]]
    Pol_Measurements_file_name = Epochs_dictionary['Pol'][Which_epochs[0]]
    Coords_file_name = Epochs_dictionary['Coord'][Which_epochs[0]]
    Res_file_name = 'Results_Epoch_' + str(Which_epochs[0]) + '.txt'
    Pico_StDev = Pico_StDev_basic * Epoch_factors[Which_epochs[0]]['Pico']
    IFM_StDev = (IFM_StDev_basic[0] * Epoch_factors[Which_epochs[0]]['IFM'],
                 IFM_StDev_basic[1] * Epoch_factors[Which_epochs[0]]['IFM'])
    ADM_StDev = ADM_StDev_basic * Epoch_factors[Which_epochs[0]]['ADM']
    Hz_StDev = Hz_StDev_basic * Epoch_factors[Which_epochs[0]]['Hz']
    V_StDev = V_StDev_basic * Epoch_factors[Which_epochs[0]]['V']
    Constraint_StDev = Constraint_StDev_basic * Epoch_factors[
                                                    Which_epochs[0]]['Con']
    Dist_StDev = (IFM_StDev[0]+ADM_StDev,IFM_StDev[1])
elif len(Which_epochs) == 2:
    LoS_Measurements_file_name = Epochs_dictionary['LoS'][Which_epochs[0]]
    LoS_Measurements_file_name_1 = Epochs_dictionary['LoS'][Which_epochs[1]]
    Pol_Measurements_file_name = Epochs_dictionary['Pol'][Which_epochs[0]]
    Pol_Measurements_file_name_1 = Epochs_dictionary['Pol'][Which_epochs[1]]
    Coords_file_name = Epochs_dictionary['Coord'][Which_epochs[0]]
    Coords_file_name_1 = Epochs_dictionary['Coord'][Which_epochs[1]]
    Res_file_name = 'Results_Epoch_' + str(Which_epochs[0]) + '.txt'
    Res_file_name_1 = 'Results_Epoch_' + str(Which_epochs[1]) + '.txt'
    Pico_StDev = Pico_StDev_basic * Epoch_factors[Which_epochs[0]]['Pico']
    IFM_StDev = (IFM_StDev_basic[0] * Epoch_factors[Which_epochs[0]]['IFM'],
                 IFM_StDev_basic[1] * Epoch_factors[Which_epochs[0]]['IFM'])
    ADM_StDev = ADM_StDev_basic * Epoch_factors[Which_epochs[0]]['ADM']
    Hz_StDev = Hz_StDev_basic * Epoch_factors[Which_epochs[0]]['Hz']
    V_StDev = V_StDev_basic * Epoch_factors[Which_epochs[0]]['V']
    Constraint_StDev = Constraint_StDev_basic * Epoch_factors[
                                                    Which_epochs[0]]['Con']
    Dist_StDev = (IFM_StDev[0]+ADM_StDev,IFM_StDev[1])
    
    Pico_StDev_E1 = Pico_StDev_basic * Epoch_factors[Which_epochs[1]]['Pico']
    IFM_StDev_E1 = (IFM_StDev_basic[0] * Epoch_factors[Which_epochs[1]]['IFM'],
                 IFM_StDev_basic[1] * Epoch_factors[Which_epochs[1]]['IFM'])
    ADM_StDev_E1 = ADM_StDev_basic * Epoch_factors[Which_epochs[1]]['ADM']
    Hz_StDev_E1 = Hz_StDev_basic * Epoch_factors[Which_epochs[1]]['Hz']
    V_StDev_E1 = V_StDev_basic * Epoch_factors[Which_epochs[1]]['V']
    Constraint_StDev_E1 = Constraint_StDev_basic * Epoch_factors[
                                                    Which_epochs[1]]['Con']
    Dist_StDev_E1 = (IFM_StDev_E1[0]+ADM_StDev_E1,IFM_StDev_E1[1])
elif len(Which_epochs) == 0:
    print('No epoch(s) were chosen to analyze. Go to config.py and '
          'change it in a variable Which_epochs')
else:
    print('Too many epochs are chosen, choose just two. Go and correct it in Which_epochs in config.py.')




Common_points = [#
									  'Girder_1', 'Girder_2', 'Girder_3', 
									  'Girder_4', 'Girder_5', 'Girder_6',
									  'Girder_7', 'Girder_8', 'Girder_9',
									  'Girder_10', 'Girder_11', 'Girder_12',
									  'Girder_13', 'Girder_14', 'Girder_15',
									  'Girder_16', 'Girder_17', 'Girder_18',
#									  'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
#									  'J', 'K', 'M', 'L', 'M'
									  ]

LSM_Excluded_measurements = {
        '-1':[],
        '0':[('Sd', 'Instrument_1', 'PQK62_8'),
             ('Sd', 'Instrument_1', 'PQK62_9')],
        '1':[],
        '2':[('Sd', 'Instrument_1', 'PQL6_9')]
        }

Lines_of_sight = {
'Hor_Left_Bottom_UP_t': ('PQK62_7','PQK62_1','Girder_5','PQL6_7',
                       'PQL6_1','Girder_13','PQK36_7','PQK36_1'),
'Hor_Left_Top_UP_t': ('PQK62_8','PQK62_2','PQL6_8','PQL6_2','PQK36_8',
                      'PQK36_2'),
'Hor_Top_Left_UP_t': ('PQK62_9','PQK62_3','PQL6_9','PQL6_3','PQK36_9',
                      'PQK36_3'),
'Hor_Top_Right_UP_t': ('PQK62_10','PQK62_4','PQL6_10',
                     'PQL6_4','PQK36_10','PQK36_4'),
'Hor_Right_Top_UP_t': ('PQK62_11','PQK62_5','PQL6_11',
                     'PQL6_5','PQK36_11','PQK36_5'),
'Hor_Right_Bottom_UP_t': ('PQK62_12','PQK62_6','Girder_6','PQL6_12',
                        'PQL6_6','Girder_14','PQK36_12','PQK36_6'),
'Hor_Left_Bottom_DN_t': ('PQK36_1','PQK36_7','Girder_13','PQL6_1',
                       'PQL6_7','Girder_5','PQK62_1','PQK62_7'),
'Hor_Left_Top_DN_t': ('PQK36_2','PQK36_8','PQL6_2','PQL6_8','PQK62_2',
                      'PQK62_8'),
'Hor_Top_Left_DN_t': ('PQK36_3','PQK36_9','PQL6_3','PQL6_9','PQK62_3',
                      'PQK62_9'),
'Hor_Top_Right_DN_t': ('PQK36_4','PQK36_10','PQL6_4',
                     'PQL6_10','PQK62_4','PQK62_10'),
'Hor_Right_Top_DN_t': ('PQK36_5','PQK36_11','PQL6_5',
                     'PQL6_11','PQK62_5','PQK62_11'),
'Hor_Right_Bottom_DN_t':('PQK36_6','PQK36_12','Girder_14','PQL6_6',
                       'PQL6_12','Girder_6','PQK62_6','PQK62_12'),
'Pico_Left_DN': ('PQL6_1','Girder_11'),
'Pico_Left_UP': ('PQL6_7','Girder_7'),
'Pico_Right_UP': ('PQL6_12','Girder_8'),
'Pico_Right_DN': ('PQL6_6','Girder_12'),
'Ver_Left_DN_DN_t': ('PQK36_2','PQK36_1','Girder_17'),
'Ver_Left_DN_UP_t': ('PQK36_8','PQK36_7','Girder_15'),
'Ver_Left_MD_DN_t': ('PQL6_2','PQL6_1','Girder_11'),
'Ver_Left_MD_UP_t': ('PQL6_8','PQL6_7','Girder_7'),
'Ver_Left_UP_DN_t': ('PQK62_2','PQK62_1','Girder_3'),
'Ver_Left_UP_UP_t': ('PQK62_8','PQK62_7','Girder_1'),
'Ver_Right_DN_DN_t': ('PQK36_5','PQK36_6','Girder_18'),
'Ver_Right_DN_UP_t': ('PQK36_11','PQK36_12','Girder_16'),
'Ver_Right_MD_DN_t': ('PQL6_5','PQL6_6','Girder_12'),
'Ver_Right_MD_UP_t': ('PQL6_11','PQL6_12','Girder_8'),
'Ver_Right_UP_DN_t': ('PQK62_5','PQK62_6','Girder_4'),
'Ver_Right_UP_UP_t': ('PQK62_11','PQK62_12','Girder_2'),
'Diag_Left_DN_t': ('PQK36_8','PQL6_1','Girder_9'),
'Diag_Left_UP_t': ('PQK62_2','PQL6_7','Girder_9'),
'Diag_Right_DN_t': ('PQK36_11','PQL6_6','Girder_10'),
'Diag_Right_UP_t': ('PQK62_5','PQL6_12','Girder_10'),
'Diag_Top_DN_t': ('PQK36_10','PQL6_3'),
'Diag_Top_UP_t': ('PQK62_4','PQL6_9'),
'Hor_Left_Bottom_UP_b': ('PQK36_1','PQK36_7','Girder_13','PQL6_1','PQL6_7',
                         'Girder_5','PQK62_1','PQK62_7'),
'Hor_Left_Top_UP_b': ('PQK36_2','PQK36_8','PQL6_2','PQL6_8','PQK62_2',
                      'PQK62_8'),  
'Hor_Top_Left_UP_b': ('PQK36_3','PQK36_9','PQL6_3','PQL6_9','PQK62_3',
                      'PQK62_9'),
'Hor_Top_Right_UP_b': ('PQK36_4','PQK36_10','PQL6_4','PQL6_10','PQK62_4',
                       'PQK62_10'),
'Hor_Right_Top_UP_b': ('PQK36_5','PQK36_11','PQL6_5','PQL6_11','PQK62_5',
                       'PQK62_11'),
'Hor_Right_Bottom_UP_b': ('PQK36_6','PQK36_12','Girder_14','PQL6_6','PQL6_12',
                          'Girder_6','PQK62_6','PQK62_12'),
'Hor_Left_Bottom_DN_b': ('PQK62_7','PQK62_1','Girder_5','PQL6_7','PQL6_1',
                         'Girder_13','PQK36_7','PQK36_1'),                        
'Hor_Left_Top_DN_b': ('PQK62_8','PQK62_2','PQL6_8','PQL6_2','PQK36_8',
                      'PQK36_2'),
'Hor_Top_Left_DN_b': ('PQK62_9','PQK62_3','PQL6_9','PQL6_3','PQK36_9',
                      'PQK36_3'),
'Hor_Top_Right_DN_b': ('PQK62_10','PQK62_4','PQL6_10','PQL6_4','PQK36_10',
                       'PQK36_4'),
'Hor_Right_Top_DN_b': ('PQK62_11','PQK62_5','PQL6_11','PQL6_5','PQK36_11',
                       'PQK36_5'),
'Hor_Right_Bottom_DN_b':('PQK62_12','PQK62_6','Girder_6','PQL6_12','PQL6_6',
                         'Girder_14','PQK36_12','PQK36_6'),
'Pico_Left_DN_b': ('Girder_11','PQL6_1'),
'Pico_Left_UP_b': ('Girder_7','PQL6_7'),
'Pico_Right_UP_b': ('Girder_8','PQL6_12'),
'Pico_Right_DN_b': ('Girder_12','PQL6_6'),
'Ver_Left_DN_DN_b': ('Girder_17','PQK36_1','PQK36_2'),
'Ver_Left_DN_UP_b': ('Girder_15','PQK36_7','PQK36_8'),
'Ver_Left_MD_DN_b': ('Girder_11','PQL6_1','PQL6_2'),
'Ver_Left_MD_UP_b': ('Girder_7','PQL6_7','PQL6_8'),
'Ver_Left_UP_DN_b': ('Girder_3','PQK62_1','PQK62_2'),
'Ver_Left_UP_UP_b': ('Girder_1','PQK62_7','PQK62_8'),
'Ver_Right_DN_DN_b': ('Girder_18','PQK36_6','PQK36_5'),
'Ver_Right_DN_UP_b': ('Girder_16','PQK36_12','PQK36_11'),
'Ver_Right_MD_DN_b': ('Girder_12','PQL6_6','PQL6_5'),
'Ver_Right_MD_UP_b': ('Girder_8','PQL6_12','PQL6_11'),
'Ver_Right_UP_DN_b': ('Girder_4','PQK62_6','PQK62_5'),
'Ver_Right_UP_UP_b': ('Girder_2','PQK62_12','PQK62_11'),
'Diag_Left_DN_b': ('Girder_9','PQL6_1','PQK36_8'),
'Diag_Left_UP_b': ('Girder_9','PQL6_7','PQK62_2'),
'Diag_Right_DN_b': ('Girder_10','PQL6_6','PQK36_11'),
'Diag_Right_UP_b': ('Girder_10','PQL6_12','PQK62_5'),
'Diag_Top_DN_b': ('PQL6_3','PQK36_10'),
'Diag_Top_UP_b': ('PQL6_9','PQK62_4')
}

