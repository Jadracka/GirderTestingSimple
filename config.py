# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:54:43 2021

@author: jbarker
"""

#import math as m
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
Which_epochs = (10,11) #the comma must stay, otherwise the variable will be int
Names_of_magnets = ['PQK36','PQL6','PQK62']
Instruments_6DoF = True

Pico_StDev_basic = 0.00000001 #mm, PicoScale Standard Deviation
IFM_StDev_basic = (0.0004,0.00015) #mm, LT IFM, based on Leica's white paper
ADM_StDev_basic = 0.010 #mm, from 2-5m based on Leica's typical errors
Hz_StDev_basic = 0.15 #mgon same like Leica TM50
V_StDev_basic = 0.15 #mgon same like Leica TM50
Constraint_StDev_basic = 0.000001#0.000000055136545 #mm


Max_diff_from_line = 1.7 #mm - maximum distance from line, report the excess
LSM_Threshold = 1e-14
LSM_Max_iterations = 15
Epsilon = 0.001
LSM_incl_IFM = True #Include IFM measurements in LSM?
LSM_incl_Cons = True #Include Constraints in LSM?


Sigma_0 = 0.000001



# Epoch specific factors
Epoch_factors = {-2:{}, -1:{}, 0:{}, 1:{}, 2:{}, 3:{}, 4:{}, 5:{}, 6:{}, 7:{},
						    8:{}, 9:{}, 10:{}, 11:{}}

Epoch_factors[-2]['Pico'] = 1
Epoch_factors[-2]['IFM'] = 1
Epoch_factors[-2]['ADM'] = 1
Epoch_factors[-2]['Hz'] = 1
Epoch_factors[-2]['V'] = 1
Epoch_factors[-2]['Con'] = 1

Epoch_factors[-1]['Pico'] = 1
Epoch_factors[-1]['IFM'] = 1
Epoch_factors[-1]['ADM'] = 1
Epoch_factors[-1]['Hz'] = 1
Epoch_factors[-1]['V'] = 1
Epoch_factors[-1]['Con'] = 1

Epoch_factors[0]['Pico'] = 1
Epoch_factors[0]['IFM'] = 5.7532
Epoch_factors[0]['ADM'] = 0.4542
Epoch_factors[0]['Hz'] = 0.738
Epoch_factors[0]['V'] = 2.6179
Epoch_factors[0]['Con'] = 1

Epoch_factors[1]['Pico'] = 1
Epoch_factors[1]['IFM'] = 6.3832
Epoch_factors[1]['ADM'] = 0.7201
Epoch_factors[1]['Hz'] = 0.9623
Epoch_factors[1]['V'] = 10
Epoch_factors[1]['Con'] = 1

Epoch_factors[2]['Pico'] = 1
Epoch_factors[2]['IFM'] = 5.5464
Epoch_factors[2]['ADM'] = 0.3309
Epoch_factors[2]['Hz'] = 1.2660
Epoch_factors[2]['V'] = 2.4926
Epoch_factors[2]['Con'] = 1

Epoch_factors[3]['Pico'] = 1
Epoch_factors[3]['IFM'] = 11.49
Epoch_factors[3]['ADM'] = 0.1
Epoch_factors[3]['Hz'] = 0.99
Epoch_factors[3]['V'] = 2.56
Epoch_factors[3]['Con'] = 1

Epoch_factors[4]['Pico'] = 1
Epoch_factors[4]['IFM'] = 6.246
Epoch_factors[4]['ADM'] = 0.0864
Epoch_factors[4]['Hz'] = 0.71
Epoch_factors[4]['V'] = 2.177
Epoch_factors[4]['Con'] = 1

Epoch_factors[5]['Pico'] = 1
Epoch_factors[5]['IFM'] = 6.8
Epoch_factors[5]['ADM'] = 0.185
Epoch_factors[5]['Hz'] = 0.5
Epoch_factors[5]['V'] = 2.177
Epoch_factors[5]['Con'] = 1

Epoch_factors[6]['Pico'] = 1
Epoch_factors[6]['IFM'] = 6.3
Epoch_factors[6]['ADM'] = 0.17
Epoch_factors[6]['Hz'] = 0.71
Epoch_factors[6]['V'] = 2.25
Epoch_factors[6]['Con'] = 1

Epoch_factors[7]['Pico'] = 1
Epoch_factors[7]['IFM'] = 9.3
Epoch_factors[7]['ADM'] = 0.03
Epoch_factors[7]['Hz'] = 0.75
Epoch_factors[7]['V'] = 2.53
Epoch_factors[7]['Con'] = 1

Epoch_factors[8]['Pico'] = 1
Epoch_factors[8]['IFM'] = 8.95
Epoch_factors[8]['ADM'] = 0.03
Epoch_factors[8]['Hz'] = 0.87
Epoch_factors[8]['V'] = 2.53
Epoch_factors[8]['Con'] = 1

Epoch_factors[9]['Pico'] = 1
Epoch_factors[9]['IFM'] = 7.05
Epoch_factors[9]['ADM'] = 0.075
Epoch_factors[9]['Hz'] = 0.86
Epoch_factors[9]['V'] = 2.955
Epoch_factors[9]['Con'] = 1

Epoch_factors[10]['Pico'] = 1
Epoch_factors[10]['IFM'] = 14.8
Epoch_factors[10]['ADM'] = 0.01
Epoch_factors[10]['Hz'] = 0.94
Epoch_factors[10]['V'] = 1.73
Epoch_factors[10]['Con'] = 1

Epoch_factors[11]['Pico'] = 1
Epoch_factors[11]['IFM'] = 6.53
Epoch_factors[11]['ADM'] = 0.075
Epoch_factors[11]['Hz'] = 0.94
Epoch_factors[11]['V'] = 2.94
Epoch_factors[11]['Con'] = 1

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
Epochs_dictionary['LoS'][-2] = "LoS_measurements_01Sep21.txt"
Epochs_dictionary['Pol'][-2] = "Testing_pol_meas_simplest.txt"
Epochs_dictionary['Coord'][-2] = "Testing_coords_simplest.txt"
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
Epochs_dictionary['LoS'][3] = "LoS_measurements_21Feb22.txt"
Epochs_dictionary['Pol'][3] = "Polar_measurements_21Feb22.txt"
Epochs_dictionary['Coord'][3] = "Point_List_21Feb22.txt"
Epochs_dictionary['LoS'][4] = "LoS_measurements_22Feb22.txt"
Epochs_dictionary['Pol'][4] = "Polar_measurements_22Feb22.txt"
Epochs_dictionary['Coord'][4] = "Point_List_22Feb22.txt"
Epochs_dictionary['LoS'][5] = "LoS_measurements_24Feb22.txt"
Epochs_dictionary['Pol'][5] = "Polar_measurements_24Feb22.txt"
Epochs_dictionary['Coord'][5] = "Point_List_24Feb22.txt"
Epochs_dictionary['LoS'][6] = "LoS_measurements_1Mar22.txt"
Epochs_dictionary['Pol'][6] = "Polar_measurements_1Mar22.txt"
Epochs_dictionary['Coord'][6] = "Point_List_1Mar22.txt"
Epochs_dictionary['LoS'][7] = "LoS_measurements_2Mar22.txt"
Epochs_dictionary['Pol'][7] = "Polar_measurements_2Mar22.txt"
Epochs_dictionary['Coord'][7] = "Point_List_2Mar22.txt"
Epochs_dictionary['LoS'][8] = "LoS_measurements_7Mar22.txt"
Epochs_dictionary['Pol'][8] = "Polar_measurements_7Mar22.txt"
Epochs_dictionary['Coord'][8] = "Point_List_7Mar22.txt"
Epochs_dictionary['LoS'][9] = "LoS_measurements_22Mar22.txt"
Epochs_dictionary['Pol'][9] = "Polar_measurements_22Mar22.txt"
Epochs_dictionary['Coord'][9] = "Point_List_22Mar22.txt"
Epochs_dictionary['LoS'][10] = "LoS_measurements_25Mar22.txt"
Epochs_dictionary['Pol'][10] = "Polar_measurements_25Mar22.txt"
Epochs_dictionary['Coord'][10] = "Point_List_25Mar22.txt"
Epochs_dictionary['LoS'][11] = "LoS_measurements_29Mar22.txt"
Epochs_dictionary['Pol'][11] = "Polar_measurements_29Mar22.txt"
Epochs_dictionary['Coord'][11] = "Point_List_29Mar22.txt"

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




Common_points = ['PQK36_1','PQK36_2','PQK36_3','PQK36_4','PQK36_5','PQK36_6',
						   'PQK36_7','PQK36_8','PQK36_9','PQK36_10','PQK36_11','PQK36_12',
						   'PQK62_1','PQK62_2','PQK62_3','PQK62_4','PQK62_5','PQK62_6',
						   'PQK62_7','PQK62_8','PQK62_9','PQK62_10','PQK62_11','PQK62_12',
						   'PQL6_1','PQL6_2','PQL6_3','PQL6_4','PQL6_5','PQL6_6'
						   'PQL6_7', 'PQL6_8','PQL6_9','PQL6_10','PQL6_11','PQL6_12']

LSM_Excluded_measurements = {
        '-2':[],
        '-1':[],
        '0':[#('Sd', 'Instrument_1', 'PQK62_8'),
             #('Sd', 'Instrument_1', 'PQK62_9')
             ],
        '1':[],
        '2':[#('Sd', 'Instrument_1', 'PQL6_9')
            ],
        '3':[],
        '4':[],
        '5':[],
        '6':[],
        '7':[],
        '8':[],
        '9':[],
        '10':[],
        '11':[]
        }

Instruments_LoS = {
'Hor_Left_Bottom_UP_t': 'Instrument_1',
'Hor_Left_Top_UP_t': 'Instrument_1',
'Hor_Top_Left_UP_t': 'Instrument_1',
'Hor_Top_Right_UP_t': 'Instrument_0',
'Hor_Right_Top_UP_t': 'Instrument_0',
'Hor_Right_Bottom_UP_t': 'Instrument_0',
'Hor_Left_Bottom_DN_t': 'Instrument_1',
'Hor_Left_Top_DN_t': 'Instrument_1',
'Hor_Top_Left_DN_t': 'Instrument_1',
'Hor_Top_Right_DN_t': 'Instrument_1',
'Hor_Right_Top_DN_t': 'Instrument_0',
'Hor_Right_Bottom_DN_t': 'Instrument_0',
'Ver_Left_DN_DN_t': 'Instrument_1',
'Ver_Left_DN_UP_t': 'Instrument_1',
'Ver_Left_MD_DN_t': 'Instrument_1',
'Ver_Left_MD_UP_t': 'Instrument_1',
'Ver_Left_UP_DN_t': 'Instrument_1',
'Ver_Left_UP_UP_t': 'Instrument_1',
'Ver_Right_DN_DN_t': 'Instrument_0',
'Ver_Right_DN_UP_t': 'Instrument_0',
'Ver_Right_MD_DN_t': 'Instrument_0',
'Ver_Right_MD_UP_t': 'Instrument_0',
'Ver_Right_UP_DN_t':'Instrument_0',
'Ver_Right_UP_UP_t': 'Instrument_0',
'Diag_Left_DN_t': 'Instrument_1',
'Diag_Left_UP_t': 'Instrument_1',
'Diag_Right_DN_t': 'Instrument_1',
'Diag_Right_UP_t': 'Instrument_1',
'Diag_Top_DN_t': 'Instrument_0',
'Diag_Top_UP_t': 'Instrument_0',
'Hor_Left_Bottom_UP_b':'Instrument_1',
'Hor_Left_Top_UP_b': 'Instrument_1',
'Hor_Top_Left_UP_b': 'Instrument_1',
'Hor_Top_Right_UP_b': 'Instrument_0',
'Hor_Right_Top_UP_b': 'Instrument_0',
'Hor_Right_Bottom_UP_b': 'Instrument_0',
'Hor_Left_Bottom_DN_b': 'Instrument_1',
'Hor_Left_Top_DN_b': 'Instrument_1',
'Hor_Top_Left_DN_b': 'Instrument_1',
'Hor_Top_Right_DN_b': 'Instrument_0',
'Hor_Right_Top_DN_b': 'Instrument_0',
'Hor_Right_Bottom_DN_b': 'Instrument_0',
'Ver_Left_DN_DN_b': 'Instrument_1',
'Ver_Left_DN_UP_b': 'Instrument_1',
'Ver_Left_MD_DN_b': 'Instrument_1',
'Ver_Left_MD_UP_b': 'Instrument_1',
'Ver_Left_UP_DN_b': 'Instrument_1',
'Ver_Left_UP_UP_b': 'Instrument_1',
'Ver_Right_DN_DN_b': 'Instrument_0',
'Ver_Right_DN_UP_b': 'Instrument_0',
'Ver_Right_MD_DN_b': 'Instrument_0',
'Ver_Right_MD_UP_b': 'Instrument_0',
'Ver_Right_UP_DN_b': 'Instrument_0',
'Ver_Right_UP_UP_b': 'Instrument_0',
'Diag_Left_DN_b': 'Instrument_1',
'Diag_Left_UP_b': 'Instrument_1',
'Diag_Right_DN_b': 'Instrument_0',
'Diag_Right_UP_b': 'Instrument_0',
'Diag_Top_DN_b': 'Instrument_0',
'Diag_Top_UP_b': 'Instrument_0'
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

# =============================================================================
# FIDS data
# =============================================================================

FIDS = {"PQL6_1":(-0.32987,0.31877,-0.26081),
			 "PQL6_2":(-0.330474,0.318568,0.259824),
			 "PQL6_3":(-0.330385,0.189728,0.380554),
			 "PQL6_4":(-0.330511,-0.190224,0.380736),
			 "PQL6_5":(-0.330219,-0.318912,0.259614),
			 "PQL6_6":(-0.33061,-0.319309,-0.259964),
			 "PQL6_7":(0.329582,0.318921,-0.260355),
			 "PQL6_8":(0.329156,0.318387,0.259552),
			 "PQL6_9":(0.329258,0.189746,0.380537),
			 "PQL6_10":(0.329284,-0.19037,0.380672),
			 "PQL6_11":(0.329272,-0.318666,0.258858),
			 "PQL6_12":(0.329065,-0.319084,-0.260201),
			 "PQL6_Aus":(0.25,0.0,0.0),
			 "PQL6_Ein":(-0.25,0.0,0.0),
			 "PQL6_Oben":(0.0,0.0,0.25),
			 "PQL6_Mitte":(0.0,0.0,0.0),
			 "PQK62_1":(-0.189870,0.319170,-0.26036),
			 "PQK62_2":(-0.190772,0.319193,0.260388),
			 "PQK62_3":(-0.190201,0.189704,0.381166),
			 "PQK62_4":(-0.189577,-0.189666,0.381471),
			 "PQK62_5":(-0.190273,-0.319128,0.25966),
			 "PQK62_6":(-0.189251,-0.318936,-0.259724),
			 "PQK62_7":(0.190031,0.31891,-0.259705),
			 "PQK62_8":(0.189403,0.31936,0.260542),
			 "PQK62_9":(0.190539,0.190083,0.38121),
			 "PQK62_10":(0.190435,-0.190107,0.38098),
			 "PQK62_11":(0.190146,-0.319081,0.260108),
			 "PQK62_12":(0.190451,-0.319175,-0.260208),
			 "PQK62_Aus":(0.25,0.0,0.0),
			 "PQK62_Ein":(-0.25,0.0,0.0),
			 "PQK62_Oben":(0.0,0.0,0.25),
			 "PQK62_Mitte":(0.0,0.0,0.0),
			 "PQK36_1":(-0.18953,0.31912,-0.2602),
			 "PQK36_2":(-0.190255,0.318954,0.259534),
			 "PQK36_3":(-0.190804,0.190681,0.381003),
			 "PQK36_4":(-0.19024,-0.189815,0.380975),
			 "PQK36_5":(-0.189777,-0.319055,0.259512),
			 "PQK36_6":(-0.189818,-0.318791,-0.259625),
			 "PQK36_7":(0.190259,0.319144,-0.259908),
			 "PQK36_8":(0.189735,0.319041,0.259833),
			 "PQK36_9":(0.189238,0.190138,0.381262),
			 "PQK36_10":(0.189378,-0.190416,0.380912),
			 "PQK36_11":(0.190297,-0.31924,0.25954),
			 "PQK36_12":(0.190255,-0.318888,-0.259811),
			 "PQK36_Aus":(0.25,0.0,0.0),
			 "PQK36_Ein":(-0.25,0.0,0.0),
			 "PQK36_Oben":(0.0,0.0,0.250),
			 "PQK36_Mitte":(0.0,0.0,0.0)
}