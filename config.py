# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:54:43 2021

@author: jbarker
"""
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
Which_epochs = (0,1) #the comma must stay, otherwise the variable will be int
Using_nominal_compare = True
Line_differences_checking = False
Print_FIDs = False
Print_typos = True
Print_real2nominal_checks = True
IFM_StDev_sys = (0.0004,0.00015) #mm, +- 0,2um based on Leica's white paper
ADM_StDev = 0.010 #mm, from 2-5m based on Leica's typical errors
Ang_StDev = (0.015,0.003) #mm, needs to be calculated for each angle separately



Dist_StDev = (IFM_StDev_sys[0]+ADM_StDev,IFM_StDev_sys[1])

"""Nominal CAD coordinates of Fiducials
Standard SA format with spaces as delimiters and no comments:
Point name 'space' X [mm] 'space' Y [mm] 'space' Z [mm] 
"""
Coords_file_name = "BetterCoords.txt"


"""Measured Lines
Standard SA format with spaces as delimiters and no comments:
Group aka Line name 'space' Point name 'space' Sd [mm] 'space' Hz [gon] 'space' V [gon] 
"""
Epochs_dictionary = {'LoS':{},'Pol':{}}
Epochs_dictionary['LoS'][0] = "Testing_measurements_Epoch0.txt"
Epochs_dictionary['LoS'][1] = "Testing_measurements_Epoch1.txt"
Epochs_dictionary['Pol'][0] = "Polar_measurements_0.txt"
Epochs_dictionary['Pol'][1] = "Polar_measurements_1.txt"

"""Which Epochs gonna be used and adding the data into the code"""
if len(Which_epochs) == 1:
    LoS_Measurements_file_name = Epochs_dictionary['LoS'][Which_epochs[0]]
    Pol_Measurements_file_name = Epochs_dictionary['Pol'][Which_epochs[0]]
elif len(Which_epochs) == 2:
    LoS_Measurements_file_name = Epochs_dictionary['LoS'][Which_epochs[0]]
    LoS_Measurements_file_name_1 = Epochs_dictionary['LoS'][Which_epochs[1]]
    Pol_Measurements_file_name = Epochs_dictionary['Pol'][Which_epochs[0]]
    Pol_Measurements_file_name_1 = Epochs_dictionary['Pol'][Which_epochs[1]]
elif len(Which_epochs) == 0:
    print('No epoch(s) were chosen to analyse. Go to config.py and change it in a variable Which_epochs')
else:
    print('Too many epochs are chosen, choose just two. Go and correct it in Which_epochs in config.py.')
    
Lines_of_sight = {
'Hor_Left_Bottom_UP': ('PQK62_7','PQK62_1','Girder_5','PQL6_7','PQL6_1','Girder_13','PQK36_7','PQK36_1'),
'Hor_Left_Top_UP': ('PQK62_8','PQK62_2','PQL6_8','PQL6_2','PQK36_8','PQK36_2'),
'Hor_Top_Left_UP': ('PQK62_9','PQK62_3','PQL6_9','PQL6_3','PQK36_9','PQK36_3'),
'Hor_Top_Right_UP': ('PQK62_10','PQK62_4','PQL6_10','PQL6_4','PQK36_10','PQK36_4'),
'Hor_Right_Top_UP': ('PQK62_11','PQK62_5','PQL6_11','PQL6_5','PQK36_11','PQK36_5'),
'Hor_Right_Bottom_UP': ('PQK62_12','PQK62_6','Girder_6','PQL6_12','PQL6_6','Girder_14','PQK36_12','PQK36_6'),
'Hor_Left_Bottom_DN': ('PQK36_1','PQK36_7','Girder_13','PQL6_1','PQL6_7','Girder_5','PQK62_1','PQK62_7'),
'Hor_Left_Top_DN': ('PQK36_2','PQK36_8','PQL6_2','PQL6_8','PQK62_2','PQK62_8'),
'Hor_Top_Left_DN': ('PQK36_3','PQK36_9','PQL6_3','PQL6_9','PQK62_3','PQK62_9'),
'Hor_Top_Right_DN': ('PQK36_4','PQK36_10','PQL6_4','PQL6_10','PQK62_4','PQK62_10'),
'Hor_Right_Top_DN': ('PQK36_5','PQK36_11','PQL6_5','PQL6_11','PQK62_5','PQK62_11'),
'Hor_Right_Bottom_DN':('PQK36_6','PQK36_12','Girder_14','PQL6_6','PQL6_12','Girder_6','PQK62_6','PQK62_12'),
'IFM_Left_DN': ('PQL6_1','Girder_11'),
'IFM_Left_UP': ('PQL6_7','Girder_7'),
'IFM_Right_UP': ('PQL6_12','Girder_8'),
'IFM_Right_DN': ('PQL6_6','Girder_12'),
'Ver_Left_DN_DN': ('PQK36_2','PQK36_1','Girder_17'),
'Ver_Left_DN_UP': ('PQK36_8','PQK36_7','Girder_15'),
'Ver_Left_MD_DN': ('PQL6_2','PQL6_1','Girder_11'),
'Ver_Left_MD_UP': ('PQL6_8','PQL6_7','Girder_7'),
'Ver_Left_UP_DN': ('PQK62_2','PQK62_1','Girder_3'),
'Ver_Left_UP_UP': ('PQK62_8','PQK62_7','Girder_1'),
'Ver_Right_DN_DN': ('PQK36_5','PQK36_6','Girder_18'),
'Ver_Right_DN_UP': ('PQK36_11','PQK36_12','Girder_16'),
'Ver_Right_MD_DN': ('PQL6_5','PQL6_6','Girder_12'),
'Ver_Right_MD_UP': ('PQL6_11','PQL6_12','Girder_8'),
'Ver_Right_UP_DN': ('PQK62_5','PQK62_6','Girder_4'),
'Ver_Right_UP_UP': ('PQK62_11','PQK62_12','Girder_2'),
'Diag_Left_DN': ('PQK36_8','PQL6_1','Girder_9'),
'Diag_Left_UP': ('PQK62_2','PQL6_7','Girder_9'),
'Diag_Right_DN': ('PQK36_11','PQL6_6','Girder_10'),
'Diag_Right_UP': ('PQK62_5','PQL6_12','Girder_10'),
'Diag_Top_DN': ('PQK36_10','PQL6_3'),
'Diag_Top_UP': ('PQK62_4','PQL6_9')
}

