# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:21:32 2021

@author: jbarker
"""

#import scipy as sp
import numpy as np
#import sys
#import string
#import math as m
import config as cg
import functions as fc
import MainCode as MC

#from operator import itemgetter
#from collections import namedtuple

   
Nominal_coords = MC.Nominal_coords   
LoS_measurements = MC.LoS_measurements
Pol_measurements = MC.Pol_measurements
Transformed_Pol_measurements = MC.Transformed_Pol_measurements
StDevs_IFM_measurements = MC.StDevs_IFM_measurements
measured_distances_in_lines = MC.measured_distances_in_lines
sorted_measured_points_in_lines = MC.sorted_measured_points_in_lines
Two_epochs = MC.Two_epochs

if Two_epochs:
    LoS_measurements_E1 = MC.LoS_measurements_E1
    Pol_measurements_E1 = MC.Pol_measurements_E1
    Transformed_Pol_measurements_E1 = MC.Transformed_Pol_measurements_E1
    StDevs_IFM_measurements_E1 = MC.StDevs_IFM_measurements_E1
    measured_distances_in_lines_E1 = MC.measured_distances_in_lines_E1
    sorted_measured_points_in_lines_E1 = MC.sorted_measured_points_in_lines_E1

count_IFM_measurements = sum([len(v) for k, v in\
                                         measured_distances_in_lines.items()])
count_Pico_measurements = 4
count_Pol_measurements = (sum([len(v) for k, v in Pol_measurements.items()]))*3
                        
count_all_observations = count_IFM_measurements + count_Pico_measurements \
                         + count_Pol_measurements
del count_IFM_measurements, count_Pico_measurements, count_Pol_measurements

unknowns,count_unknowns,count_instruments = fc.find_unknowns(
                                                  Transformed_Pol_measurements)

if Two_epochs:
    unknowns_E1 = fc.find_unknowns(Pol_measurements_E1)
    if cg.Print_epoch_checks:
        if unknowns==unknowns_E1:
            print('Unknowns in both epoch match')
        else:
            print("Unknowns in epochs don't match")

A_matrix = np.zeros([count_all_observations,count_unknowns])
L_vector = np.zeros(count_all_observations)
X_vector = np.zeros(count_unknowns)
del count_all_observations

Aproximates = fc.merge_measured_coordinates(
                        Transformed_Pol_measurements['Instrument_0'],
                        Transformed_Pol_measurements['Instrument_1'])

# Filling the X vector
for item in unknowns[:-count_instruments]:
#    print(item)
    pass
    

# IFM measurements
for line in sorted_measured_points_in_lines:
#    print(line)
    pass
    
print('End of A_matrix code')
    