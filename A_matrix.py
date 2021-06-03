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
Pol_measurements_cart = MC.Pol_measurements_cart
StDevs_IFM_measurements = MC.StDevs_IFM_measurements
measured_distances_in_lines = MC.measured_distances_in_lines
sorted_measured_points_in_lines = MC.sorted_measured_points_in_lines

if MC.Two_epochs:
    LoS_measurements_E1 = MC.LoS_measurements_E1
    Pol_measurements_E1 = MC.Pol_measurements_E1
    Pol_measurements_cart_E1 = MC.Pol_measurements_cart_E1
    StDevs_IFM_measurements_E1 = MC.StDevs_IFM_measurements_E1
    measured_distances_in_lines_E1 = MC.measured_distances_in_lines_E1
    sorted_measured_points_in_lines_E1 = MC.sorted_measured_points_in_lines_E1

#A_matrix = np.array([[], []], np.double)
#m = 5
#n = 6
#A_matrix = np.zeroes([m,n])
#print(A_matrix.shape)

count_IFM_measurements = sum([len(v)+1 for k, v in\
                                         measured_distances_in_lines.items()])\
                         - len(measured_distances_in_lines)
count_Pico_measurements = 4
count_Pol_measurements = (sum([len(v) for k, v in Pol_measurements.items()]))*3
                        
count_all_observations = count_IFM_measurements + count_Pico_measurements \
                         + count_Pol_measurements
unknowns = []
for instrument in Pol_measurements:
    keys = list(Pol_measurements[instrument].keys())
    unknowns = list(set(unknowns + keys))
unknowns.sort()
print(unknowns)

print(count_all_observations)