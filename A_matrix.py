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
A_matrixHR = {}

L_vector = np.array([])
L_vectorHR = []

del count_all_observations

Aproximates = fc.merge_measured_coordinates(Transformed_Pol_measurements)

X_vector = fc.filling_X(
                      Aproximates, unknowns, count_unknowns, count_instruments)

# Filling A and L with IFM measurements
for line in measured_distances_in_lines:
    # meas_i is index of the measurand in line, so I can pair it with the
    # Point IDs from sorted points in lines, first measured point is on the
    # same index as the distance and second point is on index + 1
    for meas_i,distance in enumerate(measured_distances_in_lines[line]):
        # I need index of the measurement in L vector, so I can put it on
        # correct corresponding row in First plan matrix A
        L_i = len(L_vector)
        L_vector = np.append(L_vector, distance)
        """Now figuring the index of points in the unknowns list, so I know
           which column of the A matrix. The original index is multiplied by 3
           because unknowns are names of points (and instrument orientations)
           not a complete list of unknowns = Point X,Y and Z => *3         """
        PointFrom = sorted_measured_points_in_lines[line][meas_i]
        PointTo = sorted_measured_points_in_lines[line][meas_i+1]
        PointFrom_i = 3*unknowns.index(PointFrom) 
        PointTo_i = 3*unknowns.index(PointTo)
        dX,dY,dZ = fc.ParD_sd(Aproximates[PointTo],Aproximates[PointFrom])
        A_matrix[L_i,PointTo_i:PointTo_i+3] = dX,dY,dZ
        # for point "From" the partial derivatives change sign
        A_matrix[L_i,PointFrom_i:PointFrom_i+3] = -dX,-dY,-dZ
        # Documenting in human readible format what are the A and L elements
        L_vectorHR.append((line, PointFrom, PointTo))
        A_matrixHR[(L_i,PointFrom_i)] = ['dX', line, PointFrom, PointTo]
        A_matrixHR[(L_i,PointFrom_i+1)] = ['dY', line, PointFrom, PointTo]
        A_matrixHR[(L_i,PointFrom_i+2)] = ['dZ', line, PointFrom, PointTo]
        A_matrixHR[(L_i,PointTo_i)] = ['dX', line, PointTo, PointFrom]
        A_matrixHR[(L_i,PointTo_i+1)] = ['dY', line, PointTo, PointFrom]
        A_matrixHR[(L_i,PointTo_i+2)] = ['dZ', line, PointTo, PointFrom]
del PointFrom,PointTo, PointFrom_i,PointTo_i,dX,dY,dZ,L_i

for instrument,points in Pol_measurements.items():
    for i,point in enumerate(points):
        print(point)
        
    
    pass

# For HZ angles
#if ((count_unknowns-2*count_instruments) <= PointFrom_i <= (
#                                             count_unknowns-count_instruments)
#or (count_unknowns-2*count_instruments) <= PointTo_i <= (
#                                            count_unknowns-count_instruments)):


print('End of A_matrix code')