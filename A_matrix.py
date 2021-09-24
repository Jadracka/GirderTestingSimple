# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:21:32 2021

@author: jbarker
"""

#import scipy as sp
import numpy as np
from numpy.linalg import inv
#import sys
#import string
#import config as cg
import math as m
import functions as fc
import MainCode as MC
import config as cg

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

#if Two_epochs:
#    LoS_measurements_E1 = MC.LoS_measurements_E1
#    Pol_measurements_E1 = MC.Pol_measurements_E1
#    Transformed_Pol_measurements_E1 = MC.Transformed_Pol_measurements_E1
#    StDevs_IFM_measurements_E1 = MC.StDevs_IFM_measurements_E1
#    measured_distances_in_lines_E1 = MC.measured_distances_in_lines_E1
#    sorted_measured_points_in_lines_E1 = MC.sorted_measured_points_in_lines_E1

count_IFM_measurements = sum([len(v) for k, v in\
                                         measured_distances_in_lines.items()])
count_Pico_measurements = 0
count_Pol_measurements = (sum([len(v) for k, v in Pol_measurements.items()]))
                        
count_all_observations = count_Pico_measurements \
                         + 3*count_Pol_measurements + count_IFM_measurements
del count_IFM_measurements, count_Pico_measurements

unknowns,count_unknowns, instruments, count_instruments = fc.find_unknowns(
                                                  Transformed_Pol_measurements)
Aproximates = fc.merge_measured_coordinates(Transformed_Pol_measurements)

for instrument in instruments:
    Aproximates['Ori_' + instrument] = 0


# Creating a list of tuple pairs filled with PointIDs of point to point pairs
# on magnets. This will be used to fill the A matrix with pseudomeasurements
# instead of constraints.
Combinations_for_constraints = []
for magnet in cg.Names_of_magnets:
    points_on_magnet = [key for key in Aproximates.keys() if magnet in key]
    magnet_combinations = [(PointA, PointB) for index, PointA in enumerate(
            points_on_magnet) for PointB in points_on_magnet[index + 1:]]
    Combinations_for_constraints.extend(magnet_combinations)
count_constraints = len(Combinations_for_constraints)

#if Two_epochs:
#    unknowns_E1 = fc.find_unknowns(Pol_measurements_E1)
#   if cg.Print_epoch_checks:
#        if unknowns==unknowns_E1:
#            print('Unknowns in both epoch match')
#        else:
#            print("Unknowns in epochs don't match")

X_vector, X_vectorHR = fc.filling_X(
                      Aproximates, unknowns, count_unknowns, count_instruments)

# filling G matrix

G_matrix = np.zeros([count_unknowns,4])

sumx=0
sumy=0
sumz=0
c=0
for i,unknown in enumerate(unknowns):
    if i < (len(unknowns) - 2*count_instruments):
        sumx = sumx + Aproximates[unknown][0]
        sumy = sumy + Aproximates[unknown][0]
        sumz = sumz + Aproximates[unknown][0]
        c=c+1
meanx = sumx/c
meany = sumy/c
meanz = sumz/c
del c,sumx,sumy,sumz

for i,unknown in enumerate(unknowns):
    if i < (len(unknowns) - 2*count_instruments):
        iii = 3*i
        G_matrix[iii,0]   = 1 
        G_matrix[iii+1,1] = 1
        G_matrix[iii+2,2] = 1
        G_matrix[iii,3] =     (Aproximates[unknown][1] - meany)/10000
        G_matrix[iii+1,3] = - (Aproximates[unknown][0] - meanx)/10000
del i, iii, unknown

    
def Filling_A_L_P_LX0(Nominal_coords,Aproximates,
                      count_all_observations,
                      count_constraints,
                      Combinations_for_constraints,
                      measured_distances_in_lines,
                      Pol_measurements,
                      ):
    A_matrix = np.zeros([count_all_observations , count_unknowns])
    """+ count_constraints"""                                                       
    A_matrixHR = {}
    
    L_vector = np.array([])
    LX0_vector = np.array([])
    
    L_vectorHR = []
    
    P_vector = np.array([])
    
    
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
            StDev = fc.StDev_sys_ppm(distance,cg.IFM_StDev)
            P_vector = np.append(P_vector,pow(cg.Sigma_0,2) / pow(StDev,2))
            """Now figuring the index of points in the unknowns list, so I know
               which column of the A matrix. The original index is multiplied 
               by 3 because unknowns are names of points (and instrument 
               orientations) 
               not a complete list of unknowns = Point X,Y and Z => *3      """
            PointFrom = sorted_measured_points_in_lines[line][meas_i]
            PointTo = sorted_measured_points_in_lines[line][meas_i+1]
            LX0_vector = np.append(LX0_vector, fc.slope_distance(
                                  Aproximates[PointTo],Aproximates[PointFrom]))
            PointFrom_i = 3*unknowns.index(PointFrom) 
            PointTo_i = 3*unknowns.index(PointTo)
            dX,dY,dZ = fc.ParD_sd(Aproximates[PointTo],Aproximates[PointFrom])
            A_matrix[L_i,PointTo_i:PointTo_i+3] = dX,dY,dZ
            # for point "From" the partial derivatives change sign
            A_matrix[L_i,PointFrom_i:PointFrom_i+3] = -dX,-dY,-dZ
            # Documenting in human readible format the A and L elements
            L_vectorHR.append((line, PointFrom, PointTo))
            A_matrixHR[(L_i,PointFrom_i)] = ['dX', line, PointFrom, PointTo]
            A_matrixHR[(L_i,PointFrom_i+1)] = ['dY', line, PointFrom, PointTo]
            A_matrixHR[(L_i,PointFrom_i+2)] = ['dZ', line, PointFrom, PointTo]
            A_matrixHR[(L_i,PointTo_i)] = ['dX', line, PointTo, PointFrom]
            A_matrixHR[(L_i,PointTo_i+1)] = ['dY', line, PointTo, PointFrom]
            A_matrixHR[(L_i,PointTo_i+2)] = ['dZ', line, PointTo, PointFrom]
    del PointFrom,PointTo, PointFrom_i,PointTo_i,dX,dY,dZ,L_i, meas_i, distance
    
    L_subv_Hz = np.array([])
    L_subv_V = np.array([])
    L_subv_sd = np.array([])
    L_subv_Hz_HR = []
    L_subv_V_HR = []
    L_subv_sd_HR = []
    
    P_subv_Hz = np.array([])
    P_subv_V = np.array([])
    P_subv_sd = np.array([])
    
    LX0_subv_Hz = np.array([])
    LX0_subv_V = np.array([])
    LX0_subv_sd = np.array([])
    
    counter = 0
    for inst,points in Pol_measurements.items():
        instrument_i = 3 * unknowns.index(inst) # Starting column of the instrument
        Ori_inst_i = X_vectorHR.index('Ori_'+inst) # Inst's Orientation index
        for point in points:
            # evaluating the partial derivatives for all polar observations
            Hz_dX, Hz_dY, Hz_dZ, Hz_dO = fc.ParD_Hz(Aproximates[point],
                                                    Aproximates[inst])
            V_dX, V_dY, V_dZ = fc.ParD_V(Aproximates[point],Aproximates[inst])
            sd_dX, sd_dY, sd_dZ = fc.ParD_sd(Aproximates[point],Aproximates[inst])
            # Returning measured values
            sd,Hz,V = Pol_measurements[inst][point][0:3] #Here are angles in gons!!!
            # Filling the L subvectors for Hz,V and sd
            L_subv_Hz = np.append(L_subv_Hz, fc.angle_normalize(
                                                        fc.gon2rad(Hz),'rad'))
            P_subv_Hz = np.append(P_subv_Hz, pow(cg.Sigma_0,2) / pow(
                                              fc.gon2rad(cg.Ang_StDev/1000),2))
            Hz_angle_from_aprox = fc.angle_normalize(
                        fc.angle_normalize(
                            fc.horizontal_angle_from_Coords(Aproximates[point],
                                                Aproximates[inst]),'rad') - \
                            fc.angle_normalize(X_vector[Ori_inst_i],'rad'),
                                            'rad')
                        
            LX0_subv_Hz = np.append(LX0_subv_Hz, Hz_angle_from_aprox)
    
            L_subv_V = np.append(L_subv_V, fc.angle_normalize(
                                                          fc.gon2rad(V),'rad'))
            P_subv_V = np.append(P_subv_V, pow(cg.Sigma_0,2) / pow(
                                              fc.gon2rad(cg.Ang_StDev/1000),2))
            LX0_subv_V = np.append(LX0_subv_V, fc.angle_normalize(
                                  fc.vertical_angle_from_Coords(
                                  Aproximates[point],Aproximates[inst]),'rad'))
            
            L_subv_sd = np.append(L_subv_sd, sd)
            P_subv_sd = np.append(P_subv_sd, pow(cg.Sigma_0,2) / pow(
                                             fc.StDev_sys_ppm(sd,cg.Dist_StDev),2))
            LX0_subv_sd = np.append(LX0_subv_sd, fc.slope_distance(
                                             Aproximates[point],Aproximates[inst]))
            
            # Filling the Human readible version of L subvectors
            L_subv_Hz_HR.append(('Hz', inst, point))
            L_subv_V_HR.append(('V', inst, point))
            L_subv_sd_HR.append(('sd', inst, point))
            # Returning the starting column of the point
            Point_i = 3*unknowns.index(point)
            # Row offsets of V and sd measurements (allows filling in same loop)
            Hz_offset = count_all_observations - 3*count_Pol_measurements
            V_offset = count_all_observations - 2*count_Pol_measurements
            sd_offset = count_all_observations - count_Pol_measurements
            # Filling A with Hz partial derivatives
            A_matrix[Hz_offset+counter,Point_i:Point_i+3] = Hz_dX, Hz_dY, Hz_dZ
            A_matrix[V_offset+counter,Point_i:Point_i+3] = V_dX, V_dY, V_dZ
            A_matrix[sd_offset+counter,Point_i:Point_i+3] = sd_dX, sd_dY, sd_dZ
            A_matrix[Hz_offset+counter,Ori_inst_i] = Hz_dO
            A_matrix[Hz_offset+counter,instrument_i:instrument_i+3] = \
                                                             -Hz_dX, -Hz_dY, -Hz_dZ
            A_matrix[V_offset+counter,instrument_i:instrument_i+3] = \
                                                                -V_dX, -V_dY, -V_dZ
            A_matrix[sd_offset+counter,instrument_i:instrument_i+3] = \
                                                             -sd_dX, -sd_dY, -sd_dZ
            # Filling Human readible A matrix
            A_matrixHR[(Hz_offset+counter,Point_i)] = ['Hz/dX', inst, point]
            A_matrixHR[(Hz_offset+counter,Point_i+1)] = ['Hz/dY', inst, point]
            A_matrixHR[(Hz_offset+counter,Point_i+2)] = ['Hz/dZ', inst, point]
            A_matrixHR[(V_offset+counter,Point_i)] = ['V/dX', inst, point]
            A_matrixHR[(V_offset+counter,Point_i+1)] = ['V/dY', inst, point]
            A_matrixHR[(V_offset+counter,Point_i+2)] = ['V/dZ', inst, point]
            A_matrixHR[(sd_offset+counter,Point_i)] = ['sd/dX', inst, point]
            A_matrixHR[(sd_offset+counter,Point_i+1)] = ['sd/dY', inst, point]
            A_matrixHR[(sd_offset+counter,Point_i+2)] = ['sd/dZ', inst, point]
            counter += 1
    L_vector = np.append(L_vector,[L_subv_Hz,L_subv_V,L_subv_sd])
    LX0_vector = np.append(LX0_vector,[LX0_subv_Hz,LX0_subv_V,LX0_subv_sd])
    L_vectorHR = L_vectorHR + L_subv_Hz_HR + L_subv_V_HR + L_subv_sd_HR
    P_vector = np.append(P_vector,[P_subv_Hz,P_subv_V,P_subv_sd])
    
    # Filling A and L, and A and L Human Readable with constraint rows
    
#    L_lenght = len(L_vector)
#    for index, const in enumerate(Combinations_for_constraints):
#        A_row_index = L_lenght + index
#        PointFrom = Combinations_for_constraints[index][0]
#        PointTo = Combinations_for_constraints[index][1]
#        #Probably needs to take the coordinates from somewhere else!
#        dX,dY,dZ = fc.ParD_sd(Nominal_coords[PointTo],Nominal_coords[PointFrom])
#        sd = fc.slope_distance(Nominal_coords[PointFrom],Nominal_coords[PointTo])
#        sd_aprox = fc.slope_distance(Aproximates[PointFrom],Aproximates[PointTo])
#        PointFrom_i = 3*unknowns.index(PointFrom) 
#        PointTo_i = 3*unknowns.index(PointTo)
#        A_matrix[A_row_index,PointTo_i:PointTo_i+3] = dX,dY,dZ
#        # for point "From" the partial derivatives change sign
#        A_matrix[A_row_index,PointFrom_i:PointFrom_i+3] = -dX,-dY,-dZ
#        # Documenting in human readible format what are the A and L elements
#        L_vector = np.append(L_vector,sd)
#        LX0_vector = np.append(LX0_vector,sd_aprox)
#        P_vector = np.append(P_vector,pow(cg.Sigma_0,2)/pow(cg.Constraint_StDev,2))
#        L_vectorHR.append(('constraint', PointFrom, PointTo))
#        A_matrixHR[(A_row_index,PointFrom_i)] = ['dX', 'distance constraint', 
#                   PointFrom, PointTo]
#        A_matrixHR[(A_row_index,PointFrom_i+1)] = ['dY', 'distance constraint', 
#                   PointFrom, PointTo]
#        A_matrixHR[(A_row_index,PointFrom_i+2)] = ['dZ', 'distance constraint', 
#                   PointFrom, PointTo]
#        A_matrixHR[(A_row_index,PointTo_i)] = ['dX', 'distance constraint', 
#                   PointTo, PointFrom]
#        A_matrixHR[(A_row_index,PointTo_i+1)] = ['dY', 'distance constraint', 
#                   PointTo, PointFrom]
#        A_matrixHR[(A_row_index,PointTo_i+2)] = ['dZ', 'distance constraint', 
#                   PointTo, PointFrom]
    P_matrix = np.diagflat(P_vector)
        
    LSM_can_be_done = (len(P_vector) ==  len(L_vector)) and (
                    len(L_vector) == A_matrix.shape[0]) and (
                                         A_matrix.shape[0] == len(LX0_vector))
    if not LSM_can_be_done:
        print("Error during filling P, L, A, and LX0, sizes are not same")
  
    if len(X_vector) != A_matrix.shape[1]:
        print("Error during filling A, X, sizes don't match.")
        LSM_can_be_done = False
    return LSM_can_be_done,A_matrix,L_vector,P_matrix,\
           LX0_vector,A_matrixHR,L_vectorHR, Hz_offset, sd_offset


LSM_can_be_done,A_matrix,L_vector,P_matrix,LX0_vector,\
A_matrixHR,L_vectorHR, Hz_offset, sd_offset \
      = Filling_A_L_P_LX0(MC.Nominal_coords,Aproximates,
                          count_all_observations,
                          count_constraints,
                          Combinations_for_constraints,
                          MC.measured_distances_in_lines,
                          MC.Pol_measurements)

Aproximates_original = Aproximates.copy()
LX0_vector_original = LX0_vector.copy()
X_vector_original = X_vector.copy()

threshold = 0.01 #fraction of basic unit
metric = threshold + 1
counter = 0
<<<<<<< Updated upstream
while (metric > threshold) and (counter < 0):
    print('\n Iteration', counter)
#    l = LX0_vector - L_vector
    l = L_vector - LX0_vector
#    for i in range(len(LX0_vector)):
#        print(i,l[i],LX0_vector[i],L_vector[i])
=======

while (metric > threshold) and (counter < 1):
    print('\n Iteration', counter)
    l = LX0_vector - L_vector
    for i in range(Hz_offset,sd_offset-1):
        l[i] = fc.angle_normalize(l[i],'rad')
#    if l[1] > 6.283:
#        l[1] = l[1]-2*m.pi
#    elif l[1] < -6.238:
#        l[1] = l[1]+2*m.pi
>>>>>>> Stashed changes
    N = A_matrix.transpose().dot(A_matrix)#.dot(P_matrix)
    print(np.linalg.det(A_matrix.transpose().dot(A_matrix)))
    O = np.zeros([4,4])
    N_extended = np.block([[N,G_matrix],[np.transpose(G_matrix),O]])
    print('Rank N_extended:', np.linalg.matrix_rank(N_extended))
    print("Determinant of G-extended N: ", np.linalg.det(N_extended))
    print("Condition number N: ", np.linalg.cond(N_extended))
    n = A_matrix.transpose().dot(l)#.dot(P_matrix)
    N_inv = inv(N_extended)[:-4,:-4]
#    print("N_inv max: ",np.amax(N_inv), " min: ",np.amin(N_inv) )#~10^7 
    dx = -N_inv.dot(n)
    print('dx',max(abs(dx)), np.argmax(abs(dx)))
    X_vector += dx
    vI = A_matrix.dot(dx) - l
    print('vI', vI[np.argmax(abs(vI))], np.argmax(abs(vI)))
    Aproximates = fc.filling_Aproximates(unknowns, X_vector, instruments)
    LSM_can_be_done,A_matrix,L_vector,P_matrix,LX0_vector,\
    A_matrixHR,L_vectorHR, Hz_offset, sd_offset \
      = Filling_A_L_P_LX0(MC.Nominal_coords,Aproximates,
                          count_all_observations,
                          count_constraints,
                          Combinations_for_constraints,
                          MC.measured_distances_in_lines,
                          MC.Pol_measurements)
    vII = - LX0_vector + L_vector
#    if vII[1] > 6.283:
#        vII[1] = vII[1]-2*m.pi
#    elif vII[1] < -6.238:
#        vII[1] = vII[1]+2*m.pi
    v = vI-vII
    print('vII', vII[np.argmax(abs(vII))], np.argmax(abs(vII)))
    print('v', v[np.argmax(abs(v))], np.argmax(abs(v)))
    metric = max(abs(dx))
    counter += 1

print('End of A_matrix code')