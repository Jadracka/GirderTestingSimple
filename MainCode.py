# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:02:05 2021

@author: jbarker
"""

#import scipy as sp
import numpy as np

#import sys
#import string
import math as m
import config as cg
import functions as fc
from numpy.linalg import inv
from angle import Angle as a
import Helmert3Dtransform as ht

from datetime import datetime as dt
date_time = dt.fromtimestamp(dt.timestamp(dt.now()))
#import Helmert3Dtransform as helm

#from operator import itemgetter
#from collections import namedtuple

# Are there two epochs to calculate?
if len(cg.Which_epochs)>1:
    Two_epochs = True
else:
    Two_epochs = False

# =============================================================================
# Loading measurement files and Coordinates, if two epochs are set, files load
# for them as well.
# =============================================================================
Nominal_coords = fc.Coords_read_in(cg.Coords_file_name)
LoS_measurements = fc.Measurements_read_in(cg.LoS_Measurements_file_name)
Pol_measurements = fc.Polar_2F_meas_read_in(cg.Pol_Measurements_file_name,
                                            Sd_StDev = cg.Dist_StDev,
                                            Hz_StDev = cg.Hz_StDev,
                                            V_StDev = cg.V_StDev)
if Two_epochs:
    LoS_measurements_E1 = fc.Measurements_read_in(
                                               cg.LoS_Measurements_file_name_1)
    Pol_measurements_E1 = fc.Polar_2F_meas_read_in(
                                             cg.Pol_Measurements_file_name_1,
                                             Sd_StDev = cg.Dist_StDev,
                                             Hz_StDev = cg.Hz_StDev,
                                             V_StDev = cg.V_StDev)
    Nominal_coords_E1 = fc.Coords_read_in(cg.Coords_file_name_1)


# =============================================================================
# Initial checks for LoS measurements.
# If you want to print, change: Print_typos to True in config.py
# =============================================================================
"""
EPOCH 0 - Pre transport
"""
# Checking point names for typos and misspells and creating the list of
# measured points in lines sorted based on the config file definition:
sorted_measured_points_in_lines = {}
all_measured_points = []
measured_lines_all_good = True
LoS_measured_points_all_good = True
#nominal_lines_all_measured = True #not checking at the moment
nominal_points_all_measured = True
all_points_in_lines_measured = True

for line in LoS_measurements:
    Epoch_num = cg.Which_epochs[0]
    if (cg.Print_typos) and (line not in cg.Lines_of_sight):
# printing which lines are in measurements input but are not in the
# default naming either due to typo or just simply missing in the nominal
# LoS decription
        print("Epoch_%s: Line %s was measured, but not expected." % (str(Epoch_num), line))
        measured_lines_all_good = False
    else:
        line_points_sorted = []
        for point in cg.Lines_of_sight[line]:
            if point in LoS_measurements[line].keys():
                line_points_sorted.append(point)
        sorted_measured_points_in_lines[line] = tuple(line_points_sorted)
        if line_points_sorted != list(cg.Lines_of_sight[line]):
            all_points_in_lines_measured = False
            if cg.Print_typos:
                print("Epoch_%s: Not all points were measured in %s line"
                      % (str(Epoch_num), line))
        del line_points_sorted
    for point in LoS_measurements[line]:
        if (cg.Print_typos) and (point not in Nominal_coords):
            print("Epoch_%s: Measured point with name %s in %s is not in the "
                  "Nominal Coordinate file." % (str(Epoch_num), point, line))
            LoS_measured_points_all_good = False
        if point not in cg.Lines_of_sight[line]:
            LoS_measured_points_all_good = False
            if cg.Print_typos:
                print("Epoch_%s: Point %s does not nominally belong to %s line"
                      % (str(Epoch_num), point, line))
        if point not in all_measured_points:
            all_measured_points.append(point)
    del line, point
for point in Nominal_coords.keys():
    if (cg.Print_typos) and (point not in all_measured_points):
        print("Epoch_%s: Point %s was not measured in any line." % (str(Epoch_num),point))
        nominal_points_all_measured = False
    del point

if (cg.Print_real2nominal_checks) and (measured_lines_all_good):
    print("Epoch_%s: All measured lines were expected, no typos found." % (str(Epoch_num)))
if (cg.Print_real2nominal_checks) and (LoS_measured_points_all_good):
    print("Epoch_%s: All measured points are correct, in correct lines, no typos"
          " found." % (str(Epoch_num)))
if (cg.Print_real2nominal_checks) and (nominal_points_all_measured):
    print("Epoch_%s: All nominal points in IFM lines were measured at least "
          "once." % (str(Epoch_num)))
if (cg.Print_real2nominal_checks) and not (all_points_in_lines_measured):
    print("Epoch_%s: Not all points in lines were measured. Continuing in "
          "analysis." % (str(Epoch_num)))


del all_measured_points, nominal_points_all_measured, \
    all_points_in_lines_measured

Pol_measurements_cart = {}
for instrument in Pol_measurements:
    Pol_measurements_cart[instrument] = {}
    for point in Pol_measurements[instrument]:
        Pol_measurements_cart[instrument][point] = fc.polar2cart3Dgon(
                                          Pol_measurements[instrument][point])
    del instrument, point

for instrument in Pol_measurements_cart:
    points = tuple(Pol_measurements_cart[instrument].keys())
    for i in range (1,len(points)):
        Measured = fc.slope_distance(Pol_measurements_cart[instrument][
                     points[i]],Pol_measurements_cart[instrument][points[i-1]])
        if points[i-1] not in Nominal_coords.keys() and cg.Print_typos:
            print('Epoch_%s: Point %s measured by %s is not in Nominals.'
                  %(str(Epoch_num), points[i-1], instrument))
        if (points[i] in Nominal_coords.keys()) and ((points[i-1]) in
                                                        Nominal_coords.keys()):
            Nominal = fc.slope_distance(Nominal_coords[points[i]],
                                        Nominal_coords[points[i-1]])
            delta = Nominal - Measured
    del instrument, i, points, Measured, Nominal

measured_distances_in_lines = {}
if cg.Using_nominal_compare:
    nominal_distances_in_line = {}
    differences_in_distances = {}
    StDev_distances_in_lines = {}
if measured_lines_all_good and LoS_measured_points_all_good:
    del measured_lines_all_good, LoS_measured_points_all_good
    # Calculating distance deltas
    for line in LoS_measurements:
        deltas = ()
        if cg.Using_nominal_compare:
            deltas_nominal = ()
        for i in range (1,len(sorted_measured_points_in_lines[line])):
            delta = (abs(LoS_measurements[line][
                       sorted_measured_points_in_lines[line][i]][0]\
                    - LoS_measurements[line][
                       sorted_measured_points_in_lines[line][i-1]][0]),)
            deltas = deltas + delta
            if cg.Using_nominal_compare:
                d = fc.slope_distance(
                            Nominal_coords[
                                sorted_measured_points_in_lines[line][i-1]],
                            Nominal_coords[
                                sorted_measured_points_in_lines[line][i]])
                deltas_nominal = deltas_nominal + (d,)
        measured_distances_in_lines[line] = deltas
        if cg.Using_nominal_compare:
            nominal_distances_in_line[line] = deltas_nominal
            del deltas_nominal
    del line, i, delta, deltas,
    if cg.Using_nominal_compare:
        for line in LoS_measurements:
            differences_in_distances[line] = np.asarray(
                                              nominal_distances_in_line[line])\
                                            - np.asarray(
                                            measured_distances_in_lines[line])
            StDev_distances_in_lines[line] = np.std(
                                               differences_in_distances[line])
            if StDev_distances_in_lines[line] == 0:
                StDev_distances_in_lines[line] = None
            differences_in_distances[line] = tuple(
                                                differences_in_distances[line])
        del line
else:
    print("Analysis for Epoch_%s cannot be performed as there are typos and "
          "errors in input data. Please correct before running the script "
          "again. To help troubleshoot, change Print_typos in config.py to "
          "True." % (str(Epoch_num)))


# Checking the distance from line for LT-IFM measurements
for line in LoS_measurements:
    average_Hz = sum(v[1] for v in LoS_measurements[line].values()) /\
                float(len(LoS_measurements[line]))
    average_V = sum(v[2] for v in LoS_measurements[line].values()) /\
                float(len(LoS_measurements[line]))
    counter = 0
    for point in LoS_measurements[line]:
        Hz_diff = fc.gon2rad(average_Hz - LoS_measurements[line][point][1])\
                    * LoS_measurements[line][point][0]
        V_diff = fc.gon2rad(average_V - LoS_measurements[line][point][2])\
                    * LoS_measurements[line][point][0]
        Diff = m.sqrt(m.pow(Hz_diff,2)+m.pow(V_diff,2))
        if cg.Max_diff_from_line < Diff:
            counter = counter + 1
            print("Line: %s, in Epoch_%s, point %s exceeds Maximum difference "
                  "from line of %1.3f. The total difference is %1.3f mm, "
                  " with horizontal component %1.3f mm and vertical component "
                  "%1.3f mm" %(line, str(Epoch_num), point, 
                               cg.Max_diff_from_line, abs(Diff),abs(Hz_diff), abs(V_diff)))

del line, average_Hz, average_V, counter, point, Hz_diff, V_diff, Diff


"""
EPOCH 1 - Post transport
"""

if Two_epochs:
    # Checking point names for typos and misspells and creating the list of
    # measured points in lines sorted based on the config file definition:
    sorted_measured_points_in_lines_E1 = {}
    all_measured_points_E1 = []
    measured_lines_all_good_E1 = True
    LoS_measured_points_all_good_E1 = True
    #nominal_lines_all_measured = True #not checking at the moment
    nominal_points_all_measured_E1 = True
    all_points_in_lines_measured_E1 = True

    for line in LoS_measurements_E1:
        Epoch_num1 = cg.Which_epochs[1]
        if (cg.Print_typos) and (line not in cg.Lines_of_sight):
    # printing which lines are in measurements input but are not in the
    # default naming either due to typo or just simply missing in the nominal
    # LoS decription
            print("Epoch_%s: Line %s was measured, but not expected." \
                              % (str(Epoch_num1),line))
            measured_lines_all_good_E1 = False
        else:
            line_points_sorted = []
            for point in cg.Lines_of_sight[line]:
                if point in LoS_measurements_E1[line].keys():
                    line_points_sorted.append(point)
            sorted_measured_points_in_lines_E1[line]= tuple(line_points_sorted)
            if line_points_sorted != list(cg.Lines_of_sight[line]):
                all_points_in_lines_measured_E1 = False
                if cg.Print_typos:
                    print("Epoch_%s: Not all points were measured in %s line"
                          % (str(Epoch_num1),line))
            del line_points_sorted
        for point in LoS_measurements_E1[line]:
            if (cg.Print_typos) and (point not in Nominal_coords):
                print("Epoch_%s: Measured point with name %s in %s is not in the"
                      " Nominal Coordinate file." % (str(Epoch_num1),point, line))
                LoS_measured_points_all_good_E1 = False
            if point not in cg.Lines_of_sight[line]:
                LoS_measured_points_all_good_E1 = False
                if cg.Print_typos:
                    print("Epoch_%s: Point %s does not nominally belong to %s "
                          "line" % (str(Epoch_num1),point, line))
            if point not in all_measured_points_E1:
                all_measured_points_E1.append(point)
        del line, point
    for point in Nominal_coords.keys():
        if (cg.Print_typos) and (point not in all_measured_points_E1):
            print("Epoch_%s: Point %s was not measured in any line." \
                              % (str(Epoch_num1),point))
            nominal_points_all_measured_E1 = False
        del point

    if (cg.Print_real2nominal_checks) and (measured_lines_all_good_E1):
        print("Epoch_%s: All measured lines were expected, no typos found." \
                          %(str(Epoch_num1)))
    if (cg.Print_real2nominal_checks) and (LoS_measured_points_all_good_E1):
        print("Epoch_%s: All measured points are correct, in correct lines, no"
              " typos found." % (str(Epoch_num1)))
    if (cg.Print_real2nominal_checks) and (nominal_points_all_measured_E1):
        print("Epoch_%s: All nominal points in IFM lines were measured at "
              "least once." % (str(Epoch_num1)))
    if (cg.Print_real2nominal_checks) and not (all_points_in_lines_measured_E1):
        print("Epoch_%s: Not all points in lines were measured. Continuing in "
              "analysis." % (str(Epoch_num1)))

    del all_measured_points_E1, nominal_points_all_measured_E1,\
        all_points_in_lines_measured_E1

    Pol_measurements_cart_E1 = {}
    for instrument in Pol_measurements_E1:
        Pol_measurements_cart_E1[instrument] = {}
        for point in Pol_measurements_E1[instrument]:
            Pol_measurements_cart_E1[instrument][point] = fc.polar2cart3Dgon(
                                        Pol_measurements_E1[instrument][point])
        del instrument, point

    for instrument in Pol_measurements_cart_E1:
        points = tuple(Pol_measurements_cart_E1[instrument].keys())
        for i in range(1, len(points)):
            Measured = fc.slope_distance(Pol_measurements_cart_E1[instrument][
                         points[i]],
                         Pol_measurements_cart_E1[instrument][points[i-1]])
            if points[i-1] not in Nominal_coords.keys() and\
                           cg.Print_real2nominal_checks:
                print('Epoch_%s: Point %s measured by %s is not in Nominals.'
                      % (str(Epoch_num1), points[i-1], instrument))
            if (points[i] in Nominal_coords.keys()) and ((points[i-1]) in
                                                        Nominal_coords.keys()):
                Nominal = fc.slope_distance(Nominal_coords[points[i]],
                                            Nominal_coords[points[i-1]])
                delta = Nominal - Measured
        del instrument, i, points, Measured, Nominal

    measured_distances_in_lines_E1 = {}
    if cg.Using_nominal_compare:
        nominal_distances_in_line_E1 = {}
        differences_in_distances_E1 = {}
        StDev_distances_in_lines_E1 = {}
    if measured_lines_all_good_E1 and LoS_measured_points_all_good_E1:
        del measured_lines_all_good_E1, LoS_measured_points_all_good_E1
        # Calculating distance deltas
        for line in LoS_measurements_E1:
            deltas = ()
            if cg.Using_nominal_compare:
                deltas_nominal = ()
            for i in range (1,len(sorted_measured_points_in_lines_E1[line])):
                delta = (abs(LoS_measurements_E1[line][
                           sorted_measured_points_in_lines_E1[line][i]][0]
                        - LoS_measurements_E1[line][
                           sorted_measured_points_in_lines_E1[line][i-1]][0]),)
                deltas = deltas + delta
                if cg.Using_nominal_compare:
                    d = fc.slope_distance(
                             Nominal_coords[
                                sorted_measured_points_in_lines_E1[line][i-1]],
                             Nominal_coords[
                                sorted_measured_points_in_lines_E1[line][i]])
                    deltas_nominal = deltas_nominal + (d,)
            measured_distances_in_lines_E1[line] = deltas
            if cg.Using_nominal_compare:
                nominal_distances_in_line_E1[line] = deltas_nominal
                del deltas_nominal
        del line, i, delta, deltas
        if cg.Using_nominal_compare:
            for line in LoS_measurements_E1:
                differences_in_distances_E1[line] = \
                      np.asarray(nominal_distances_in_line_E1[line]) \
                    - np.asarray(measured_distances_in_lines_E1[line])

                StDev_distances_in_lines_E1[line] = np.std(
                                             differences_in_distances_E1[line])
                if StDev_distances_in_lines_E1[line] == 0:
                    StDev_distances_in_lines_E1[line] = None
                differences_in_distances_E1[line] = tuple(
                                             differences_in_distances_E1[line])
            del line
    else:
        print("Analysis for Epoch_%s cannot be performed as there are typos and"
              " errors in input data. Please correct before running the script"
              " again. To help troubleshoot, change Print_typos in config.py "
              "to True." % (str(Epoch_num1)))

# Checking the distance from line for LT-IFM measurements
if Two_epochs:
    for line in LoS_measurements_E1:
        average_Hz = sum(v[1] for v in LoS_measurements_E1[line].values()) /\
                    float(len(LoS_measurements_E1[line]))
        average_V = sum(v[2] for v in LoS_measurements_E1[line].values()) /\
                    float(len(LoS_measurements_E1[line]))
        counter = 0
        for point in LoS_measurements_E1[line]:
            Hz_diff = fc.gon2rad(average_Hz - LoS_measurements_E1[line][
            point][1])\
                        * LoS_measurements_E1[line][point][0]
            V_diff = fc.gon2rad(average_V - LoS_measurements_E1[line][point][2])\
                        * LoS_measurements_E1[line][point][0]
            Diff = m.sqrt(m.pow(Hz_diff,2)+m.pow(V_diff,2))
            if cg.Max_diff_from_line < Diff:
                counter = counter + 1
                print("Line: %s, in Epoch_%s, point %s exceeds Maximum differen"
                      "ce from line of %1.3f. The total difference is %1.3f mm"
                      ", with horizontal component %1.3f mm and vertical compo"
                      "nent %1.3f mm" %(line, str(Epoch_num1), point,
                      cg.Max_diff_from_line, abs(Diff), abs(Hz_diff),
                      abs(V_diff)))
    #    print(line, counter)
    del line, average_Hz, average_V, counter, point, Hz_diff, V_diff, Diff

# =============================================================================
# Standard Deviations calculations
# =============================================================================

# Calculating StDevs for Laser Tracker IFM measurements
StDevs_IFM_measurements = {}
for line in measured_distances_in_lines:
    stdev_distance = ()
    for distance in measured_distances_in_lines[line]:
        std = fc.StDev_sys_ppm(distance,cg.IFM_StDev)
        stdev_distance = stdev_distance + (std,)
    StDevs_IFM_measurements[line] = stdev_distance
    del line, stdev_distance, std, distance

if Two_epochs:
    StDevs_IFM_measurements_E1 = {}
    for line in measured_distances_in_lines_E1:
        stdev_distance = ()
        for distance in measured_distances_in_lines_E1[line]:
            std = fc.StDev_sys_ppm(distance,cg.IFM_StDev)
            stdev_distance = stdev_distance + (std,)
        StDevs_IFM_measurements_E1[line] = stdev_distance
    del line, stdev_distance, std, distance



# =============================================================================
# EPOCH comparisons, only happens if there are 2 Epochs
# =============================================================================
all_lines_measured_same = True
if Two_epochs:
    for line in sorted_measured_points_in_lines:
        # Checking if same Lines of Sight were measured in both Epochs
        if sorted_measured_points_in_lines[line] != \
                                      sorted_measured_points_in_lines_E1[line]:
            all_lines_measured_same = False
    if not all_lines_measured_same:
        print("Lines weren't measured in the same manner, please correct!")
    if all_lines_measured_same and cg.Print_epoch_checks:
        print('All lines measured in Epoch_%s were also measured in Epoch_%s.'
              % (str(Epoch_num),str(Epoch_num1)))
    del line, all_lines_measured_same

# =============================================================================
# Calculating Helmert transformations for measured cartesian coordinates
# =============================================================================
Transformed_Pol_measurements = fc.Helmert_calc_for_PolMeas(
                                          Pol_measurements_cart,Nominal_coords)

count_IFM_measurements = sum([len(v) for k, v in\
                                         measured_distances_in_lines.items()])

count_Pol_measurements = (sum([len(v) for k, v in Pol_measurements.items()]))

count_all_observations = 3*count_Pol_measurements + count_IFM_measurements

del count_IFM_measurements

unknowns,count_unknowns, instruments, count_instruments = fc.find_unknowns(
                                                  Transformed_Pol_measurements)
Aproximates = fc.merge_measured_coordinates(Transformed_Pol_measurements)
#Aproximates_copy = Aproximates.copy()

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
        sumx += Aproximates[unknown][0]
        sumy += Aproximates[unknown][1]
        sumz += Aproximates[unknown][2]
        c+=1
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
del i, iii, unknown, meany, meanx, meanz


# =============================================================================
# Least Square Method for pre-transport epoch
# =============================================================================

LSM_can_be_done,A_matrix,L_vector,P_matrix,LX0_vector,A_matrixHR,L_vectorHR \
      = fc.Filling_A_L_P_LX0(Nominal_coords,Aproximates,
                             Combinations_for_constraints,
                             measured_distances_in_lines,
                             sorted_measured_points_in_lines,
                             instruments, count_instruments,
                             Pol_measurements,unknowns,count_unknowns,
                             X_vector, X_vectorHR
                             )


metric = cg.LSM_Threshold + 1
counter = 0
print('\n Processing Epoch_%s:' %(str(Epoch_num)))
while (metric > cg.LSM_Threshold) and (counter <= cg.LSM_Max_iterations):
    print('\n Iteration', counter)
    l = np.ndarray(L_vector.size)
    for i, lelement in enumerate(L_vector):
        if L_vectorHR[i][0] == "Hz":
            l[i] = (a(LX0_vector[i] - L_vector[i],a.T_RAD,True).angle)
        else:
            l[i] = LX0_vector[i] - L_vector[i]
    del lelement, i
    N = A_matrix.transpose() @ P_matrix @ A_matrix
    print(np.linalg.det(A_matrix.transpose() @ A_matrix))
    O = np.zeros([4,4])
    N_extended = np.block([[N,G_matrix],[np.transpose(G_matrix),O]])
    print('Rank N_extended:', np.linalg.matrix_rank(N_extended))
    try:
        print("Log-Determinant of G-extended N: ", np.linalg.slogdet(N_extended))
    except:
        print("Log-Determinant of G-extended N gives overflow ")
    print("Condition number N: ", np.linalg.cond(N_extended))
    n = A_matrix.transpose() @ P_matrix @ l
    N_inv = inv(N_extended)[:-4,:-4]
    dx = N_inv @ n
    print('dx',max(abs(dx)), np.argmax(abs(dx)))
    X_vector += dx
    v = A_matrix @ dx - l
    print('v max ', v[np.argmax(abs(v))], np.argmax(abs(v)))
    print('[vv]  ', sum(v*v))
    print('[v]  ', sum(v))

    Aproximates = fc.filling_Aproximates(unknowns, X_vector, instruments)
    LSM_can_be_done,A_matrix,L_vector,P_matrix,LX0_vector,A_matrixHR, \
    L_vectorHR = fc.Filling_A_L_P_LX0(Nominal_coords,Aproximates,
                                      Combinations_for_constraints,
                                      measured_distances_in_lines,
                                      sorted_measured_points_in_lines,
                                      instruments, count_instruments,
                                      Pol_measurements,unknowns,count_unknowns,
                                      X_vector, X_vectorHR
                                      )

    metric = max(abs(dx))
    counter += 1
del metric, counter

# =============================================================================
# LSM - results writing into file - pre-transport Epoch
# =============================================================================

Results = open(cg.Res_file_name, "w")

Header = ["Results from Epoch" + str(
        cg.Which_epochs[0]) + " [RHCS]\n", "created:" + str(date_time)
        + "\nUsing source files:\n" + '-' + str(
                cg.LoS_Measurements_file_name) + '\n', '-' + str(
                cg.Pol_Measurements_file_name) + '\n', '-' + str(
                cg.Coords_file_name) + '\n']

Results.writelines(Header)

for point in Aproximates:
    if 'Ori' not in point:
        fill =  '\n' + point + '\t' + str(
                                    Aproximates[point][0]) + '\t' + str(
                                    -Aproximates[point][1]) + '\t' + str(
                                    Aproximates[point][2])

        Results.write(str(fill))

# Closing file
Results.close()
del point, Header, fill

if Two_epochs:
# =============================================================================
# Calculating Helmert transformations for measured cartesian coordinates
# =============================================================================

    Transformed_Pol_measurements_E1 = fc.Helmert_calc_for_PolMeas(
                                   Pol_measurements_cart_E1,Nominal_coords_E1)
    
    count_IFM_measurements_E1 = sum([len(v) for k, v in\
                                    measured_distances_in_lines_E1.items()])
    
    count_Pol_measurements_E1 = (sum([len(v) for k, v in \
                                   Pol_measurements_E1.items()]))
    
    count_all_observations_E1 = 3*count_Pol_measurements_E1 + \
                                                 count_IFM_measurements_E1
    
    del count_IFM_measurements_E1
    
    unknowns_E1,count_unknowns_E1, instruments_E1, \
    count_instruments_E1 = fc.find_unknowns(Transformed_Pol_measurements_E1)
            
    Aproximates_E1 = fc.merge_measured_coordinates(
											 Transformed_Pol_measurements_E1)
    #Aproximates_E1_copy = Aproximates_E1.copy()
    
    Combinations_for_constraints_E1 = []
    for magnet in cg.Names_of_magnets:
        points_on_magnet = [key for key in Aproximates_E1.keys() if magnet in key]
        magnet_combinations = [(PointA, PointB) for index, PointA in enumerate(
            points_on_magnet) for PointB in points_on_magnet[index + 1:]]
        Combinations_for_constraints_E1.extend(magnet_combinations)
    count_constraints = len(Combinations_for_constraints_E1)
    
    
    X_vector_E1, X_vectorHR_E1 = fc.filling_X(Aproximates_E1, unknowns_E1,
                                              count_unknowns_E1, 
                                              count_instruments_E1)
    
    # filling G matrix
    
    G_matrix_E1 = np.zeros([count_unknowns_E1,4])
    
    sumx=0
    sumy=0
    sumz=0
    c=0
    for i,unknown in enumerate(unknowns_E1):
        if i < (len(unknowns_E1) - 2*count_instruments_E1):
            sumx += Aproximates_E1[unknown][0]
            sumy += Aproximates_E1[unknown][1]
            sumz += Aproximates_E1[unknown][2]
            c+=1
    meanx = sumx/c
    meany = sumy/c
    meanz = sumz/c
    del c,sumx,sumy,sumz
    
    for i,unknown in enumerate(unknowns_E1):
        if i < (len(unknowns_E1) - 2*count_instruments_E1):
            iii = 3*i
            G_matrix_E1[iii,0]   = 1
            G_matrix_E1[iii+1,1] = 1
            G_matrix_E1[iii+2,2] = 1
            G_matrix_E1[iii,3] =     (Aproximates_E1[unknown][1] - meany)/10000
            G_matrix_E1[iii+1,3] = - (Aproximates_E1[unknown][0] - meanx)/10000
    del i, iii, unknown, meany, meanx, meanz
    
    
# =============================================================================
# Least Square Method for pre-transport epoch
# =============================================================================
    
    LSM_can_be_done_E1,A_matrix_E1,L_vector_E1,P_matrix_E1,LX0_vector_E1,\
    A_matrixHR_E1,L_vectorHR_E1 = fc.Filling_A_L_P_LX0(
                                 Nominal_coords_E1,Aproximates_E1,
                                 Combinations_for_constraints_E1,
                                 measured_distances_in_lines_E1,
                                 sorted_measured_points_in_lines_E1,
                                 instruments_E1, count_instruments_E1,
                                 Pol_measurements_E1,unknowns_E1,
                                 count_unknowns_E1, X_vector_E1, X_vectorHR_E1
                                 )
	
    
    metric = cg.LSM_Threshold + 1
    counter = 0
    print('\n Processing Epoch_%s:' %(str(Epoch_num1)))
    while (metric > cg.LSM_Threshold) and (counter <= cg.LSM_Max_iterations):
        print('\n Iteration', counter)
        l_E1 = np.ndarray(L_vector_E1.size)
        for i, lelement in enumerate(L_vector_E1):
            if L_vectorHR_E1[i][0] == "Hz":
                l_E1[i] = (a(LX0_vector_E1[i] - \
						 L_vector_E1[i],a.T_RAD,True).angle)
            else:
                l_E1[i] = LX0_vector_E1[i] - L_vector_E1[i]
        del lelement, i
        N_E1 = A_matrix_E1.transpose() @ P_matrix_E1 @ A_matrix_E1
        print(np.linalg.det(A_matrix_E1.transpose() @ A_matrix_E1))
        O_E1 = np.zeros([4,4])
        N_extended_E1 = np.block([[N_E1,G_matrix_E1],
                                  [np.transpose(G_matrix_E1),O_E1]])
        print('Rank N_extended:', np.linalg.matrix_rank(N_extended_E1))
        try:
            print("Log-Determinant of G-extended N: ",
				   np.linalg.slogdet(N_extended_E1))
        except:
            print("Log-Determinant of G-extended N gives overflow ")
        print("Condition number N: ", np.linalg.cond(N_extended_E1))
        n_E1 = A_matrix_E1.transpose() @ P_matrix_E1 @ l_E1
        N_inv_E1 = inv(N_extended_E1)[:-4,:-4]
        dx_E1 = N_inv_E1 @ n_E1
        print('dx',max(abs(dx_E1)), np.argmax(abs(dx_E1)))
        X_vector_E1 += dx_E1
        v_E1 = A_matrix_E1 @ dx_E1 - l_E1
        print('v max ', v_E1[np.argmax(abs(v_E1))], np.argmax(abs(v_E1)))
        print('[vv]  ', sum(v_E1*v_E1))
        print('[v]  ', sum(v_E1))
		

		
    
        Aproximates_E1 = fc.filling_Aproximates(unknowns_E1, X_vector_E1, 
                                                instruments_E1)
		
        LSM_can_be_done_E1,A_matrix_E1,L_vector_E1,P_matrix_E1,LX0_vector_E1,\
        A_matrixHR_E1, L_vectorHR_E1 = fc.Filling_A_L_P_LX0(
                                          Nominal_coords_E1,Aproximates_E1,
                                          Combinations_for_constraints_E1,
                                          measured_distances_in_lines_E1,
                                          sorted_measured_points_in_lines_E1,
                                          instruments_E1, count_instruments_E1,
                                          Pol_measurements_E1,unknowns_E1,
                                          count_unknowns_E1, X_vector_E1, 
                                          X_vectorHR_E1
                                          )
    
        metric = max(abs(dx_E1))
        counter += 1
    del metric, counter
# =============================================================================
# LSM - results writing into file - pre-transport Epoch
# =============================================================================
    
    Results_E1 = open(cg.Res_file_name_1, "w")
    
    Header = ["Results from Epoch" + str(
            cg.Which_epochs[0]) + " [RHCS]\n", "created:" + str(date_time)
            + "\nUsing source files:\n" + '-' + str(
                    cg.LoS_Measurements_file_name_1) + '\n', '-' + str(
                    cg.Pol_Measurements_file_name_1) + '\n', '-' + str(
                    cg.Coords_file_name_1) + '\n']
    
    Results_E1.writelines(Header)
    
    for point in Aproximates_E1:
        if 'Ori' not in point:
            fill =  '\n' + point + '\t' + str(
                                        Aproximates_E1[point][0]) + '\t' + str(
                                       -Aproximates_E1[point][1]) + '\t' + str(
                                        Aproximates_E1[point][2])

            Results_E1.write(str(fill))
    
    # Closing file
    Results_E1.close()

print('End of LSM \n')

if Two_epochs:
	Results = Aproximates.copy()
	Results_E1 = Aproximates_E1.copy()
	Cov_matrix = N_inv.copy()
	Cov_matrix_E1 = N_inv_E1.copy()
# =============================================================================
# Calculating standard deviations for unknowns from CovMatrix, 
# getting rid of orientation in Result dictionary
# =============================================================================
	for i, unknown in enumerate(unknowns[:-count_instruments]):
		iii = 3*i
		point = Results[unknown]
		st_dev_X = m.sqrt(Cov_matrix[iii,iii])
		st_dev_Y = m.sqrt(Cov_matrix[iii+1,iii+1])
		st_dev_Z = m.sqrt(Cov_matrix[iii+2,iii+2])
		update = point + (st_dev_X, st_dev_Y, st_dev_Z)
		Results[unknown] = update
	del update, point, st_dev_X, st_dev_Y, st_dev_Z
	Ori_keys = ()
	for key in Results.keys():
		if "Ori_" in key:
			Ori_keys += (key,)
	for key in Ori_keys:
		Results.pop(key)
	del i, unknown, iii, Ori_keys, key
	
	for i, unknown in enumerate(unknowns_E1[:-count_instruments_E1]):
		iii = 3*i
		point = Results_E1[unknown]
		st_dev_X = m.sqrt(Cov_matrix_E1[iii,iii])
		st_dev_Y = m.sqrt(Cov_matrix_E1[iii+1,iii+1])
		st_dev_Z = m.sqrt(Cov_matrix_E1[iii+2,iii+2])
		update = point + (st_dev_X, st_dev_Y, st_dev_Z)
		Results_E1[unknown] = update
	del update, point, st_dev_X, st_dev_Y, st_dev_Z
	Ori_keys = ()
	for key in Results_E1.keys():
		if "Ori_" in key:
			Ori_keys += (key,)
	for key in Ori_keys:
		Results_E1.pop(key)
	del i, unknown, iii, Ori_keys, key
# =============================================================================
# Building dictionary of points wanted to do the transform on
# =============================================================================
	Identical_points = {}
	Identical_points_E1 = {}
	for point in cg.Common_points:
		if point in Results.keys():
			Identical_points[point] = Results[point]
		if point in Results_E1.keys():
			Identical_points_E1[point] = Results_E1[point]
			
	x = ht.Helmert_transform(Identical_points_E1, Identical_points)
	Transformed_Results_E1 = ht.Transformation(x, Results_E1)
	
	Res_set = set(Results)
	Res_E1_set = set(Transformed_Results_E1)
	Movements = {}
	for point in Res_set.intersection(Res_E1_set):
		diffs = tuple(np.array(Results[point][:3]) - \
					np.array(Transformed_Results_E1[point][:3]))
		mag = fc.slope_distance(Results[point], Transformed_Results_E1[point])
		result = diffs + (mag,)
		# delta X, Y, Z and magnitude:
		Movements[point] = result
	del point, diffs, mag, result

		

print('End of MainCode')
