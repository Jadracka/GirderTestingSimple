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
                                               cg.Pol_Measurements_file_name_1)

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
    if (cg.Print_typos) and (line not in cg.Lines_of_sight):
# printing which lines are in measurements input but are not in the
# default naming either due to typo or just simply missing in the nominal
# LoS decription
        print("Epoch0: Line %s was measured, but not expected." % (line))
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
                print("Epoch0: Not all points were measured in %s line"
                      % (line))
        del line_points_sorted
    for point in LoS_measurements[line]:
        if (cg.Print_typos) and (point not in Nominal_coords):
            print("Epoch0: Measured point with name %s in %s is not in the "
                  "Nominal Coordinate file." % (point, line))
            LoS_measured_points_all_good = False
        if point not in cg.Lines_of_sight[line]:
            LoS_measured_points_all_good = False
            if cg.Print_typos:
                print("Epoch0: Point %s does not nominally belong to %s line"
                      % (point, line))
        if point not in all_measured_points:
            all_measured_points.append(point)
    del line, point
for point in Nominal_coords.keys():
    if (cg.Print_typos) and (point not in all_measured_points):
        print("Epoch0: Point %s was not measured in any line." % (point))
        nominal_points_all_measured = False
    del point

if (cg.Print_real2nominal_checks) and (measured_lines_all_good):
    print("Epoch0: All measured lines were expected, no typos found.")
if (cg.Print_real2nominal_checks) and (LoS_measured_points_all_good):
    print("Epoch0: All measured points are correct, in correct lines, no typos"
          " found.")
if (cg.Print_real2nominal_checks) and (nominal_points_all_measured):
    print("Epoch0: All nominal points in IFM lines were measured at least "
          "once.")
if (cg.Print_real2nominal_checks) and not (all_points_in_lines_measured):
    print("Epoch0: Not all points in lines were measured. Continuing in "
          "analysis.")


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
            print('Epoch0: Point %s measured by %s is not in Nominals.'
                  %(points[i-1], instrument))
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
    print("Analysis for Epoch 0 cannot be performed as there are typos and "
          "errors in input data. Please correct before running the script "
          "again. To help troubleshoot, change Print_typos in config.py to "
          "True.")

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
            print("Line: %s, in Epoch 0, point %s exceeds Maximum difference "
                  "from line of %1.3f. The total difference is %1.3f mm, "
                  " with horizontal component %1.3f mm and vertical component "
                  "%1.3f mm" %(line, point, cg.Max_diff_from_line, abs(Diff),
                               abs(Hz_diff), abs(V_diff)))
#    print(line, counter)
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
        if (cg.Print_typos) and (line not in cg.Lines_of_sight):
    # printing which lines are in measurements input but are not in the
    # default naming either due to typo or just simply missing in the nominal
    # LoS decription
            print("Epoch1: Line %s was measured, but not expected." % (line))
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
                    print("Epoch 1: Not all points were measured in %s line"
                          % (line))
            del line_points_sorted
        for point in LoS_measurements_E1[line]:
            if (cg.Print_typos) and (point not in Nominal_coords):
                print("Epoch1: Measured point with name %s in %s is not in the"
                      " Nominal Coordinate file." % (point, line))
                LoS_measured_points_all_good_E1 = False
            if point not in cg.Lines_of_sight[line]:
                LoS_measured_points_all_good_E1 = False
                if cg.Print_typos:
                    print("Epoch1: Point %s does not nominally belong to %s "
                          "line" % (point, line))
            if point not in all_measured_points_E1:
                all_measured_points_E1.append(point)
        del line, point
    for point in Nominal_coords.keys():
        if (cg.Print_typos) and (point not in all_measured_points_E1):
            print("Epoch1: Point %s was not measured in any line." % (point))
            nominal_points_all_measured_E1 = False
        del point

    if (cg.Print_real2nominal_checks) and (measured_lines_all_good_E1):
        print("Epoch1: All measured lines were expected, no typos found.")
    if (cg.Print_real2nominal_checks) and (LoS_measured_points_all_good_E1):
        print("Epoch1: All measured points are correct, in correct lines, no "
              "typos found.")
    if (cg.Print_real2nominal_checks) and (nominal_points_all_measured_E1):
        print("Epoch1: All nominal points in IFM lines were measured at least "
              "once.")
    if (cg.Print_real2nominal_checks) and not (all_points_in_lines_measured_E1):
        print("Epoch1: Not all points in lines were measured. Continuing in "
              "analysis.")

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
                print('Epoch1: Point %s measured by %s is not in Nominals.'
                      % (points[i-1], instrument))
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
                delta = (LoS_measurements_E1[line][
                           sorted_measured_points_in_lines_E1[line][i]][0]
                        - LoS_measurements_E1[line][
                           sorted_measured_points_in_lines_E1[line][i-1]][0],)
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
        print("Analysis for Epoch 1 cannot be performed as there are typos and"
              " errors in input data. Please correct before running the script"
              " again. To help troubleshoot, change Print_typos in config.py "
              "to True.")

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
                print("Line: %s, in Epoch 1, point %s exceeds Maximum differen"
						        "ce from line of %1.3f. The total difference is %1.3f mm"
                      ", with horizontal component %1.3f mm and vertical compo"
                      "nent %1.3f mm" %(line, point, cg.Max_diff_from_line, 
										abs(Diff), abs(Hz_diff), abs(V_diff)))
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

# Calculating XYZ StDevs for polar measurements, adding the results to the
# Pol_measurements_cart(_E1) as extension of the existing tuple format:
# X, Y, Z, StDev_X, StDev_Y, StDev_Z

""" vvv DOESN'T HAVE TO BE HERE - I DO THIS IN THE DATA INPUT FUNCTION! vvv """

#for instrument in Pol_measurements_cart:
#    for point in Pol_measurements_cart[instrument]:
#        # Leica's strange way of describing angular precision to normal
#        StDev_HZ_Z = fc.gon2rad(cg.Ang_StDev/1000)
#        # Combining ADM and IFM precision for polar measurements
#        StDev_S = cg.ADM_StDev + fc.StDev_sys_ppm(Pol_measurements[instrument][
#                                                        point][0],cg.IFM_StDev)
#        # Polar measurements precisions appended to the measured dictionary
#        StDev_meas = (StDev_S,StDev_HZ_Z,StDev_HZ_Z)
#        Pol_measurements[instrument][point] = Pol_measurements[instrument][
#                                                            point] + StDev_meas
#        # X, Y, Z StDevs calculated and appended to the _cart measured data
#        StDevXYZ = fc.StDev_XYZ_from_Polar(Pol_measurements[instrument][
#                                          point],StDev_S,StDev_HZ_Z,StDev_HZ_Z)
#        Pol_measurements_cart[instrument][point] = Pol_measurements_cart[
#                instrument][point] + StDevXYZ
#    del point, instrument, StDevXYZ, StDev_HZ_Z, StDev_S, StDev_meas
#
#if Two_epochs:
#    for instrument in Pol_measurements_cart_E1:
#        for point in Pol_measurements_cart_E1[instrument]:
#            # Leica's strange way of describing angular precision to normal
#            StDev_HZ_Z = fc.StDev_angle(
#                    Pol_measurements_E1[instrument][point][0],cg.Ang_StDev)
#            # Combining ADM and IFM precision for polar measurements
#            StDev_S = cg.ADM_StDev + fc.StDev_sys_ppm(Pol_measurements_E1[
#                                            instrument][point][0],cg.IFM_StDev)
#            # Polar measurements precisions appended to the measured dictionary
#            StDev_meas = (StDev_S,StDev_HZ_Z,StDev_HZ_Z)
#            Pol_measurements_E1[instrument][point] = Pol_measurements_E1[
#                                                instrument][point] + StDev_meas
#            # X, Y, Z StDevs calculated and appended to the _cart measured data
#            StDevXYZ = fc.StDev_XYZ_from_Polar(Pol_measurements_E1[instrument][
#                                          point],StDev_S,StDev_HZ_Z,StDev_HZ_Z)
#            Pol_measurements_cart_E1[instrument][point] = \
#                         Pol_measurements_cart_E1[instrument][point] + StDevXYZ
#        del point, instrument, StDevXYZ, StDev_HZ_Z, StDev_S, StDev_meas

""" ^^^ AFTER MAKING SURE IT DOES THE SAME, CAN BE DELETED ^^^ """

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
        print('All lines measured in Epoch 0 were also measured in Epoch 1.')
    del line, all_lines_measured_same

# =============================================================================
# Calculating Helmert transformations for measured cartesian coordinates
# =============================================================================

Transformed_Pol_measurements = fc.Helmert_calc_for_PolMeas(
                                          Pol_measurements_cart,Nominal_coords)
if Two_epochs:
    Transformed_Pol_measurements_E1 = fc.Helmert_calc_for_PolMeas(
                                       Pol_measurements_cart_E1,Nominal_coords)



print('End of MainCode')
