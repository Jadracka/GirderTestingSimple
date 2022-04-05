# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:02:05 2021

@author: jbarker
"""

import scipy.stats as spst
import numpy as np

import sys
#import string
import math as m
import config as cg
import functions as fc
import Helmert3Dtransform as ht
from angle import Angle as a

from datetime import datetime as dt
date_time = dt.fromtimestamp(dt.timestamp(dt.now()))

# Are there two epochs to calculate?
if len(cg.Which_epochs)>1:
    Two_epochs = True
    Epoch_num = cg.Which_epochs[0]
    Epoch_num1 = cg.Which_epochs[1]
else:
    Two_epochs = False
    Epoch_num = cg.Which_epochs[0]


# =============================================================================
# Loading measurement files and Coordinates, if two epochs are set, files load
# for them as well.
# Measurements indicated for exclusion at config will be eliminated.
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
                                             Sd_StDev = cg.Dist_StDev_E1,
                                             Hz_StDev = cg.Hz_StDev_E1,
                                             V_StDev = cg.V_StDev_E1)
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
        print("Epoch_%s: Point %s was not measured in any line." 
			  % (str(Epoch_num),point))
        nominal_points_all_measured = False
    del point

if (cg.Print_real2nominal_checks) and (measured_lines_all_good):
    print("Epoch_%s: All measured lines were expected, no typos found." 
		  % (str(Epoch_num)))
if (cg.Print_real2nominal_checks) and (LoS_measured_points_all_good):
    print("Epoch_%s: All measured points are correct, in correct lines, no ty"
          "pos found." % (str(Epoch_num)))
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
                np.longdouble(len(LoS_measurements[line]))
    average_V = sum(v[2] for v in LoS_measurements[line].values()) /\
                np.longdouble(len(LoS_measurements[line]))
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
                               cg.Max_diff_from_line, abs(Diff),
							   abs(Hz_diff), abs(V_diff)))

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
                print("Epoch_%s: Measured point with name %s in %s is not in "
                      "the Nominal Coordinate file." 
					  % (str(Epoch_num1),point, line))
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
                    np.longdouble(len(LoS_measurements_E1[line]))
        average_V = sum(v[2] for v in LoS_measurements_E1[line].values()) /\
                    np.longdouble(len(LoS_measurements_E1[line]))
        counter = 0
        for point in LoS_measurements_E1[line]:
            Hz_diff = fc.gon2rad(average_Hz - LoS_measurements_E1[line][
            point][1])\
                        * LoS_measurements_E1[line][point][0]
            V_diff = fc.gon2rad(average_V - \
								LoS_measurements_E1[line][point][2])\
								* LoS_measurements_E1[line][point][0]
            Diff = m.sqrt(m.pow(Hz_diff,2)+m.pow(V_diff,2))
            if cg.Max_diff_from_line < Diff:
                counter = counter + 1
                print("Line: %s, in Epoch_%s, point %s exceeds Maximum differ"
                      "ence from line of %1.3f. The total difference is %1.3f "
                      "mm, with horizontal component %1.3f mm and vertical "
                      "component %1.3f mm" %(line, str(Epoch_num1), point,
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
        std = (fc.StDev_sys_ppm(distance,cg.IFM_StDev))/1000
        stdev_distance = stdev_distance + (std,)
    StDevs_IFM_measurements[line] = stdev_distance
    del line, stdev_distance, std, distance

if Two_epochs:
    StDevs_IFM_measurements_E1 = {}
    for line in measured_distances_in_lines_E1:
        stdev_distance = ()
        for distance in measured_distances_in_lines_E1[line]:
            std = (fc.StDev_sys_ppm(distance,cg.IFM_StDev_E1))/1000
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
print("Comparisons and initial imports done.")
# =============================================================================
# Calculating Helmert transformations for measured cartesian coordinates
# =============================================================================
Transformed_Pol_measurements, Trans_par = fc.Helmert_calc_for_PolMeas(
                                          Pol_measurements_cart,Nominal_coords)

for meas in cg.LSM_Excluded_measurements[str(Epoch_num)]:
    Pol_measurements[meas[1]][meas[2]].pop(meas[0])
try:
    del meas
except NameError:
    pass

count_IFM_measurements = sum([len(v) for k, v in\
                                         measured_distances_in_lines.items()])

count_Sd = fc.Count_meas_types(Pol_measurements, 'Sd')
count_Hz = fc.Count_meas_types(Pol_measurements, 'Hz')
count_V = fc.Count_meas_types(Pol_measurements, 'V')

count_Pol_measurements = (sum([len(v) for k, v in Pol_measurements.items()]))
excluded_count = len(cg.LSM_Excluded_measurements[str(Epoch_num)])

if count_Sd + count_Hz + count_V == 3*count_Pol_measurements - excluded_count:
    count_all_observations = count_Sd + count_Hz + count_V + count_IFM_measurements
else:
    sys.exit("Counts of measurements don't agree.")
    

unknowns,count_unknowns, instruments, count_instruments = fc.find_unknowns(
                            Transformed_Pol_measurements, cg.Instruments_6DoF)
Aproximates = fc.merge_measured_coordinates(Transformed_Pol_measurements)

if cg.LSM_incl_Cons:
    for magnet in cg.Names_of_magnets:
        try:
            magnet_fids = {k: v for k, v in cg.FIDS.items() if k.startswith(
																														magnet)}
            xAp = ht.Helmert_transform(magnet_fids,Aproximates)
            magnet_fids_trans_Ap = ht.Transformation(xAp, magnet_fids)
            xNc = ht.Helmert_transform(magnet_fids,Nominal_coords)
            magnet_fids_trans_Nc = ht.Transformation(xNc, magnet_fids)
            for key in magnet_fids:
                if key not in Aproximates: #Assert only triggres on virtual poitns.
                    Aproximates[key] = magnet_fids_trans_Ap[key]
                    Nominal_coords[key] = magnet_fids_trans_Nc[key]
                    unknowns.insert(-2*count_instruments, key)
                    count_unknowns += 3
            Beam_axis_analysis = True
        except KeyError:
            print("Magnet %s is not in FIDS. Analysis of beam axis cannot be"
					      "performed." %(magnet))
            Beam_axis_analysis = False

if cg.Instruments_6DoF:
    for instrument in Trans_par:
        Angles = Trans_par[instrument][-3:]
        Rx = (a(-Angles[0],a.T_RAD, True).angle)
        Ry = (a(-Angles[1],a.T_RAD, True).angle)
        Rz = (a(-Angles[2],a.T_RAD, True).angle)
        Aproximates['Ori_'+instrument] = (Rx, Ry, Rz)
    del Rx, Ry, Rz, Angles
else:
    for instrument in Trans_par:
        Aproximates['Ori_'+instrument] = (a(-Trans_par[instrument][-1],
                                          a.T_RAD, True).angle)
            
print("Initial Helmert transform pretransport epoch, unknown counts and Aproximates filling done.")
# =============================================================================
# Least Square Method for pre-transport epoch
# =============================================================================
P_matrix, Results, Qxx, Qvv, Cov_matrix, s02, dof, w, s02_IFM, s02_Hz, s02_V,  \
s02_Sd, s02_con = fc.LSM(Epoch_num, 
										Nominal_coords, Aproximates,
									   measured_distances_in_lines,				
									   sorted_measured_points_in_lines,
									   instruments, count_instruments,
									   Pol_measurements,count_Pol_measurements,
									   count_IFM_measurements,
									   unknowns,count_unknowns, cg.IFM_StDev, 
                       cg.Instruments_6DoF, Trans_par, cg.Epsilon)

CovM_filename = str("Cov_matrix_" + str(Epoch_num) + ".txt")

if not cg.LSM_incl_Cons:
    np.savetxt(CovM_filename, Cov_matrix)

print("LSM for pre transport epoch done.")


# =============================================================================
# Calculating Helmert transformations for measured cartesian coordinates
# =============================================================================

if Two_epochs:


    Transformed_Pol_measurements_E1, Trans_par_E1 = fc.Helmert_calc_for_PolMeas(
                                   Pol_measurements_cart_E1,Nominal_coords_E1)
    
    count_IFM_measurements_E1 = sum([len(v) for k, v in\
                                    measured_distances_in_lines_E1.items()])
    
    count_Pol_measurements_E1 = (sum([len(v) for k, v in \
                                   Pol_measurements_E1.items()]))
    for meas in cg.LSM_Excluded_measurements[str(Epoch_num1)]:
        Pol_measurements_E1[meas[1]][meas[2]].pop(meas[0])
    try:
        del meas
    except NameError:
        pass
    excluded_count_E1 = len(cg.LSM_Excluded_measurements[str(Epoch_num1)])
    count_Sd_E1 = fc.Count_meas_types(Pol_measurements_E1, 'Sd')
    count_Hz_E1 = fc.Count_meas_types(Pol_measurements_E1, 'Hz')
    count_V_E1 = fc.Count_meas_types(Pol_measurements_E1, 'V')
 
    if count_Sd + count_Hz + count_V == \
                            3*count_Pol_measurements - excluded_count:
        count_all_observations = count_Sd + count_Hz + count_V + \
                                 count_IFM_measurements
    else:
        sys.exit("Counts of measurements don't agree.")
    
    unknowns_E1,count_unknowns_E1, instruments_E1, \
    count_instruments_E1 = fc.find_unknowns(Transformed_Pol_measurements_E1, 
                                            cg.Instruments_6DoF)
            
    Aproximates_E1 = fc.merge_measured_coordinates(
											 Transformed_Pol_measurements_E1)

    if cg.LSM_incl_Cons:
        for magnet in cg.Names_of_magnets:
            try:
                magnet_fids_E1 = {k: v for k, v in cg.FIDS.items() if \
																							  k.startswith(magnet)}
                xAp_E1 = ht.Helmert_transform(magnet_fids_E1,Aproximates_E1)
                magnet_fids_trans_Ap_E1 = ht.Transformation(xAp_E1, 
																											magnet_fids_E1)
                xNc_E1 = ht.Helmert_transform(magnet_fids_E1,Nominal_coords_E1)
                magnet_fids_trans_Nc_E1 = ht.Transformation(xNc_E1, 
																											magnet_fids_E1)
                for key in magnet_fids_E1:
                    if key not in Aproximates_E1:
                        Aproximates_E1[key] = magnet_fids_trans_Ap_E1[key]
                        Nominal_coords_E1[key] = magnet_fids_trans_Nc_E1[key]
                        unknowns_E1.insert(-2*count_instruments_E1, key)
                        count_unknowns_E1 += 3
                Beam_axis_analysis = True
            except KeyError:
                print("Magnet %s is not in FIDS. Analysis of beam axis cannot"
						        " be performed." %(magnet))
                Beam_axis_analysis = False


    if cg.Instruments_6DoF:
        for instrument in Trans_par:
            Aproximates_E1['Ori_'+instrument] = Trans_par_E1[instrument][-3:]
    else:
        for instrument in Trans_par:
            Aproximates_E1['Ori_'+instrument] = Trans_par_E1[instrument][-1]
    
# =============================================================================
# Least Square Method for post-transport epoch
# =============================================================================
    P_matrix_E1,Results_E1, Qxx_E1, Qvv_E1, Cov_matrix_E1, s02_E1, dof_E1, w_E1, s02_IFM_E1,\
	s02_Hz_E1, s02_V_E1, s02_Sd_E1, s02_con_E1 = fc.LSM(Epoch_num1, 
										Nominal_coords_E1, Aproximates_E1, 
										measured_distances_in_lines_E1,				
									   sorted_measured_points_in_lines_E1,
									   instruments_E1, count_instruments_E1,
									   Pol_measurements_E1,
									   count_Pol_measurements_E1,
									   count_IFM_measurements_E1, unknowns_E1, 
                       count_unknowns_E1, cg.IFM_StDev_E1, 
                       cg.Instruments_6DoF, Trans_par_E1, cg.Epsilon) 
    CovM_filename_E1 = str("Cov_matrix_" + str(Epoch_num1) + ".txt")
    if not cg.LSM_incl_Cons:
        np.savetxt(CovM_filename_E1, Cov_matrix_E1)
    del Cov_matrix_E1
del Cov_matrix

print('End of LSM \n')

# =============================================================================
# Calculating standard deviations for unknowns from CovMatrix, 
# getting rid of orientation in Result dictionary
# The StD are calculated from Cov Matrix (real cov m) without constraints!
# =============================================================================
try:
    Cov_matrix = np.array(np.loadtxt(CovM_filename, dtype=float))
    covmfile_exists = True
except:
    print(CovM_filename, "doesn't exist, run the program first without"
									    " constraints")
    covmfile_exists = False
    
if Two_epochs:	
	try:
	    Cov_matrix_E1 = np.array(np.loadtxt(CovM_filename_E1, dtype=float))
	except:
	    print(CovM_filename_E1, "doesn't exist, run the program first without"
										    " constraints")
	    covmfile_exists = False
    
if covmfile_exists and cg.LSM_incl_Cons:
	for i, unknown in enumerate(unknowns[:-(count_instruments*2 + 12)]):
	    iii = 3*i
	    point = Results[unknown]
	    st_dev_X = m.sqrt(Cov_matrix[iii,iii])
	    st_dev_Y = m.sqrt(Cov_matrix[iii+1,iii+1])
	    st_dev_Z = m.sqrt(Cov_matrix[iii+2,iii+2])
	    st_dev_XYZ = m.sqrt(pow(st_dev_X,2)+pow(st_dev_Y,2)+pow(st_dev_Z,2))
	    update = point + (st_dev_X, st_dev_Y, st_dev_Z, st_dev_XYZ)
	    Results[unknown] = update
	del update, point, st_dev_X, st_dev_Y, st_dev_Z, st_dev_XYZ
	
#	if cg.Instruments_6DoF:
#	    for i, unknown in enumerate(unknowns[-count_instruments:]):
#	        iii = 3*i
#	        point = Results[unknown]
#	        st_dev_Rx = m.sqrt(Cov_matrix[iii,iii])
#	        st_dev_Ry = m.sqrt(Cov_matrix[iii+1,iii+1])
#	        st_dev_Rz = m.sqrt(Cov_matrix[iii+2,iii+2])
#	        update = point + (st_dev_Rx, st_dev_Ry, st_dev_Rz)
#	        Results[unknown] = update
#	    del update, point, st_dev_Rx, st_dev_Ry, st_dev_Rz
	
#	if not cg.Instruments_6DoF:
#	    Ori_keys = ()
#	    for key in Results.keys():
#	        if "Ori_" in key:
#	            Ori_keys += (key,)
#	        for key in Ori_keys:
#	            Results.pop(key)
#	        del i, unknown, iii, Ori_keys, key
	
	if Two_epochs:	
	    for i, unknown in enumerate(unknowns_E1[:-(count_instruments_E1*2 + 12)]):
	        iii = 3*i
	        point = Results_E1[unknown]
	        st_dev_X = m.sqrt(Cov_matrix_E1[iii,iii])
	        st_dev_Y = m.sqrt(Cov_matrix_E1[iii+1,iii+1])
	        st_dev_Z = m.sqrt(Cov_matrix_E1[iii+2,iii+2])
	        st_dev_XYZ = m.sqrt(pow(st_dev_X,2)+pow(st_dev_Y,2)+pow(st_dev_Z,2))
	        update = point + (st_dev_X, st_dev_Y, st_dev_Z, st_dev_XYZ)
	        Results_E1[unknown] = update
	    del update, point, st_dev_X, st_dev_Y, st_dev_Z, st_dev_XYZ
	    if cg.Instruments_6DoF:
	        for i, unknown in enumerate(unknowns_E1[-count_instruments_E1:]):
	            iii = 3*i
	            point = Results[unknown]
	            st_dev_Rx = m.sqrt(Cov_matrix_E1[iii,iii])
	            st_dev_Ry = m.sqrt(Cov_matrix_E1[iii+1,iii+1])
	            st_dev_Rz = m.sqrt(Cov_matrix_E1[iii+2,iii+2])
	            update = point + (st_dev_Rx, st_dev_Ry, st_dev_Rz)
	            Results_E1[unknown] = update
	        del update, point, st_dev_Rx, st_dev_Ry, st_dev_Rz
	
	    if not cg.Instruments_6DoF:
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
	        diffs = tuple(np.array(Results[point][:3]) - np.array(
																				Transformed_Results_E1[point][:3]))
	        mag = round(fc.slope_distance(Results[point],
																			   Transformed_Results_E1[point]),10)
	        result = diffs + (mag,)
			# delta X, Y, Z and magnitude:
	        Movements[point] = result
	    del point, diffs, mag, result
	# =============================================================================
	# Deformation comparison between two epochs
	# =============================================================================
	Results_2 = {}
	Results_E1_2 = {}
	Results_def_2 = {}
	exclude_points = set(['Aus', 'Oben', 'Ein', 'Mitte', 'Instrument'])
	if Two_epochs:
	    identicals = list(set(unknowns) & set(unknowns_E1))
	    for point in identicals:
	        if not any(name in point for name in exclude_points):
	            F = spst.f.ppf(q=1-0.05, dfn=3, dfd=dof-4)
	            F_E1 = spst.f.ppf(q=1-0.05, dfn=3, dfd=dof_E1-4)
	            F_def = spst.f.ppf(q=1-0.05, dfn=3, dfd=dof_E1+dof-8)
	            iii = 3*(unknowns.index(point)) 
	            iii_E1 = 3*(unknowns_E1.index(point))
	            sub_matrix = Cov_matrix[iii:iii+3, iii:iii+3]
	            sub_matrix_E1 = Cov_matrix_E1[iii_E1:iii_E1+3, iii_E1:iii_E1+3]
	            def_matrix = sub_matrix + sub_matrix_E1
	            eigen_val, eigen_vec = np.linalg.eig(sub_matrix)
	            eigen_val_E1, eigen_vec_E1 = np.linalg.eig(sub_matrix_E1)
	            eigen_val_def, eigen_vec_def = np.linalg.eig(def_matrix)
	            # For first epoch:
	            AoA = m.atan2(eigen_vec[0,1],eigen_vec[0,0])
	            ZoA = m.acos(eigen_vec[0,2]/np.linalg.norm((eigen_vec[:,0])))
	            AoB = m.atan2(eigen_vec[1,1],eigen_vec[1,0])
	            ZoB = m.acos(eigen_vec[1,2]/np.linalg.norm((eigen_vec[:,1])))
	            AoC = m.atan2(eigen_vec[2,1],eigen_vec[2,0])
	            ZoC = m.acos(eigen_vec[2,2]/np.linalg.norm((eigen_vec[:,2])))
	            Ac = m.sqrt(3*F*eigen_val[0])
	            Bc = m.sqrt(3*F*eigen_val[1])
	            Cc = m.sqrt(3*F*eigen_val[2])
	            Results_2[point] = (AoA, ZoA, AoB, ZoB, AoC, ZoC, Ac, Bc, Cc)
				     # For second epoch:
	            AoA_E1 = m.atan2(eigen_vec_E1[0,1],eigen_vec_E1[0,0])
	            ZoA_E1 = m.acos(eigen_vec_E1[0,2]/np.linalg.norm((
																									eigen_vec_E1[:,0])))
	            AoB_E1 = m.atan2(eigen_vec_E1[1,1],eigen_vec_E1[1,0])
	            ZoB_E1 = m.acos(eigen_vec_E1[1,2]/np.linalg.norm((
																									eigen_vec_E1[:,1])))
	            AoC_E1 = m.atan2(eigen_vec_E1[2,1],eigen_vec_E1[2,0])
	            ZoC_E1 = m.acos(eigen_vec_E1[2,2]/np.linalg.norm((
																									eigen_vec_E1[:,2])))
	            Ac_E1 = m.sqrt(3*F_E1*eigen_val_E1[0])
	            Bc_E1 = m.sqrt(3*F_E1*eigen_val_E1[1])
	            Cc_E1 = m.sqrt(3*F_E1*eigen_val_E1[2])
	            Results_E1_2[point] = (AoA_E1, ZoA_E1, AoB_E1, ZoB_E1, AoC_E1, 
															  ZoC_E1, Ac_E1, Bc_E1, Cc_E1)
						# For deformations:
	            AoA_def = m.atan2(eigen_vec_def[0,1],eigen_vec_def[0,0])
	            ZoA_def = m.acos(eigen_vec_def[0,2]/np.linalg.norm((
																									eigen_vec_def[:,0])))
	            AoB_def = m.atan2(eigen_vec_def[1,1],eigen_vec_def[1,0])
	            ZoB_def = m.acos(eigen_vec_def[1,2]/np.linalg.norm((
																									eigen_vec_def[:,1])))
	            AoC_def = m.atan2(eigen_vec_def[2,1],eigen_vec_def[2,0])
	            ZoC_def = m.acos(eigen_vec_def[2,2]/np.linalg.norm((
																									eigen_vec_def[:,2])))
	            Ac_def = m.sqrt(3*F_def*eigen_val_def[0])
	            Bc_def = m.sqrt(3*F_def*eigen_val_def[1])
	            Cc_def = m.sqrt(3*F_def*eigen_val_def[2])
				
	            xA, yA, zA = fc.polar2cart3Drad(1, AoA_def, ZoA_def)
	            xB, yB, zB = fc.polar2cart3Drad(1, AoB_def, ZoB_def)
	            xC, yC, zC = fc.polar2cart3Drad(1, AoC_def, ZoC_def)
	
	            Mov_vec = np.asarray(Movements[point][:-1])
	            A_cart = np.array([xA,yA,zA])
	            dotpA = np.dot(A_cart,Mov_vec)
	            B_cart = np.array([xB,yB,zB])
	            dotpB = np.dot(B_cart,Mov_vec)
	            C_cart = np.array([xC,yC,zC])
	            dotpC = np.dot(C_cart,Mov_vec)
				
	            np.dot(A_cart, B_cart)
	            print(point, m.pow(dotpA/Ac_def,2) + m.pow(dotpB/Bc_def,2) + m.pow(dotpC/Cc_def,2))
	            is_out = m.pow(dotpA/Ac_def,2) + m.pow(dotpB/Bc_def,2) + m.pow(
																									dotpC/Cc_def,2) > 1
	
	            Movements[point] += (is_out,)#np.dot(A_cart, B_cart)/min(Ac_def,
						            #Bc_def), np.dot(A_cart, C_cart)/min(Ac_def, Cc_def),
											#np.dot(B_cart, C_cart)/min(Bc_def, Cc_def))
	            Results_def_2[point] = (AoA_def, ZoA_def, AoB_def, ZoB_def, 
															   AoC_def, ZoC_def, Ac_def, Bc_def, Cc_def)
#	        else:
#	            F = spst.f.ppf(q=1-0.05, dfn=1, dfd=dof)
#	            F_E1 = spst.f.ppf(q=1-0.05, dfn=1, dfd=dof_E1)
#	            if cg.Instruments_6DoF:
#	                print(point)
#	                iii = 3 * (unknowns.index(point)-12)
#	                Rxc = m.sqrt(F * s02 * Cov_matrix[iii,iii])
#	                Ryc = m.sqrt(F * s02 * Cov_matrix[iii+1,iii+1])
#	                Rzc = m.sqrt(F * s02 * Cov_matrix[iii+2,iii+2])
#	                Results_2[point] = (Rxc, Ryc, Rzc)
#	            else:
#	                i = len(unknowns) - unknowns.index(point)
#	                Oric = m.sqrt(F * s02 * Cov_matrix[-i, -i])
#	                Results_2[point] = (Oric)
	    del point, iii, iii_E1       




# =============================================================================
# LSM - results writing into file - pre-transport Epoch
# =============================================================================

Results_file = open(cg.Res_file_name, "w")

Header = ["Results from Epoch" + str(
        cg.Which_epochs[0]) + " [RHCS]\n", "created:" + str(date_time)
        + "\nUsing source files:\n" + '-' + str(
                cg.LoS_Measurements_file_name) + '\n', '-' + str(
                cg.Pol_Measurements_file_name) + '\n', '-' + str(
                cg.Coords_file_name) + '\n']

Results_file.writelines(Header)

for point in Aproximates:
    if 'Ori' in point:
        fill =  '\n' + point + '\t' + str(
                                    Results[point][0]) + '\t' + str(
                                   -Results[point][1]) + '\t' + str(
                                    Results[point][2])
        Results_file.write(str(fill))
    elif 'Aus' or 'Oben' or 'Mitte' or 'Aus' in point:
        fill =  '\n' + point + '\t' + str(
                                    Results[point][0]) + '\t' + str(
                                   -Results[point][1]) + '\t' + str(
                                    Results[point][2])
        Results_file.write(str(fill))
    else:
        fill =  '\n' + point + '\t' + str(
                                    Results[point][0]) + '\t' + str(
                                   -Results[point][1]) + '\t' + str(
                                    Results[point][2]) + '\t' + str(
                                    Results[point][3]) + '\t' + str(
                                    Results[point][4]) + '\t' + str(
                                    Results[point][5]) + '\t' + str(
                                    Results[point][6])

        Results_file.write(str(fill))

# Closing file
Results_file.close()
del point, Header, fill

print('End of MainCode')

# =============================================================================
# LSM - results writing into file - post-transport Epoch
# =============================================================================
if Two_epochs:    
    Results_file_E1 = open(cg.Res_file_name_1, "w")
    
    Header = ["Results from Epoch" + str(
            cg.Which_epochs[0]) + " [RHCS]\n", "created:" + str(date_time)
            + "\nUsing source files:\n" + ' -' + str(
                    cg.LoS_Measurements_file_name_1) + '\n', ' -' + str(
                    cg.Pol_Measurements_file_name_1) + '\n', ' -' + str(
                    cg.Coords_file_name_1) + '\n'
                    ]
    
    Results_file_E1.writelines(Header)
    
    for point in Aproximates_E1:
        if 'Ori' in point:
            fill =  '\n' + point + '\t' + str(
                                        Results_E1[point][0]) + '\t' + str(
                                       -Results_E1[point][1]) + '\t' + str(
                                        Results_E1[point][2])
            Results_file_E1.write(str(fill))
        elif 'Aus' or 'Oben' or 'Mitte' or 'Aus' in point:
            fill =  '\n' + point + '\t' + str(
                                        Results_E1[point][0]) + '\t' + str(
                                       -Results_E1[point][1]) + '\t' + str(
                                        Results_E1[point][2])
            Results_file_E1.write(str(fill))
        else:
            fill =  '\n' + point + '\t' + str(
                                        Results_E1[point][0]) + '\t' + str(
                                       -Results_E1[point][1]) + '\t' + str(
                                        Results_E1[point][2]) + '\t' + str(
                                        Results_E1[point][3]) + '\t' + str(
                                        Results_E1[point][4]) + '\t' + str(
                                        Results_E1[point][5]) + '\t' + str(
                                        Results_E1[point][6])

            Results_file_E1.write(str(fill))
    
    # Closing file
    Results_file_E1.close()
    del Header, fill, point