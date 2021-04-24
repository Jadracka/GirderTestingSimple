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

#from operator import itemgetter
#from collections import namedtuple

# =============================================================================
# Loading measurement files and Coordinates, if two epochs are set, files load
# for them as well.
# =============================================================================
Nominal_coords = fc.Coords_read_in(cg.Coords_file_name)
LoS_measurements = fc.Measurements_read_in(cg.LoS_Measurements_file_name)
Pol_measurements = fc.Measurements_read_in(cg.Pol_Measurements_file_name)
if len(cg.Which_epochs) > 1:
    LoS_measurements_1 = fc.Measurements_read_in(
                                               cg.LoS_Measurements_file_name_1)
    Pol_measurements_1 = fc.Measurements_read_in(
                                               cg.Pol_Measurements_file_name_1)

# =============================================================================
# Initial checks for LoS measurements.
# If you want to print, change: Print_typos to True
# =============================================================================

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
        print("Line %s was measured, but not expected." % (line))
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
                print("Not all points were measured in %s line"% (line))
        del line_points_sorted
    for point in LoS_measurements[line]:
        if (cg.Print_typos) and (point not in Nominal_coords):
            print("Measured point with name %s in %s is not in the Nominal" \
                  " Coordinate file." % (point, line))
            LoS_measured_points_all_good = False
        if point not in cg.Lines_of_sight[line]:
            LoS_measured_points_all_good = False
            if cg.Print_typos:
                print("Point %s does not nominally belong to %s line"
                      % (point, line))
        if point not in all_measured_points:
            all_measured_points.append(point)
    del line, point
for point in Nominal_coords.keys():
    if (cg.Print_typos) and (point not in all_measured_points):
        print("Point %s was not measured in any line." % (point))
        nominal_points_all_measured = False
    del point

if (cg.Print_typos) and (measured_lines_all_good):
    print("All measured lines were expected, no typos found.")
if (cg.Print_typos) and (LoS_measured_points_all_good):
    print("All measured points are correct, in correct lines, no typos found.")
if (cg.Print_typos) and (nominal_points_all_measured):
    print("All nominal points in IFM lines were measured at least once.")
if (cg.Print_typos) and not (all_points_in_lines_measured):
    print("Not all points in lines were measured. Continuing in analysis.")


del all_measured_points, nominal_points_all_measured, \
    all_points_in_lines_measured

Pol_measurements_cart = {}
for instrument in Pol_measurements:
    Pol_measurements_cart[instrument] = {}
    for point in Pol_measurements[instrument]:
        Pol_measurements_cart[instrument][point] = fc.polar2cart3D(
                                          Pol_measurements[instrument][point])
    del instrument, point

for instrument in Pol_measurements_cart:
    points = tuple(Pol_measurements_cart[instrument].keys())
    for i in range (1,len(points)):
        Measured = fc.slope_distance(Pol_measurements_cart[instrument][
                     points[i]],Pol_measurements_cart[instrument][points[i-1]])
        if points[i-1] not in Nominal_coords.keys() and cg.Print_typos:
            print('Point %s measured by %s is not in Nominals.' %(points[i-1],
                  instrument))
        if (points[i] in Nominal_coords.keys()) and ((points[i-1]) in
                                                        Nominal_coords.keys()):
            Nominal = fc.slope_distance(Nominal_coords[points[i]],
                                        Nominal_coords[points[i-1]])
            delta = Nominal - Measured
    del instrument, i

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
            delta = (LoS_measurements[line][
                       sorted_measured_points_in_lines[line][i]][0]\
                    - LoS_measurements[line][
                       sorted_measured_points_in_lines[line][i-1]][0],)
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
    del line, i, delta, deltas, d
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
    print("Analysis cannot be performed as there are typos and errors in "
          "input data. Please correct before running the script again. "
          "To help troubleshoot, change Print_typos in config.py to True.")

G17 = fc.StDev_XYZ_from_Polar(Nominal_coords['Girder_17'],cg.Dist_StDev,
                            cg.Ang_StDev,cg.Ang_StDev)
G18 = fc.StDev_XYZ_from_Polar(Nominal_coords['Girder_18'],cg.Dist_StDev,
                            cg.Ang_StDev,cg.Ang_StDev)
StDev1817 = fc.StDev_distance(Nominal_coords['Girder_17'],
                    Nominal_coords['Girder_18'],G17,G18)
print(StDev1817)
print('End of the script')
