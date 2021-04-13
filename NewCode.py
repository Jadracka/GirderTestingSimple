# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:02:05 2021

@author: jbarker
"""

#import scipy as sp
import re

import numpy as np

#import sys
#import string
#import math as m

import config as cg
import functions as fc

#from operator import itemgetter
#from collections import namedtuple


Coords_file = open(cg.Coords_file_name,'r')
Nominal_coords = {}

for line in Coords_file.readlines():
    """Reads in and parses coordinates file. Ignores point notes at the end.
       It can handle multiple occurences of delimeters,
       but not a combination of them."""
    words = re.split(';+|,+|\t+| +',line.strip())
    Nominal_coords[words[0]] = (float(words[1]), 
                                float(words[2]), 
                                float(words[3]))
    del line
    
Coords_file.close()

LoS_Meas_file = open(cg.Measurements_file_name,'r')
LoS_measurements = {}

for row in LoS_Meas_file.readlines():
    """
    - Reads in and parses Meas_file. 
    - The format is string Line name 'space' 
      string Point name 'space' float Distance [mm] 'space' 
      float Hz angle [gon] 'space' float Z angle [gon].
    - Ignores point notes at the end.
    - It can handle multiple occurences of delimeters but not a 
      combination of them.
        - LoS_measurements - a Dictionary of lines which contains a 
                                    Dictionary of points, where Point
                                    name is a key and measured values
                                    are triplet tuple"""
    words = re.split(';+|,+|\t+| +',row.strip())
    if words[0] not in LoS_measurements.keys():
        LoS_measurements[words[0]] = {}
        LoS_measurements[words[0]][words[1]] =                         \
          (float(words[2]), float(words[3]), float(words[4]))
    else: 
        LoS_measurements[words[0]][words[1]] =                         \
          (float(words[2]), float(words[3]), float(words[4]))
    del words, row
LoS_Meas_file.close()

# =============================================================================
# Initial checks for LoS measurements.
# If you want to print, change: Print_typos to True
# =============================================================================

# Checking point names for typos and misspells and creating the list of
# measured points in lines sorted based on the config file definition:
sorted_measured_points_in_lines = {}
all_measured_points = []
measured_lines_all_good = True
measured_points_all_good = True
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
            measured_points_all_good = False
        if point not in cg.Lines_of_sight[line]:
            measured_points_all_good = False
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
if (cg.Print_typos) and (measured_points_all_good):
    print("All measured points are correct, in correct lines, no typos found.")   
if (cg.Print_typos) and (nominal_points_all_measured):
    print("All nominal points were measured at least once.")
if (cg.Print_typos) and not (all_points_in_lines_measured):
    print("Not all points in lines were measured. Continueing in analysis.")


del all_measured_points, nominal_points_all_measured,\
    all_points_in_lines_measured

#sorted_measured_points_in_lines

if measured_lines_all_good and measured_points_all_good:
    del measured_lines_all_good, measured_points_all_good
    # Calculating distance deltas
    for line in LoS_measurements:
        for i in range (1,len(sorted_measured_points_in_lines[line])):
#            delta = fc.slope_distance
            print(i) #, delta
            pass
        #fc.slope_distance
    if cg.Using_nominal_compare0:
        pass
else:
    print("Analysis cannot be performed as there are typos and errors in " 
          "input data. Please correct before running the script again. "
          "To help troubleshoot, change Print_typos in config.py to True.")

print('End of the script')