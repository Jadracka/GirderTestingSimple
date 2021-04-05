# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 10:50:35 2021

@author: jbarker
"""

#import scipy as sp
import numpy as np
#import matplotlib as mp
#import sys
#import string
#import math as m
import re
import config as cg
import functions as fc
from operator import itemgetter
#from collections import namedtuple

Coords_file = open(cg.Coords_file_name,'r')
Nominal_coords = {}

for line in Coords_file.readlines():
    """Reads in and parses coordinates file. Ignores point notes at the end.
       It can handle multiple occurences of delimeters,
       but not a combination of them."""
    words = re.split(';+|,+|\t+| +',line.strip())
    Nominal_coords[words[0]] = (float(words[1]), float(words[2]), float(words[3]))
    del line

Coords_file.close()

if cg.Using_CAD_compare:
#Calculates nominal distances in ideal lines as seen in config file
    nominal_distances_in_lines = {}
    for line, Tupple_of_points in cg.Lines_of_sight.items():
        First_point_Coords = Nominal_coords[Tupple_of_points[0]]
        distances_list = []
        for point in Tupple_of_points[1:]:
            distance = fc.slope_distance(First_point_Coords,Nominal_coords[point])
            distances_list.append(distance)
            nominal_distances_in_lines[line] = distances_list


Meas_file = open(cg.Measurements_file_name,'r')
Meas = {}
Measured_lines = []
Points_on_measured_lines = {}
measured_distances = {}
full_measurement_matrix = {}

for row in Meas_file.readlines():
    """
    - Reads in and parses Meas_file. 
    - The format is string Group name (aka Line of sight) 'space' 
      string Point name 'space' float Distance [mm] 'space' 
      float Hz angle [gon] 'space' float V angle [gon].
    - Ignores point notes at the end.
    - It can handle multiple occurences of delimeters but not a combination of 
      them.
    - It creates
        - Meas - Dictionary of list of Tuples where the name of the Line of 
                 sight is the key and tuples are measured values to the points
        - full measurement matrix - a Dictionary of lines which contains a 
                                    Dictionary of points, where Point name is 
                                    a key and measured values are triplet tuple
        - measured_distances - Dictionary of measured disances"""
    words = re.split(';+|,+|\t+| +',row.strip())
    if str(words[0]) not in Measured_lines:
        Measured_lines.append(words[0])
    if words[0] not in Meas:
        Meas[words[0]] = [(float(words[2]), float(words[3]), float(words[4]))]
        measured_distances[words[0]] = (words[2],)
        full_measurement_matrix[words[0]] = {}
        full_measurement_matrix[words[0]][words[1]] =                         \
          (float(words[2]), float(words[3]), float(words[4]))
    else: 
        Meas[words[0]]= Meas[words[0]]+[(float(words[2]), float(words[3]),    \
          float(words[4]))]
        measured_distances[words[0]] = measured_distances[words[0]] +         \
          (words[2],)
        full_measurement_matrix[words[0]][words[1]] =                         \
          (float(words[2]), float(words[3]), float(words[4]))
    del words, row

# =============================================================================
# Finfing the minimal distance measured in each line and storing it, together
# in dictionary of lists. Dictionary keys = line names, values = list of pairs:
# Point_name, measured distance (first of the triplet tuple of measurements)
# =============================================================================

minimal_distances_in_line = {}

for line in full_measurement_matrix:
    """itemgetter(which position I want to consider)[returns what is needed 
       - without [] whole tuple]"""
    minimal_distances_in_line[line] = []
    distance = min(Meas[line],key=itemgetter(0))[0]
    """some pretty cool mumbo jumbo! Using the .values() and .keys()
       for a nested dictionary [line], searching for key of the nested 
       dictionary (point name) using .index(), where the index value is the 
       minimum distance found from the value=tuple of measured triplet of the 
       key=point"""
    minimal_distance_point_name = \
    (list(full_measurement_matrix[line].keys())    \
    [list(full_measurement_matrix[line].values()).index(min(Meas[line],\
     key=itemgetter(0)))])
    minimal_distances_in_line[line] = [minimal_distance_point_name, distance]
    del minimal_distance_point_name, distance

for dictionary in full_measurement_matrix:
    # Points on measured lines dictionary of key = line name and list = points 
    Points_on_measured_lines[dictionary] =                                    \
    list(full_measurement_matrix[dictionary].keys())
    del dictionary

# =============================================================================
# Calculating the nominal and real differential distances for each point
# of the line except the one closest to the tracker (=beginning of the line).
# Storing those in dictionary of lists. Key is the line's name.
# =============================================================================
nominal_distances_measured = {}
for line in measured_distances:
    distances_list = []
    for point in Points_on_measured_lines[line]:
        if point != minimal_distances_in_line[line][0]:
            delta_distance = fc.slope_distance(Nominal_coords                 \
                [minimal_distances_in_line[line][0]],Nominal_coords[point])
            distances_list.append(delta_distance)
    nominal_distances_measured[line] = distances_list
    del distances_list, delta_distance, line, point

real_distances_measured = {}
for line in measured_distances:
    distances_list = []
    for distance in measured_distances[line]:
        if float(distance) != float(minimal_distances_in_line[line][1]):
            delta_distance = float(distance) -                                \
                float(minimal_distances_in_line[line][1])
            distances_list.append(delta_distance)
    real_distances_measured[line] = distances_list
    del distances_list, delta_distance, line, distance

differences_real_nominal_distances = {}
for line in nominal_distances_measured:
    diff = np.subtract((nominal_distances_measured[line]),   \
                      (real_distances_measured[line]))
    differences_real_nominal_distances[line] = list(diff)
    if cg.Line_differences_checking:
        if diff[0] > 1:
            print('discrapancy in:',line,'\n')
        else:
            print('It looks ok')
    del line, diff

