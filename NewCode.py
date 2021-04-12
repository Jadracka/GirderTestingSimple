# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:02:05 2021

@author: jbarker
"""

#import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
#import sys
#import string
#import math as m
import re
import config as cg
import functions as fc
from mpl_toolkits.mplot3d import Axes3D
#from operator import itemgetter
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

if cg.Print_FIDs:
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.axis('equal')
    xdata = []
    ydata = []
    zdata = []
    PointID = []
    for key in Nominal_coords:
        PointID.append(key)
        xdata.append(Nominal_coords[key][0])
        ydata.append(Nominal_coords[key][1])
        zdata.append(Nominal_coords[key][2])
#    ax.scatter3D(xdata, ydata, zdata);
    for i in range(len(PointID)): #plot each point + it's index as text above
        ax.scatter(xdata[i],ydata[i],zdata[i],color='b') 
        ax.text(xdata[i],ydata[i],zdata[i],  '%s' % (str(PointID[i])), size=9, zorder=1, color='k')
    plt.show()

LoS_Meas_file = open(cg.Measurements_file_name,'r')
LoS_measurements = {}

for row in LoS_Meas_file.readlines():
    """
    - Reads in and parses Meas_file. 
    - The format is string Line name 'space' 
      string Point name 'space' float Distance [mm] 'space' 
      float Hz angle [gon] 'space' float Z angle [gon].
    - Ignores point notes at the end.
    - It can handle multiple occurences of delimeters but not a combination of 
      them.
        - LoS_measurements - a Dictionary of lines which contains a 
                                    Dictionary of points, where Point name is 
                                    a key and measured values are triplet tuple"""
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

# Checking point names for typos and misspells and creating the list of measured
# points in lines sorted based on the config file definition:
sorted_LoS_measurements = {}
for line in LoS_measurements:
    if (cg.Print_typos) and (line not in cg.Lines_of_sight):
# printing which lines are in measurements input but are not in the default
# naming either due to typo or just simply missing in the nominal LoS decription
        print("Line %s was measured, but not expected." % (line))
    else:
        line_points_sorted = []
        for point in cg.Lines_of_sight[line]:
            line_points_sorted.append(point)
        sorted_LoS_measurements[line] = tuple(line_points_sorted)
        del line_points_sorted    
    for point in LoS_measurements[line]:
        if (cg.Print_typos) and (point not in Nominal_coords):
            print("Measured point with name %s in %s is not in the Nominal Coordinate file." % (point, line))
    del line, point
for point in Nominal_coords.keys():
    if (cg.Print_typos) and (point not in LoS_measurements):
        print("Point %s was not measured in any line." % (point))
    del point

#for line in LoS_measurements:
#    print(list(LoS_measurements[line].keys()))
print('End')