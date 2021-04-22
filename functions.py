# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 12:53:33 2021

@author: jbarker
"""
import math as m
from numpy import pi
import re

"""
______                _   _                 
|  ___|              | | (_)                
| |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
|  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
| | | |_| | | | | (__| |_| | (_) | | | \__ \
\_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
"""
def slope_distance(X1,Y1,Z1,X2,Y2,Z2):
	"""Function calculates slope distance from point 1 to point 2
	   from given X1, Y1, Z1, X2, Y2, Z2
       THIS IS MY FIRST OVERLOADED FUNCTION!!!!"""
	dX = X2 - X1
	dY = Y2 - Y1
	dZ = Z2 - Z1
	sd = m.sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2))
	return sd

def slope_distance(Point_From,Point_To):
	"""Function calculates slope distance from point 1 to point 2
	   from given Tupples (X1, Y1, Z1) and (X2, Y2, Z2)"""
	dX = Point_To[0] - Point_From[0]
	dY = Point_To[1] - Point_From[1]
	dZ = Point_To[2] - Point_From[2]
	sd = m.sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2))
	return sd

def gon2rad(gons):
    """Function takes an angle in gons, transforms to a float and converts
       to radians"""
    rads = float(gons)*pi/200
    return rads

def rad2gon(rads):
    """Function takes an angle in gons, transforms to a float and converts
       to radians"""
    gons = float(rads)*200/pi
    return gons

def cosg(angle):
    """Function takes angle in gons and calculates cosinus"""
    result = m.cos(gon2rad(angle))
    return result

def sing(angle):
    """Function takes angle in gons and calculates sinus"""
    result = m.sin(gon2rad(angle))
    return result

def polar2cart3D(S, Hz, Z):
    return [
         S * sing(Hz) * cosg(Z),
         S * sing(Hz) * sing(Z),
         S * cosg(Hz)
    ]
    
def polar2cart3D(PointID):
    return [
         PointID[0] * sing(PointID[1]) * cosg(PointID[2]),
         PointID[0] * sing(PointID[1]) * sing(PointID[2]),
         PointID[0] * cosg(PointID[1])
    ]

def Measurements_read_in(Meas_filename):
    Meas_file = open(Meas_filename,'r')
    Measurements = {}
    for row in Meas_file.readlines():
        """
        - Reads in and parses Meas_file. 
        - The format is string Line name 'space' 
          string Point name 'space' float Distance [mm] 'space' 
          float Hz angle [gon] 'space' float Z angle [gon].
        - Ignores point notes at the end.
        - It can handle multiple occurences of delimeters but not a 
          combination of them.
            - LoS_measurements_1 - a Dictionary of lines which contains a 
                                        Dictionary of points, where Point
                                        name is a key and measured values
                                        are triplet tuple"""
        words = re.split(';+|,+|\t+| +',row.strip())
        if words[0] not in Measurements.keys():
            Measurements[words[0]] = {}
            Measurements[words[0]][words[1]] =                         \
              (float(words[2]), float(words[3]), float(words[4]))
        else: 
            Measurements[words[0]][words[1]] =                         \
              (float(words[2]), float(words[3]), float(words[4]))
    del words, row
    Meas_file.close()
    return Measurements
    
def Coords_read_in(Coords_file_name):
    Coords_file = open(Coords_file_name,'r')
    Coords = {}
    for line in Coords_file.readlines():
        """Reads in and parses coordinates file. Ignores point notes at the end.
           It can handle multiple occurences of delimeters,
           but not a combination of them."""
        words = re.split(';+|,+|\t+| +',line.strip())
        Coords[words[0]] = (float(words[1]), 
                                    float(words[2]), 
                                    float(words[3]))
        del line, words
    Coords_file.close()
    return Coords