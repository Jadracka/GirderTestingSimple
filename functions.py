# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 12:53:33 2021

@author: jbarker
"""
import math as m
from numpy import pi
import numpy as np
import re
import Helmert3Dtransform as helmt

"""
______                _   _                 
|  ___|              | | (_)                
| |_ _   _ _ __   ___| |_ _  ___  _ __  ___ 
|  _| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
| | | |_| | | | | (__| |_| | (_) | | | \__ \
\_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
"""
def slope_distance(Point_From,Point_To):
	"""Function calculates slope distance from point 1 to point 2
	   from given Tupples (X1, Y1, Z1) and (X2, Y2, Z2)"""
	dX = Point_To[0] - Point_From[0]
	dY = Point_To[1] - Point_From[1]
	dZ = Point_To[2] - Point_From[2]
	sd = m.sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2))
	return sd

def horizontal_distance(Point_From,Point_To):
    dX = Point_To[0] - Point_From[0]
    dY = Point_To[1] - Point_From[1]
    hd = m.sqrt(pow(dX,2)+pow(dY,2))
    return hd

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

def arctang(angle):
    """Function takes angle in gons and calculates cosinus"""
    result = m.atan(gon2rad(angle))
    return result

def polar2cart3Dgon(Point):
    return (
         Point[0] * sing(Point[2]) * cosg(Point[1]),
         Point[0] * sing(Point[2]) * sing(Point[1]),
         Point[0] * cosg(Point[2])
    )

def cart2polal3Dgon(Point):
    return (
         m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2) + m.pow(Point[2],2)),
         rad2gon(
            arctang(m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2))/Point[2])),
         rad2gon(arctang(Point[1]/Point[0]))
    )
         
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

def StDev_sys_ppm(Value,StDev_tuple):
    return float(StDev_tuple[0] + Value * StDev_tuple[1]/1000)

def Helmert_calc_for_PolMeas(From,To):
    Transformed_From = {}
    for instrument in From:
        x = helmt.Helmert_transform(From[instrument],To)
        Transformed_From[instrument] = helmt.Transformation(x,From[instrument])
        Transformed_From[instrument][instrument] = (tuple(x[:3]))
    return Transformed_From

def StDev_XYZ_from_Polar(Point, StDev_S, StDev_Hz, StDev_Z):
    '''takes in Point measured in Polar coordinates and outputs tuple of 
    standard deviations for the XYZ components using config.py's values for 
    standard deviations of the measurements.'''
    S = Point[0]/1000
    Hz = Point[1]
    Z = Point[2]
    StDev_X = m.sqrt(m.pow(StDev_S,2) * m.pow(cosg(Hz) * sing(Z),2)\
              + m.pow(StDev_Hz,2) * m.pow(-S * sing(Hz) * sing(Z),2)\
              + m.pow(StDev_Z,2) * m.pow(S * cosg(Hz) * cosg(Z),2))
    StDev_Y = m.sqrt(m.pow(StDev_S,2) * m.pow(sing(Hz) * sing(Z),2)\
              + m.pow(StDev_Hz,2) * m.pow(S * cosg(Hz) * sing(Z),2)\
              + m.pow(StDev_Z,2) * m.pow(S * sing(Hz) * cosg(Z),2))
    StDev_Zz = m.sqrt(m.pow(StDev_S,2) * m.pow(cosg(Z),2)\
              + m.pow(StDev_Z,2) * m.pow(-S * sing(Z),2))
    return (StDev_X,StDev_Y,StDev_Zz)

def StDev_distance(Point_From, Point_To, StDevXYZ_From, StDevXYZ_To):
    sd = slope_distance(Point_From, Point_To)
    StDev_S = m.sqrt(m.pow(StDevXYZ_From[0],2) * m.pow((Point_To[0] \
                                                    - Point_From[0])/sd,2)
                 + m.pow(StDevXYZ_To[0],2) * m.pow((Point_To[0] \
                                                    - Point_From[0])/-sd,2)
                 + m.pow(StDevXYZ_From[1],2) * m.pow((Point_To[1] \
                                                    - Point_From[1])/sd,2)
                 + m.pow(StDevXYZ_To[1],2) * m.pow((Point_To[1] \
                                                    - Point_From[1])/-sd,2)
                 + m.pow(StDevXYZ_From[2],2) * m.pow((Point_To[2] \
                                                    - Point_From[2])/sd,2)
                 + m.pow(StDevXYZ_To[2],2) * m.pow((Point_To[2] \
                                                    - Point_From[2])/-sd,2)
        )
    return StDev_S

def ParD_Hz(PointTo, PointFrom):
    # This function returns derivatives of the horizontal angle with respect to
    # all unknowns, X, Y, Z, O (orientation) for point
    # For PointFrom add '-' in front of dX and dY
    dX = -(PointTo[1] - PointFrom[1]) / (pow(PointTo[0] - PointFrom[0],2) \
            + pow(PointTo[1] - PointFrom[1],2))
    dY = (PointTo[0] - PointFrom[0]) / (pow(PointTo[0] - \
             PointFrom[0],2) + pow(PointTo[1] - PointFrom[1],2))
    dZ = 0
    dO = -1
    return dX, dY, dZ, dO

def ParD_V(PointTo, PointFrom):
    # This function returns derivatives of the zenith angle with respect to
    # all unknowns, X, Y, Z, O (orientation) for PointTo
    # For PointFrom add '-' in front of dX, dY and dZ
    dist_squared = pow(PointTo[0] - PointFrom[0],2) 
    + pow(PointTo[1] - PointFrom[1],2) + pow(PointTo[2] - PointFrom[2],2)
    h_distance = horizontal_distance(PointTo, PointFrom)
    dX = ((PointTo[0]-PointFrom[0]) * (PointTo[2]-PointFrom[2])) \
        / (dist_squared * h_distance)
    dY = ((PointTo[1]-PointFrom[1]) * (PointTo[2]-PointFrom[2])) \
        / (dist_squared * h_distance)
    dZ = - h_distance / dist_squared
    return dX, dY, dZ


def ParD_sd(PointTo, PointFrom):
    # For PointToFrom add '-' in front of dX, dY and dZ
    dist = slope_distance(PointTo, PointFrom)
    if dist == 0:
        print('ParD_sd: PointTo is identical with PointFrom.')
    dX = (PointTo[0]-PointFrom[0]) / dist
    dY = (PointTo[1]-PointFrom[1]) / dist
    dZ = (PointTo[2]-PointFrom[2]) / dist
    return dX, dY, dZ
    
def horizontal_angle_from_Coords(PointTo,PointFrom):
    Hz = m.atan2(PointTo[1]-PointFrom[1],PointTo[0]-PointFrom[0])
    return Hz

def vertical_angle_from_Coords(PointTo,PointFrom):
    dist = slope_distance(PointFrom, PointTo)
    V = m.acos((PointTo[2]-PointFrom[2])/dist)
    return V

def find_unknowns(Dict_of_measurements):
    # Dictionary of transformed measurements is input (must include instrument
    # station so it gets counted to unknowns!)
    unknowns = []
    unknown_points = []
    unknown_instrument_stations = []
    unknown_instrument_orientations = []
    for instrument in Dict_of_measurements:
        new = list(Dict_of_measurements[instrument].keys())
        new.pop(new.index(instrument))
        unknown_points = list(set(unknown_points + new))
        unknown_instrument_stations.append(instrument)
        unknown_instrument_orientations.append(('Ori_' + instrument))
    unknowns = unknown_points + unknown_instrument_stations +\
        unknown_instrument_orientations
    number_of_instruments = len(unknown_instrument_stations)
    number_of_unknowns = len(unknown_points)*3 + number_of_instruments*4
    return unknowns, number_of_unknowns, number_of_instruments


def merge_measured_coordinates(Dictionary):
    # Takes in dictionary with transformed cartesian coordinates with structure
    # Dictionary: Instrument: dictionary of measured points
    result = {}
    if len(Dictionary.keys()) == 2:
        keys = list(Dictionary.keys())
        result = dict(Dictionary[keys[0]].items() | Dictionary[keys[1]].items())
        intersection = Dictionary[keys[0]].keys() & Dictionary[keys[1]].keys()
        for point in intersection:
            point_a = np.array(Dictionary[keys[0]][point])
            point_b = np.array(Dictionary[keys[1]][point])
            result[point] = tuple((point_a + point_b)/2)
    elif len(Dictionary.keys()) == 1:
        result = Dictionary.values()
    else:
        for d in Dictionary.values():
            for PointID, Point in d.items():
                if PointID in result:
                    result[PointID].append(Point)
                else:
                    result[PointID] = [Point,]
        for PointID in result:
            result[PointID] = tuple(sum(ele) / len(result[
                                    PointID]) for ele in zip(*result[PointID]))

    return result

def filling_X(Aproximates, unknowns, count_unknowns, count_instruments):
    X = np.zeros(count_unknowns)
    XHR = []
    for i, unknown in enumerate(unknowns[:-count_instruments]):
        iii = 3*i 
        X[iii:iii+3] = np.array(Aproximates[unknown])
        # Human readible version of X
        XHR.append(('X ' + unknown))
        XHR.append(('Y ' + unknown))
        XHR.append(('Z ' + unknown))
    for unknown in unknowns[-count_instruments:]:
        XHR.append(unknown)
    return X, XHR

def filling_Aproximates(unknowns, X_vector, count_instruments):
    updated_Aproximates = dict()
    for i, point in enumerate(unknowns[:-count_instruments]):
        iii = 3*i
        updated_Aproximates[point] = tuple(X_vector[iii:iii+3])
    return updated_Aproximates