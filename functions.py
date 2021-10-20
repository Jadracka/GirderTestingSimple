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
import config as cg


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
	Sd = m.sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2))
	return Sd

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

def atang(arg):
    """Function takes argument and calculates arcus tangens"""
    result = rad2gon(m.atan(arg))
    return result

def atan2g(dx,dy):
    """Function takes two arguments and calculates arcut tangens"""
    result = rad2gon(m.atan2(dy,dx))
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
         atan2g(m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2)),Point[2]),
         atan2g(Point[1],Point[0])
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
    
def Polar_2F_meas_read_in(Meas_filename,
                          Sd_StDev=(0.010,0), Hz_StDev=0.15, V_StDev=0.15):
    """Standard deviations for standard measurements with Leica AT960/AT930)"""
    Meas_file = open(Meas_filename, 'r')
    max_diff_Hz = (Hz_StDev * m.sqrt(2))/1000
    max_index_err = (V_StDev * m.sqrt(2))/1000
    Measurements = {}
    Diffs = {}
    for row in Meas_file.readlines():
        """
        - Reads in and parses Meas_file.
        - The format is string Group name (should be instrument name) 'space'
          string Point name 'space' float Hz angle [gon] 'space'
          float Z angle [gon] 'space' float Distance [mm].
        - Ignores point notes at the end.
        - It can handle multiple occurences of delimeters but not a
          combination of them.
        - handles possible 2 face measurements and their statistics"""
        words = re.split(';+|,+|\t+| +', row.strip())
        if words[0] not in Measurements.keys():
            Measurements[words[0]] = {}
            Measurements[words[0]][words[1]] =                         \
                (float(words[4]), float(words[2]), float(words[3]), '1F')
            Diffs[words[0]] = {'Hz': (), 'V': (), 'Sd': ()}
        else:
            if words[1] not in Measurements[words[0]].keys():
                Measurements[words[0]][words[1]] =                         \
                (float(words[4]), float(words[2]), float(words[3]), '1F')
            else:
                Meas1 = Measurements[words[0]][words[1]][:3]
                Meas2 = (float(words[4]), float(words[2]), float(words[3]))
                # checking difference of distance measurements for both faces
                new_Sd = (Meas1[0] + Meas2[0])/2
                max_diff_Sd = m.sqrt(2) * StDev_sys_ppm(new_Sd, Sd_StDev)
                diff_Sd = Meas2[0] - Meas1[0]
                Diffs[words[0]]['Sd'] = Diffs[words[0]]['Sd'] + (diff_Sd/2,)
                if (diff_Sd > max_diff_Sd) and cg.Print_2F_checks:
                    print("Point: %s, measured by %s fails 2Face check in "
                          "distance. Maximum difference is %1.4f mm and "
                          "measured difference is %1.4f mm.\n"
                          % (words[1], words[0], max_diff_Sd, diff_Sd))
                # assigning face one and face two:
                if Meas2[2] >= 200 and Meas1[2] < 200:
                    Face1 = Meas1
                    Face2 = Meas2
                elif Meas1[2] >= 200 and Meas2[2] < 200:
                    Face1 = Meas2
                    Face2 = Meas1
                else:
                    print("First and second face measurements check failed,"
                          "vertical angles don't make sense for point %s "
                          "measured by %s.\n" % (words[1], words[0]))
                # averaging V 1,2
                index_err = -(Face1[2] + Face2[2] - 400)/2
                Diffs[words[0]]['V'] = Diffs[words[0]]['V'] + (index_err,)
                if (abs(index_err) > max_index_err) and cg.Print_2F_checks:
                    print("Index error is higher than expected. Point %s, meas"
                          "ured by %s has an index error of %1.4f, and maximum"
                          " error is %1.4f.\n"
                          % (words[1], words[0], index_err, max_index_err))
                new_V = Face1[2] + index_err
                # averaging Hz
                if Face1[1] > 200:
                    # for angles above 200g is better to switch to -200,200,
                    # so there is no 400 overflow problem:
                    Face2Hz = Face2[1] + 200
                    conv_Face1 = Face1[1] - 400
                    conv_Face2 = Face2Hz - 400
                    diff_Hz = conv_Face2 - conv_Face1
                    Diffs[words[0]]['Hz'] = Diffs[words[0]]['Hz'] + (diff_Hz/2,)
                    if (abs(diff_Hz) > max_diff_Hz) and cg.Print_2F_checks:
                        print("Two face difference in Hz angle is higher than "
                              "expected. Point %s, measured by %s has an index"
                              " error of %1.4f, and maximum error is %1.4f.\n"
                              % (words[1], words[0], index_err, max_index_err))
                    new_Hz = ((conv_Face1 + conv_Face2)/2) + 400
                else:
                    # Otherwise just simple averaging is ok:
                    Face2Hz = Face2[1] - 200
                    new_Hz = (Face1[1] + Face2Hz)/2
                    diff_Hz = Face2Hz - Face1[1]
                    Diffs[words[0]]['Hz'] = Diffs[words[0]]['Hz'] + (diff_Hz/2,)
                    if (abs(diff_Hz) > max_diff_Hz) and cg.Print_2F_checks:
                        print("Two face difference in Hz angle is higher than "
                              "expected. Point %s, measured by instrument %s "
                              "has an index error of %1.4f, and maximum error "
                              "is %1.4f.\n"
                              % (words[1], words[0], index_err, max_index_err))
# Now updating the original dictionary entry with averaged values:
                Measurements[words[0]].update({words[1]: (new_Sd, new_Hz,
                            new_V, max_diff_Sd/m.sqrt(2), Hz_StDev, V_StDev)})
    for instrument in Measurements:
        Diffs[instrument]['Corr_median'] = (np.nanmedian(Diffs[instrument]['Sd']),
        np.nanmedian(Diffs[instrument]['Hz']),np.nanmedian(Diffs[instrument]['V']))
        if np.isnan(Diffs[instrument]['Corr_median']).all():
            Diffs[instrument]['Corr_median'] = (0,0,0)
        for point in Measurements[instrument]:
            if '1F' in Measurements[instrument][point]:
                # Calculating StDevs for 1F Measurements
                Sd_StDev_p = m.sqrt(2) * (StDev_sys_ppm(
                                Measurements[instrument][point][0],Sd_StDev))
                """ tuple(map(sum, zip(a,b)) returns a tuple with element-wise
                addition.
                Here adding median corrections to the original measured values,
                for only 1F measured """
                Measurements[instrument].update({point: tuple(map(sum,
                    zip(Measurements[instrument][point][:3],
                    Diffs[instrument]['Corr_median']))) + (
                    Sd_StDev_p, max_diff_Hz, max_index_err)})

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
                                    -float(words[2]), 
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
    Sd = slope_distance(Point_From, Point_To)
    StDev_S = m.sqrt(m.pow(StDevXYZ_From[0],2) * m.pow((Point_To[0] \
                                                    - Point_From[0])/Sd,2)
                 + m.pow(StDevXYZ_To[0],2) * m.pow((Point_To[0] \
                                                    - Point_From[0])/-Sd,2)
                 + m.pow(StDevXYZ_From[1],2) * m.pow((Point_To[1] \
                                                    - Point_From[1])/Sd,2)
                 + m.pow(StDevXYZ_To[1],2) * m.pow((Point_To[1] \
                                                    - Point_From[1])/-Sd,2)
                 + m.pow(StDevXYZ_From[2],2) * m.pow((Point_To[2] \
                                                    - Point_From[2])/Sd,2)
                 + m.pow(StDevXYZ_To[2],2) * m.pow((Point_To[2] \
                                                    - Point_From[2])/-Sd,2)
        )
    return StDev_S

def ParD_Hz(PointTo, PointFrom):
    # This function returns derivatives of the horizontal angle with respect to
    # all unknowns, X, Y, Z, O (orientation) for point
    # For PointFrom add '-' in front of dX and dY
#    dX =-(PointTo[1] - PointFrom[1]) / (pow(PointTo[0] - PointFrom[0],2) \
#            + pow(PointTo[1] - PointFrom[1],2))
#    dY =+(PointTo[0] - PointFrom[0]) / (pow(PointTo[0] - PointFrom[0],2) \
#            + pow(PointTo[1] - PointFrom[1],2))
#    dZ = 0
#    dO = -1
    dX =+(PointTo[1] - PointFrom[1]) / (pow(PointTo[0] - PointFrom[0],2) \
            + pow(PointTo[1] - PointFrom[1],2))
    dY =-(PointTo[0] - PointFrom[0]) / (pow(PointTo[0] - PointFrom[0],2) \
            + pow(PointTo[1] - PointFrom[1],2))
    dZ = 0
    dO = 1
    return dX, dY, dZ, dO

def ParD_V(PointTo, PointFrom):
    # This function returns derivatives of the zenith angle with respect to
    # all unknowns, X, Y, Z, O (orientation) for PointTo
    # For PointFrom add '-' in front of dX, dY and dZ
    dist_squared = pow(PointTo[0] - PointFrom[0],2) \
    + pow(PointTo[1] - PointFrom[1],2) + pow(PointTo[2] - PointFrom[2],2)
    h_distance = horizontal_distance(PointTo, PointFrom)
#    dX = -((PointTo[0]-PointFrom[0]) * (PointTo[2]-PointFrom[2])) \
#        / (dist_squared * h_distance)
#    dY = ((PointTo[1]-PointFrom[1]) * (PointTo[2]-PointFrom[2])) \
#        / (dist_squared * h_distance)
#    dZ = - h_distance / dist_squared
    dX = -((PointTo[0]-PointFrom[0]) * (PointTo[2]-PointFrom[2])) \
        / (dist_squared * h_distance)
    dY = -((PointTo[1]-PointFrom[1]) * (PointTo[2]-PointFrom[2])) \
        / (dist_squared * h_distance)
    dZ = + h_distance / dist_squared
    return dX, dY, dZ


def ParD_Sd(PointTo, PointFrom):
    # For PointToFrom add '-' in front of dX, dY and dZ
    dist = slope_distance(PointTo, PointFrom)
    if dist == 0:
        print('ParD_Sd: PointTo is identical with PointFrom.')
#    dX = (PointTo[0]-PointFrom[0]) / dist
#    dY = (PointTo[1]-PointFrom[1]) / dist
#    dZ = (PointTo[2]-PointFrom[2]) / dist
    dX = -(PointTo[0]-PointFrom[0]) / dist
    dY = -(PointTo[1]-PointFrom[1]) / dist
    dZ = -(PointTo[2]-PointFrom[2]) / dist
    return dX, dY, dZ
    
def horizontal_angle_from_Coords(PointTo,PointFrom):
    Hz = m.atan2(PointFrom[1]-PointTo[1],PointFrom[0]-PointTo[0])
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
    instruments = []
    for instrument in Dict_of_measurements:
        instruments.append(instrument)
        new = list(Dict_of_measurements[instrument].keys())
        new.pop(new.index(instrument))
        unknown_points = list(set(unknown_points + new))
        unknown_instrument_stations.append(instrument)
        unknown_instrument_orientations.append(('Ori_' + instrument))
    unknowns = unknown_points + unknown_instrument_stations +\
        unknown_instrument_orientations
    count_instruments = len(unknown_instrument_stations)
    number_of_unknowns = len(unknown_points)*3 + count_instruments*4
    return unknowns, number_of_unknowns, instruments, count_instruments


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
### Hier muss noch die Berechnung des Abrisses für die O-Unbek rein!        
        # Human readible version of X
        XHR.append(('X ' + unknown))
        XHR.append(('Y ' + unknown))
        XHR.append(('Z ' + unknown))
    for unknown in unknowns[-count_instruments:]:
        XHR.append(unknown)
    return X, XHR

def filling_Aproximates(unknowns, X_vector, instruments):
    updated_Aproximates = dict()
    inst_count = len(instruments)
    for i, point in enumerate(unknowns[:-inst_count]):
        iii = 3*i
        updated_Aproximates[point] = tuple(X_vector[iii:iii+3])
    for i,instrument in enumerate(instruments):
        updated_Aproximates['Ori_' + instrument] = X_vector[-inst_count+i]
    return updated_Aproximates