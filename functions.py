# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 12:53:33 2021

@author: jbarker
"""
import math as m
from numpy import pi
import numpy as np
#import scipy as sp
from scipy import linalg
from numpy.linalg import inv
import re
import Helmert3Dtransform as helmt
import config as cg
from angle import Angle as a
import sys

def pretty_print (x):
    for i in x:
        print("{:7.2f} ".format(i),end='')
    print()
    return()

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

def slope_distance_6DoF(Aproximates,point,inst):
	Rx, Ry, Rz = Aproximates['Ori_'+inst]
	dX = Aproximates[point][0] - Aproximates[inst][0]
	dY = Aproximates[point][1] - Aproximates[inst][1]
	dZ = Aproximates[point][2] - Aproximates[inst][2]
	Rxc = m.cos(Rx)
	Rxs = m.sin(Rx)
	Ryc = m.cos(Ry)
	Rys = m.sin(Ry)
	Rzc = m.cos(Rz)
	Rzs = m.sin(Rz)
	Sd = m.sqrt(pow((Ryc*Rzc*dX + Ryc*Rzs*dY + Rys*dZ),2)+\
               pow((-Rxs*Rys*Rzc - Rxc*Rzs)*dX + \
						(-Rxs*Rys*Rzs + Rxc*Rzc)*dY + (Rxs*Ryc*dZ),2) + \
               pow((-Rxc*Rys*Rzc + Rxs*Rzs)*dX + (-Rxc*Rys*Rzs - Rxs*Rzc)*dY +\
				      (Rxc*Ryc*dZ),2))
	return Sd

def vertical_angle_6DoF(Aproximates,point,inst):
	Rx, Ry, Rz = Aproximates['Ori_'+inst]
	dX = Aproximates[point][0] - Aproximates[inst][0]
	dY = Aproximates[point][1] - Aproximates[inst][1]
	dZ = Aproximates[point][2] - Aproximates[inst][2]
	Rxc = m.cos(Rx)
	Rxs = m.sin(Rx)
	Ryc = m.cos(Ry)
	Rys = m.sin(Ry)
	Rzc = m.cos(Rz)
	Rzs = m.sin(Rz)
	V = a(m.acos(((-Rxc*Rys*Rzc + Rxs*Rzs)*dX + (-Rxc*Rys*Rzs - Rxs*Rzc)*dY + \
			    (Rxc*Ryc*dZ))/ \
			    m.sqrt(pow((Ryc*Rzc*dX + Ryc*Rzs*dY + Rys*dZ),2)+\
	                 pow((-Rxs*Rys*Rzc - Rxc*Rzs)*dX + \
							  (-Rxs*Rys*Rzs + Rxc*Rzc)*dY + (Rxs*Ryc*dZ),2) + \
					    pow((-Rxc*Rys*Rzc + Rxs*Rzs)*dX + \
				           (-Rxc*Rys*Rzs - Rxs*Rzc)*dY + (Rxc*Ryc*dZ),2))),
			a.T_RAD,True).angle
	return V

def horizontal_angle_6DoF(Aproximates,point,inst):
	Rx, Ry, Rz = Aproximates['Ori_'+inst]
	dX = Aproximates[point][0] - Aproximates[inst][0]
	dY = Aproximates[point][1] - Aproximates[inst][1]
	dZ = Aproximates[point][2] - Aproximates[inst][2]
	Rxc = m.cos(Rx)
	Rxs = m.sin(Rx)
	Ryc = m.cos(Ry)
	Rys = m.sin(Ry)
	Rzc = m.cos(Rz)
	Rzs = m.sin(Rz)
	Hz = a(m.atan2((-Rxs*Rys*Rzc - Rxc*Rzs)*dX + (-Rxs*Rys*Rzs + Rxc*Rzc)*dY + \
			       (Rxs*Ryc*dZ),(Ryc*Rzc*dX + Ryc*Rzs*dY + Rys*dZ)),a.T_RAD,True).angle
	return Hz

def horizontal_distance(Point_From,Point_To):
    dX = Point_To[0] - Point_From[0]
    dY = Point_To[1] - Point_From[1]
    hd = m.sqrt(pow(dX,2)+pow(dY,2))
    return hd

def gon2rad(gons):
    """Function takes an angle in gons, transforms to a longdouble and converts
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

def StDev_sys_ppm(Value,StDev_tuple):
    return float(StDev_tuple[0] + Value * StDev_tuple[1]/1000)

def polar2cart3Dgon(Point):
    return (
         Point['Sd'] * sing(Point['V']) * cosg(Point['Hz']),
         Point['Sd'] * sing(Point['V']) * sing(Point['Hz']),
         Point['Sd'] * cosg(Point['V'])
    )

def polar2cart3Drad(r, theta, phi):
    X = r * m.cos(theta) * m.sin(phi)
    Y = r * m.sin(theta) * m.sin(phi)
    Z = r * m.cos(phi)
    return X,Y,Z

def cart2polar3Dgon(Point):
    return (
         m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2) + m.pow(Point[2],2)),
         atan2g(m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2)),Point[2]),
         atan2g(Point[1],Point[0])
         )

def cart2polar3Drad(Point):
    return (
         m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2) + m.pow(Point[2],2)),
         m.atan2(m.sqrt(m.pow(Point[0],2) + m.pow(Point[1],2)),Point[2]),
         m.atan2(Point[1],Point[0])
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
              (float(words[4])/1000, float(words[2]), float(words[3]))
        else: 
            Measurements[words[0]][words[1]] =                         \
              (float(words[4])/1000, float(words[2]), float(words[3]))
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
        - handles possible 2 face measurements and their statistics
        - outputs distances in meters"""
        words = re.split(';+|,+|\t+| +', row.strip())
        if words[0] not in Measurements.keys():
            Measurements[words[0]] = {}
            Measurements[words[0]][words[1]] = {'Sd': float(words[4])/1000, 
                        'Hz': float(words[2]), 'V': float(words[3]), 'Face':1}
            Diffs[words[0]] = {'Hz': (), 'V': (), 'Sd': ()}
        else:
            if words[1] not in Measurements[words[0]].keys():
                Measurements[words[0]][words[1]] = {'Sd':float(words[4])/1000, 
                 'Hz': float(words[2]), 'V': float(words[3]), 'Face':1}
            else:
                Meas1 = Measurements[words[0]][words[1]]
                Meas2 = {'Sd':float(words[4])/1000, 'Hz': float(words[2]), 
                         'V':float(words[3])}
                # checking difference of distance measurements for both faces
                new_Sd = (Meas1['Sd'] + Meas2['Sd'])/2
                max_diff_Sd = (m.sqrt(2) * StDev_sys_ppm(new_Sd, Sd_StDev))/1000
                diff_Sd = Meas2['Sd'] - Meas1['Sd']
                Diffs[words[0]]['Sd'] = Diffs[words[0]]['Sd'] + (diff_Sd/2,)
                if (abs(diff_Sd) > max_diff_Sd) and cg.Print_2F_checks:
                    print("Point: %s, measured by %s fails 2Face check in "
                          "distance. Maximum difference is %1.4f m and "
                          "measured difference is %1.4f m.\n"
                          % (words[1], words[0], max_diff_Sd, abs(diff_Sd)))
                # assigning face one and face two:
                if Meas2['V'] >= 200 and Meas1['V'] < 200:
                    Face1 = Meas1
                    Face2 = Meas2
                elif Meas1['V'] >= 200 and Meas2['V'] < 200:
                    Face1 = Meas2
                    Face2 = Meas1
                else:
                    print("First and second face measurements check failed,"
                          "vertical angles don't make sense for point %s "
                          "measured by %s.\n" % (words[1], words[0]))
                # averaging V 1,2
                index_err = -(Face1['V'] + Face2['V'] - 400)/2
                Diffs[words[0]]['V'] = Diffs[words[0]]['V'] + (index_err,)
                if (abs(index_err) > max_index_err) and cg.Print_2F_checks:
                    print("Index error is higher than expected. Point %s, meas"
                          "ured by %s has an index error of %1.4f, and maximum"
                          " error is %1.4f.\n"
                          % (words[1], words[0], index_err, max_index_err))
                new_V = Face1['V'] + index_err
                # averaging Hz
                if Face1['Hz'] > 200:
                    # for angles above 200g is better to switch to -200,200,
                    # so there is no 400 overflow problem:
                    Face2Hz = Face2['Hz'] + 200
                    conv_Face1 = Face1['Hz'] - 400
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
                    Face2Hz = Face2['Hz'] - 200
                    new_Hz = (Face1['Hz'] + Face2Hz)/2
                    diff_Hz = Face2Hz - Face1['Hz']
                    Diffs[words[0]]['Hz'] = Diffs[words[0]]['Hz'] + (diff_Hz/2,)
                    if (abs(diff_Hz) > max_diff_Hz) and cg.Print_2F_checks:
                        print("Two face difference in Hz angle is higher than "
                              "expected. Point %s, measured by instrument %s "
                              "has an index error of %1.4f, and maximum error "
                              "is %1.4f.\n"
                              % (words[1], words[0], index_err, max_index_err))
# Now updating the original dictionary entry with averaged values:
                Measurements[words[0]].update({words[1]: {'Sd': new_Sd, 
                            'Hz': new_Hz, 'V': new_V, 'Face':2,
                            'StDev_Sd': max_diff_Sd/m.sqrt(2), 
                            'StDev_Hz': Hz_StDev, 'StDev_V': V_StDev}})
    for instrument in Measurements:
        if len(Diffs[instrument]['Sd']) > 0 and len(
            Diffs[instrument]['Hz']) > 0 and len(Diffs[instrument]['V']) > 0:
			
            Diffs[instrument]['Corr_median'] = {
                    'Med_Sd': np.median(Diffs[instrument]['Sd']),
                    'Med_Hz': np.median(Diffs[instrument]['Hz']), 
                    'Med_V': np.median(Diffs[instrument]['V'])}
        else:
            Diffs[instrument]['Corr_median'] = {'Med_Sd': 0,'Med_Hz': 0, 
                                                'Med_V': 0}
        for point in Measurements[instrument]:
            if Measurements[instrument][point]['Face'] == 1:
                # Calculating StDevs for 1F Measurements
                Sd_StDev_p = (m.sqrt(2) * (StDev_sys_ppm(
                        Measurements[instrument][point]['Sd'],Sd_StDev)))/1000
                """ tuple(map(sum, zip(a,b)) returns a tuple with element-wise
                addition.
                Here adding median corrections to the original measured values,
                for only 1F measured """
                Measurements[instrument].update({point: {
                        'Hz': Measurements[instrument][point]['Hz'] + 
                        Diffs[instrument]['Corr_median']['Med_Hz'],
                        'Sd': Measurements[instrument][point]['Sd'] + 
                        Diffs[instrument]['Corr_median']['Med_Sd'],
                        'V': Measurements[instrument][point]['V'] + 
                        Diffs[instrument]['Corr_median']['Med_V'],
                        'StDev_Sd': Sd_StDev_p,'StDev_Hz': Hz_StDev*m.sqrt(2), 
                        'StDev_V':V_StDev*m.sqrt(2)}})

    del words, row
    Meas_file.close()
    return Measurements

def Count_meas_types(d, meas_type):
    c = int(meas_type in d)
    for v in d.values():
        if isinstance(v, dict):
            c += Count_meas_types(v, meas_type)
    return c

#counter_Hz = 0
#counter_Sd = 0
#counter_V = 0
#for instrument in Measurements:
#    for point in Measurements[instrument]:
#        for types in Measurements[instrument][point].keys():
#            if types == 'Hz':
#                counter_Hz +=1
#            if types == 'Sd':
#                counter_Sd +=1
#            if types == 'V':
#                counter_V +=1

           
def Coords_read_in(Coords_file_name):
    Coords_file = open(Coords_file_name,'r')
    Coords = {}
    for line in Coords_file.readlines():
        """Reads in and parses coordinates file. Ignores point notes at the end.
           It can handle multiple occurences of delimeters,
           but not a combination of them."""
        words = re.split(';+|,+|\t+| +',line.strip())
        Coords[words[0]] = (float(words[1])/1000, 
                            -float(words[2])/1000, 
                            float(words[3])/1000)
        del line, words
    Coords_file.close()
    return Coords

def Helmert_calc_for_PolMeas(From,To):
    Transformed_From = {}
    Trans_par = {}
    for instrument in From:
        x = helmt.Helmert_transform(From[instrument],To)
        Trans_par[instrument] = tuple(x)
        Transformed_From[instrument] = helmt.Transformation(x,From[instrument])
        Transformed_From[instrument][instrument] = (tuple(x[:3]))
    return Transformed_From, Trans_par

def ParD_Hz(PointTo, PointFrom):
    # This function returns derivatives of the horizontal angle with respect to
    # all unknowns, X, Y, Z, O (orientation) for point
    # For PointFrom add '-' in front of dX and dY
    dX = (PointTo[1] - PointFrom[1]) / (pow(PointTo[0] - PointFrom[0],2) \
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
    dX = -(PointTo[0]-PointFrom[0]) / dist
    dY = -(PointTo[1]-PointFrom[1]) / dist
    dZ = -(PointTo[2]-PointFrom[2]) / dist
    return dX, dY, dZ

def Sd_6Dof(dX, dY, dZ, Rxc, Rxs, Ryc, Rys, Rzc, Rzs):
    d_Sd = m.sqrt(pow((Ryc*Rzc*dX + Ryc*Rzs*dY + Rys*dZ),2)+\
                  pow((-Rxs*Rys*Rzc - Rxc*Rzs)*dX + (-Rxs*Rys*Rzs + Rxc*Rzc)*dY + (Rxs*Ryc*dZ),2) + \
                  pow((-Rxc*Rys*Rzc + Rxs*Rzs)*dX + (-Rxc*Rys*Rzs - Rxs*Rzc)*dY + (Rxc*Ryc*dZ),2))
    return d_Sd

def Sd_6Dof_noRot(dX, dY, dZ):
    d_Sd = m.sqrt(pow(dX,2)+pow(dY,2)+pow(dZ,2))
    return d_Sd

def Hz_6Dof(dX, dY, dZ, Rxc, Rxs, Ryc, Rys, Rzc, Rzs):
    d_Hz = m.atan2(((-Rxs*Rys*Rzc - Rxc*Rzs)*dX + (-Rxs*Rys*Rzs + Rxc*Rzc)*dY + (Rxs*Ryc)*dZ), (Ryc*Rzc*dX + Ryc*Rzs*dY + Rys*dZ))
    return d_Hz

def V_6Dof(dX, dY, dZ, Rxc, Rxs, Ryc, Rys, Rzc, Rzs):
    d_V = m.acos(((-Rxc*Rys*Rzc + Rxs*Rzs)*dX + (-Rxc*Rys*Rzs - Rxs*Rzc)*dY + (Rxc*Ryc)*dZ) / \
                m.sqrt(pow((Ryc*Rzc*dX + Ryc*Rzs*dY + Rys*dZ),2)+\
                       pow((-Rxs*Rys*Rzc - Rxc*Rzs)*dX + (-Rxs*Rys*Rzs + Rxc*Rzc)*dY + (Rxs*Ryc*dZ),2) + \
                       pow((-Rxc*Rys*Rzc + Rxs*Rzs)*dX + (-Rxc*Rys*Rzs - Rxs*Rzc)*dY + (Rxc*Ryc*dZ),2)))
    return d_V

def Par_6DoF_noRot(PointTo, PointFrom, Aproximates, epsilon):
    dX = Aproximates[PointTo][0] - Aproximates[PointFrom][0]
    dY = Aproximates[PointTo][1] - Aproximates[PointFrom][1]
    dZ = Aproximates[PointTo][2] - Aproximates[PointFrom][2]
    dXe = Aproximates[PointTo][0] - Aproximates[PointFrom][0] + epsilon
    dYe = Aproximates[PointTo][1] - Aproximates[PointFrom][1] + epsilon
    dZe = Aproximates[PointTo][2] - Aproximates[PointFrom][2] + epsilon
    dX_e = Aproximates[PointTo][0] - Aproximates[PointFrom][0] - epsilon
    dY_e = Aproximates[PointTo][1] - Aproximates[PointFrom][1] - epsilon
    dZ_e = Aproximates[PointTo][2] - Aproximates[PointFrom][2] - epsilon
    dX_Sd = (Sd_6Dof_noRot(dXe, dY, dZ) - Sd_6Dof_noRot(dX_e, dY, dZ))\
            /(2*epsilon)
    dY_Sd = (Sd_6Dof_noRot(dX, dYe, dZ) - Sd_6Dof_noRot(dX, dY_e, dZ))\
            /(2*epsilon)
    dZ_Sd = (Sd_6Dof_noRot(dX, dY, dZe) - Sd_6Dof_noRot(dX, dY, dZ_e))\
            /(2*epsilon)
    return dX_Sd, dY_Sd, dZ_Sd

def Par_6Dof(PointTo, Instrument, Aproximates, epsilon):
#    print(Instrument, PointTo)
    Rx, Ry, Rz = Aproximates['Ori_'+Instrument]
    dX = Aproximates[PointTo][0] - Aproximates[Instrument][0]
    dY = Aproximates[PointTo][1] - Aproximates[Instrument][1]
    dZ = Aproximates[PointTo][2] - Aproximates[Instrument][2]
    dXe = Aproximates[PointTo][0] - Aproximates[Instrument][0] + epsilon
    dYe = Aproximates[PointTo][1] - Aproximates[Instrument][1] + epsilon
    dZe = Aproximates[PointTo][2] - Aproximates[Instrument][2] + epsilon
    dX_e = Aproximates[PointTo][0] - Aproximates[Instrument][0] - epsilon
    dY_e = Aproximates[PointTo][1] - Aproximates[Instrument][1] - epsilon
    dZ_e = Aproximates[PointTo][2] - Aproximates[Instrument][2] - epsilon
    Rxc = m.cos(Rx)
    Rxs = m.sin(Rx)
    Ryc = m.cos(Ry)
    Rys = m.sin(Ry)
    Rzc = m.cos(Rz)
    Rzs = m.sin(Rz)
    Rxce = m.cos(Rx + epsilon)
    Rxse = m.sin(Rx + epsilon)
    Ryce = m.cos(Ry + epsilon)
    Ryse = m.sin(Ry + epsilon)
    Rzce = m.cos(Rz + epsilon)
    Rzse = m.sin(Rz + epsilon)
    Rxc_e = m.cos(Rx - epsilon)
    Rxs_e = m.sin(Rx - epsilon)
    Ryc_e = m.cos(Ry - epsilon)
    Rys_e = m.sin(Ry - epsilon)
    Rzc_e = m.cos(Rz - epsilon)
    Rzs_e = m.sin(Rz - epsilon)
#    print(Rxc, Rxs, Rys, Rzc, )
    dX_Sd =  (Sd_6Dof(dXe,  dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX_e, dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dY_Sd =  (Sd_6Dof(dX,   dYe,  dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY_e, dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dZ_Sd =  (Sd_6Dof(dX,   dY,   dZe,  Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY,   dZ_e, Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRx_Sd = (Sd_6Dof(dX,   dY,   dZ,   Rxce,  Rxse,  Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY,   dZ,   Rxc_e, Rxs_e, Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRy_Sd = (Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryce,  Ryse,  Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc_e, Rys_e, Rzc,   Rzs))   /(2*epsilon)
    dRz_Sd = (Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzce,  Rzse) -
              Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc_e, Rzs_e)) /(2*epsilon)
    dX_Hz =  (Hz_6Dof(dXe,  dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Hz_6Dof(dX_e, dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dY_Hz =  (Hz_6Dof(dX,   dYe,  dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Hz_6Dof(dX,   dY_e, dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dZ_Hz =  (Hz_6Dof(dX,   dY,   dZe,  Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Hz_6Dof(dX,   dY,   dZ_e, Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRx_Hz = (Hz_6Dof(dX,   dY,   dZ,   Rxce,  Rxse,  Ryc,   Rys,   Rzc,   Rzs)  -
              Hz_6Dof(dX,   dY,   dZ,   Rxc_e, Rxs_e, Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRy_Hz = (Hz_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryce,  Ryse,  Rzc,   Rzs)  -
              Hz_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc_e, Rys_e, Rzc,   Rzs))   /(2*epsilon)
    dRz_Hz = (Hz_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzce,  Rzse) -
              Hz_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc_e, Rzs_e)) /(2*epsilon)
    dX_V =  (V_6Dof(dXe,  dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
             V_6Dof(dX_e, dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dY_V =  (V_6Dof(dX,   dYe,  dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
             V_6Dof(dX,   dY_e, dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dZ_V =  (V_6Dof(dX,   dY,   dZe,  Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
             V_6Dof(dX,   dY,   dZ_e, Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRx_V = (V_6Dof(dX,   dY,   dZ,   Rxce,  Rxse,  Ryc,   Rys,   Rzc,   Rzs)  -
             V_6Dof(dX,   dY,   dZ,   Rxc_e, Rxs_e, Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRy_V = (V_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryce,  Ryse,  Rzc,   Rzs)  -
             V_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc_e, Rys_e, Rzc,   Rzs))   /(2*epsilon)
    dRz_V = (V_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzce,  Rzse) -
             V_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc_e, Rzs_e)) /(2*epsilon)
    
    return dX_Sd, dY_Sd, dZ_Sd, dRx_Sd, dRy_Sd, dRz_Sd, \
           dX_Hz, dY_Hz, dZ_Hz, dRx_Hz, dRy_Hz, dRz_Hz, \
           dX_V, dY_V, dZ_V, dRx_V, dRy_V, dRz_V

def Par_6Dof_IFM(PointTo, PointFrom, line, Aproximates, epsilon):
    Instrument = cg.Instruments_LoS[line]
    Rx, Ry, Rz = Aproximates['Ori_'+Instrument]
    dX = Aproximates[PointTo][0] - Aproximates[PointFrom][0]
    dY = Aproximates[PointTo][1] - Aproximates[PointFrom][1]
    dZ = Aproximates[PointTo][2] - Aproximates[PointFrom][2]
    dXe = Aproximates[PointTo][0] - Aproximates[PointFrom][0] + epsilon
    dYe = Aproximates[PointTo][1] - Aproximates[PointFrom][1] + epsilon
    dZe = Aproximates[PointTo][2] - Aproximates[PointFrom][2] + epsilon
    dX_e = Aproximates[PointTo][0] - Aproximates[PointFrom][0] - epsilon
    dY_e = Aproximates[PointTo][1] - Aproximates[PointFrom][1] - epsilon
    dZ_e = Aproximates[PointTo][2] - Aproximates[PointFrom][2] - epsilon
    Rxc = m.cos(Rx)
    Rxs = m.sin(Rx)
    Ryc = m.cos(Ry)
    Rys = m.sin(Ry)
    Rzc = m.cos(Rz)
    Rzs = m.sin(Rz)
    Rxce = m.cos(Rx + epsilon)
    Rxse = m.sin(Rx + epsilon)
    Ryce = m.cos(Ry + epsilon)
    Ryse = m.sin(Ry + epsilon)
    Rzce = m.cos(Rz + epsilon)
    Rzse = m.sin(Rz + epsilon)
    Rxc_e = m.cos(Rx - epsilon)
    Rxs_e = m.sin(Rx - epsilon)
    Ryc_e = m.cos(Ry - epsilon)
    Rys_e = m.sin(Ry - epsilon)
    Rzc_e = m.cos(Rz - epsilon)
    Rzs_e = m.sin(Rz - epsilon)
    dX_Sd =  (Sd_6Dof(dXe,  dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX_e, dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dY_Sd =  (Sd_6Dof(dX,   dYe,  dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY_e, dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dZ_Sd =  (Sd_6Dof(dX,   dY,   dZe,  Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY,   dZ_e, Rxc,   Rxs,   Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRx_Sd = (Sd_6Dof(dX,   dY,   dZ,   Rxce,  Rxse,  Ryc,   Rys,   Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY,   dZ,   Rxc_e, Rxs_e, Ryc,   Rys,   Rzc,   Rzs))   /(2*epsilon)
    dRy_Sd = (Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryce,  Ryse,  Rzc,   Rzs)  -
              Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc_e, Rys_e, Rzc,   Rzs))   /(2*epsilon)
    dRz_Sd = (Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzce,  Rzse) -
              Sd_6Dof(dX,   dY,   dZ,   Rxc,   Rxs,   Ryc,   Rys,   Rzc_e, Rzs_e)) /(2*epsilon)
    
    return dX_Sd, dY_Sd, dZ_Sd, dRx_Sd, dRy_Sd, dRz_Sd, Instrument

def horizontal_angle_from_Coords(PointTo,PointFrom):
    Hz = a(m.atan2(PointTo[1]-PointFrom[1],PointTo[0]-PointFrom[0]), 
           a.T_RAD,True).angle
    return Hz

def vertical_angle_from_Coords(PointTo,PointFrom):
    dist = slope_distance(PointFrom, PointTo)
    V = m.acos((PointTo[2]-PointFrom[2])/dist)
    return V

def find_unknowns(Dict_of_measurements, Instruments_6DoF):
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
    if Instruments_6DoF:
        number_of_unknowns = len(unknown_points)*3 + count_instruments*6
    else:
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
        key = tuple(Dictionary.keys())
        result = Dictionary[key[0]]
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

def filling_X(Aproximates, unknowns, count_unknowns, count_instruments, 
              Instruments_6DoF):
    X = np.zeros(count_unknowns)
    XHR = []
    for i, unknown in enumerate(unknowns[:-count_instruments]):
        iii = 3*i
        X[iii:iii+3] = np.array(Aproximates[unknown])
        # Human readible version of X
        XHR.append(('X ' + unknown))
        XHR.append(('Y ' + unknown))
        XHR.append(('Z ' + unknown))
    if Instruments_6DoF:
        for i, unknown in enumerate(unknowns[-count_instruments:]):
            iii = 3*(count_instruments - i)
            Angles = Aproximates[unknown]
            Rx = (a(Angles[0],a.T_RAD, True).angle)
            Ry = (a(Angles[1],a.T_RAD, True).angle)
            Rz = (a(Angles[2],a.T_RAD, True).angle)
            X[-iii] = Rx
            X[-iii+1] = Ry
            X[-iii+2] = Rz
            XHR.append(('RX ' + unknown))
            XHR.append(('RY ' + unknown))
            XHR.append(('RZ ' + unknown))
    else:
        for i, unknown in enumerate(unknowns[-count_instruments:]):
            index = count_instruments - i
            X[-index] = np.array(a(Aproximates[unknown],a.T_RAD, True).angle)
            XHR.append(unknown)
    return X, XHR

def filling_G(number_of_unknowns, unknowns, Aproximates, 
						    number_of_instruments, Instruments_6DoF):
    if Instruments_6DoF:
        G_matrix = np.zeros([number_of_unknowns,6])
        sumx=0
        sumy=0
        sumz=0
        c=0
        for i,unknown in enumerate(unknowns):
            if i < (len(unknowns) - 2*number_of_instruments):
                sumx += Aproximates[unknown][0]
                sumy += Aproximates[unknown][1]
                sumz += Aproximates[unknown][2]
                c+=1
        meanx = sumx/c
        meany = sumy/c
        meanz = sumz/c
        del c,sumx,sumy,sumz
        
        for i,unknown in enumerate(unknowns):
            if i < (len(unknowns) - 2*number_of_instruments):
                iii = 3*i
                G_matrix[iii,2]   = 1
                G_matrix[iii+1,3] = 1
                G_matrix[iii+2,4] = 1
                G_matrix[iii,1] =   - (Aproximates[unknown][2] - meanz)/100
                G_matrix[iii,5] =     (Aproximates[unknown][1] - meany)/100
                G_matrix[iii+1,0] =   (Aproximates[unknown][2] - meanz)/100
                G_matrix[iii+1,5] = - (Aproximates[unknown][0] - meanx)/100
                G_matrix[iii+2,0] = - (Aproximates[unknown][1] - meany)/100
                G_matrix[iii+2,1] =   (Aproximates[unknown][0] - meanx)/100
        del i, iii, unknown, meany, meanx, meanz
    else:
        G_matrix = np.zeros([number_of_unknowns,4])
        
        sumx=0
        sumy=0
        sumz=0
        c=0
        for i,unknown in enumerate(unknowns):
            if i < (len(unknowns) - 2*number_of_instruments):
                sumx += Aproximates[unknown][0]
                sumy += Aproximates[unknown][1]
                sumz += Aproximates[unknown][2]
                c+=1
        meanx = sumx/c
        meany = sumy/c
        meanz = sumz/c
        del c,sumx,sumy,sumz
        
        for i,unknown in enumerate(unknowns):
            if i < (len(unknowns) - 2*number_of_instruments):
                iii = 3*i
                G_matrix[iii,0]   = 1
                G_matrix[iii+1,1] = 1
                G_matrix[iii+2,2] = 1
                G_matrix[iii,3] =     (Aproximates[unknown][1] - meany)/100
                G_matrix[iii+1,3] = - (Aproximates[unknown][0] - meanx)/100
        del i, iii, unknown, meany, meanx, meanz
    return G_matrix


def filling_Aproximates(unknowns, X_vector, instruments, Instruments_6DoF):
    updated_Aproximates = dict()
    inst_count = len(instruments)
    if Instruments_6DoF:
        for i, unknown in enumerate(unknowns):
            iii = 3*i
            updated_Aproximates[unknown] = tuple(X_vector[iii:iii+3])
    else:
        for i, point in enumerate(unknowns[:-inst_count]):
            iii = 3*i
            updated_Aproximates[point] = tuple(X_vector[iii:iii+3])
        for i,instrument in enumerate(instruments):
            updated_Aproximates['Ori_' + instrument] = X_vector[-inst_count+i]
    return updated_Aproximates

def Filling_A_L_P_LX0(Nominal_coords,Aproximates, Trans_par,
                      Combinations_for_constraints,
                      measured_distances_in_lines,
                      sorted_measured_points_in_lines,
                      instruments, count_instruments,
                      Pol_measurements,unknowns,count_unknowns,
                      X_vector, X_vectorHR, IFM_StDev, Instruments_6DoF, 
                      epsilon
                      ):
    count_Sd = Count_meas_types(Pol_measurements, 'Sd')
    count_Hz = Count_meas_types(Pol_measurements, 'Hz')
    count_V = Count_meas_types(Pol_measurements, 'V')
    
    count_all_observations = count_Sd + count_Hz + count_V
    count_A_rows = count_all_observations
	
    if cg.LSM_incl_IFM:
        count_IFM_measurements = sum([len(v) for k, v in\
                                         measured_distances_in_lines.items()])
        count_A_rows += count_IFM_measurements
		
    if cg.LSM_incl_Cons:
        count_constraints = len(Combinations_for_constraints)
        count_A_rows += count_constraints
	
    A_matrix = np.zeros([count_A_rows, count_unknowns], dtype=float)
    
    A_matrixHR = {}
    
    L_vector = np.array([], dtype=float)
    LX0_vector = np.array([], dtype=float)
    
    L_vectorHR = []
    
    P_vector = np.array([],dtype=float)
    Q_vector = np.array([],dtype=float)
	
    if cg.LSM_incl_IFM:
	    # Filling A and L with IFM measurements
	    for line in measured_distances_in_lines:
	        # meas_i is index of the measurand in line, so I can pair it with the
	        # Point IDs from sorted points in lines, first measured point is on the
	        # same index as the distance and second point is on index + 1
	        for meas_i,distance in enumerate(measured_distances_in_lines[line]):
	            # I need index of the measurement in L vector, so I can put it on
	            # correct corresponding row in First plan matrix A
	            L_i = len(L_vector)
	            L_vector = np.append(L_vector, distance)
	            StDev = (StDev_sys_ppm(distance,IFM_StDev))/1000
	            P_vector = np.append(P_vector,pow(cg.Sigma_0,2) / pow(StDev,2))
	            Q_vector = np.append(Q_vector,1/pow(cg.Sigma_0,2) * pow(StDev,2))
	            """Now figuring the index of points in the unknowns list, so I know
	               which column of the A matrix. The original index is multiplied 
	               by 3 because unknowns are names of points (and instrument 
	               orientations) 
	               not a complete list of unknowns = Point X,Y and Z => *3     """
	            PointFrom = sorted_measured_points_in_lines[line][meas_i]
	            PointTo = sorted_measured_points_in_lines[line][meas_i+1]
	            LX0_vector = np.append(LX0_vector, slope_distance(
	                                  Aproximates[PointTo],Aproximates[PointFrom]))
	            PointFrom_i = 3*unknowns.index(PointFrom) 
	            PointTo_i = 3*unknowns.index(PointTo)
	            if Instruments_6DoF:
	                dX_Sd, dY_Sd, dZ_Sd, dRx_Sd, dRy_Sd, dRz_Sd, inst = \
	                Par_6Dof_IFM(PointTo, PointFrom, line, Aproximates, epsilon)
	                Ori_inst_i = X_vectorHR.index('RX Ori_'+inst)
	                A_matrix[L_i,PointTo_i:PointTo_i+3] = dX_Sd, dY_Sd, dZ_Sd
	                # for point "From" the partial derivatives change sign
	                A_matrix[L_i,PointFrom_i:PointFrom_i+3] = -dX_Sd, -dY_Sd, -dZ_Sd
	                A_matrix[L_i,Ori_inst_i:Ori_inst_i+3] = -dRx_Sd, -dRy_Sd, -dRz_Sd

	            else:
	                dX,dY,dZ = ParD_Sd(Aproximates[PointTo],Aproximates[PointFrom])
	                A_matrix[L_i,PointTo_i:PointTo_i+3] = dX,dY,dZ
	                # for point "From" the partial derivatives change sign
	                A_matrix[L_i,PointFrom_i:PointFrom_i+3] = -dX,-dY,-dZ
	            # Documenting in human readible format the A and L elements
	            L_vectorHR.append((line, PointFrom, PointTo))
	            A_matrixHR[(L_i,PointFrom_i)] = ['dX', line, PointFrom, PointTo]
	            A_matrixHR[(L_i,PointFrom_i+1)] = ['dY', line, PointFrom, PointTo]
	            A_matrixHR[(L_i,PointFrom_i+2)] = ['dZ', line, PointFrom, PointTo]
	            A_matrixHR[(L_i,PointTo_i)] = ['dX', line, PointTo, PointFrom]
	            A_matrixHR[(L_i,PointTo_i+1)] = ['dY', line, PointTo, PointFrom]
	            A_matrixHR[(L_i,PointTo_i+2)] = ['dZ', line, PointTo, PointFrom]
	    del PointFrom,PointTo, PointFrom_i,PointTo_i,L_i, meas_i, distance
	    try:
	        del dX, dY, dZ
	    except UnboundLocalError:
	        pass
	    try:
	        del dX_Sd, dY_Sd, dZ_Sd
	    except UnboundLocalError:
	        pass
 
    L_subv_Hz = np.array([])
    L_subv_V = np.array([])
    L_subv_Sd = np.array([])
    L_subv_Hz_HR = []
    L_subv_V_HR = []
    L_subv_Sd_HR = []
    
    P_subv_Hz = np.array([])
    P_subv_V = np.array([])
    P_subv_Sd = np.array([])

    Q_subv_Hz = np.array([])
    Q_subv_V = np.array([])
    Q_subv_Sd = np.array([])
    
    LX0_subv_Hz = np.array([])
    LX0_subv_V = np.array([])
    LX0_subv_Sd = np.array([])

    # Estimate orientation unknowns
#    for instrument in instruments:
#        ori_sum = 0
#        count = 0
#        for target in Pol_measurements[instrument]:
#            ori_local = a(Pol_measurements[instrument][target][1],a.T_GON) - \
#                        a(horizontal_angle_from_Coords(Aproximates[target],
#														  Aproximates[instrument]),a.T_RAD)
#            ori_sum += ori_local
#            count += 1
#        Aproximates['Ori_' + instrument] = a(ori_sum/count,a.T_GON,True).angle
#        X_vector[X_vectorHR.index('Ori_' + instrument)] = a(-ori_sum/\
#				count,a.T_GON,True).angle

    counter = 0
    for inst,points in Pol_measurements.items():
        instrument_i = 3 * unknowns.index(inst) # Starting column of the instrument
        if Instruments_6DoF:
            Ori_inst_i = X_vectorHR.index('RZ Ori_'+inst)
        else:
            Ori_inst_i = X_vectorHR.index('Ori_'+inst) # Inst's Orientation index
        for point in points:
            # Row offsets of V and Sd measurements (allows filling in same loop)
            Hz_offset = count_all_observations - (count_Hz + count_Sd + count_V)
            V_offset = count_all_observations - (count_Hz + count_Sd)
            Sd_offset = count_all_observations - count_Sd
            if cg.LSM_incl_IFM:
                Hz_offset += count_IFM_measurements
                V_offset += count_IFM_measurements
                Sd_offset += count_IFM_measurements
            # Returning the starting column of the point
            Point_i = 3*unknowns.index(point)
                        
            if Instruments_6DoF:
                dX_Sd, dY_Sd, dZ_Sd, dRx_Sd, dRy_Sd, dRz_Sd, \
                dX_Hz, dY_Hz, dZ_Hz, dRx_Hz, dRy_Hz, dRz_Hz, \
                dX_V, dY_V, dZ_V, dRx_V, dRy_V, dRz_V = Par_6Dof(point, inst, 
                                                                 Aproximates, 
                                                                 epsilon)
            else:
                 Hz_dX, Hz_dY, Hz_dZ, Hz_dO = ParD_Hz(Aproximates[point],
                                                    Aproximates[inst])
                 V_dX, V_dY, V_dZ = ParD_V(Aproximates[point],
                                           Aproximates[inst])
                 Sd_dX, Sd_dY, Sd_dZ = ParD_Sd(Aproximates[point],
                                               Aproximates[inst])   
            
            # Filling the L subvectors for Hz,V and Sd
            
            # handeling of excluded measurements using KeyError for each type:
            try:
                Pol_measurements[inst][point]['Sd']
            except KeyError:
                pass#print(inst, point + " Sd")
            else:
                Sd = Pol_measurements[inst][point]['Sd']
                StDev_Sd = Pol_measurements[inst][point]['StDev_Sd']
                L_subv_Sd = np.append(L_subv_Sd, Sd)
                P_subv_Sd = np.append(P_subv_Sd, pow(cg.Sigma_0,2) / pow(
																StDev_Sd,2))
                Q_subv_Sd = np.append(Q_subv_Sd, 1/pow(cg.Sigma_0,2) * pow(
																StDev_Sd,2))
                if Instruments_6DoF:
                    LX0_subv_Sd = np.append(LX0_subv_Sd, slope_distance_6DoF(
                                         Aproximates,point,inst))
                else:
                    LX0_subv_Sd = np.append(LX0_subv_Sd, slope_distance(
                                         Aproximates[point],Aproximates[inst]))
                L_subv_Sd_HR.append(('Sd', inst, point))
                if Instruments_6DoF:
                    A_matrix[Sd_offset+counter,Point_i] = dX_Sd
                    A_matrix[Sd_offset+counter,Point_i+1] = dY_Sd
                    A_matrix[Sd_offset+counter,Point_i+2] = dZ_Sd
                    A_matrix[Sd_offset+counter,instrument_i] = -dX_Sd
                    A_matrix[Sd_offset+counter,instrument_i+1] = -dY_Sd
                    A_matrix[Sd_offset+counter,instrument_i+2] = -dZ_Sd
#                    A_matrix[Sd_offset+counter,Ori_inst_i-2] = dRx_Sd
#                    A_matrix[Sd_offset+counter,Ori_inst_i-1] = dRy_Sd
#                    A_matrix[Sd_offset+counter,Ori_inst_i] = dRz_Sd
                    A_matrix[Sd_offset+counter,Ori_inst_i-2] = -dRx_Sd   # Markus
                    A_matrix[Sd_offset+counter,Ori_inst_i-1] = -dRy_Sd   # Markus
                    A_matrix[Sd_offset+counter,Ori_inst_i]   = -dRz_Sd   # Markus
                    A_matrixHR[(Sd_offset+counter,Point_i)] = ['Sd/dX', inst, 
                                                                   point]
                    A_matrixHR[(Sd_offset+counter,Point_i+1)] = ['Sd/dY', inst, 
                                                                   point]
                    A_matrixHR[(Sd_offset+counter,Point_i+2)] = ['Sd/dZ', inst, 
                                                                         point]
                    A_matrixHR[(Sd_offset+counter,instrument_i)] = ['Sd/dX',
                                                               point, inst]
                    A_matrixHR[(Sd_offset+counter,instrument_i+1)] = ['Sd/dY',
                                                                    point,inst]
                    A_matrixHR[(Sd_offset+counter,instrument_i+2)] = ['Sd/dZ',
                                                                    point,inst]
                    A_matrixHR[(Sd_offset+counter,Ori_inst_i-2)] = ['Sd/dRx',
                                                               point, inst]
                    A_matrixHR[(Sd_offset+counter,Ori_inst_i-1)] = ['Sd/dRy',
                                                                    point,inst]
                    A_matrixHR[(Sd_offset+counter,Ori_inst_i)] = ['Sd/dRz',
                                                                    point,inst]
                else:
                    A_matrix[Sd_offset+counter,Point_i:Point_i+3] = Sd_dX,  \
                                                                Sd_dY,Sd_dZ
                    A_matrix[Sd_offset+counter,instrument_i:instrument_i+3] = \
                                                        -Sd_dX, -Sd_dY, -Sd_dZ
                
                    A_matrixHR[(Sd_offset+counter,Point_i)] = ['Sd/dX', inst, 
                                                                   point]
                    A_matrixHR[(Sd_offset+counter,Point_i+1)] = ['Sd/dY', inst, 
                                                                   point]
                    A_matrixHR[(Sd_offset+counter,Point_i+2)] = ['Sd/dZ', inst, 
                                                                   point]
                    A_matrixHR[(Sd_offset+counter,instrument_i)] = ['Sd/dX',
                                                               point, inst]
                    A_matrixHR[(Sd_offset+counter,instrument_i+1)] = ['Sd/dY',
                                                                 point,inst]
                    A_matrixHR[(Sd_offset+counter,instrument_i+2)] = ['Sd/dZ',
                                                                   point,inst]
            try:
                Pol_measurements[inst][point]['Hz']
            except KeyError:
                pass#print(inst, point + " Hz")
            else:
                Hz = Pol_measurements[inst][point]['Hz']
                StDev_Hz = Pol_measurements[inst][point]['StDev_Hz']
                L_subv_Hz = np.append(L_subv_Hz, a(Hz,a.T_GON,True).angle)
                P_subv_Hz = np.append(P_subv_Hz, pow(cg.Sigma_0,2) / pow(
                                                    gon2rad(StDev_Hz/1000),2))
                Q_subv_Hz = np.append(Q_subv_Hz, 1/pow(cg.Sigma_0,2) * pow(
                                                    gon2rad(StDev_Hz/1000),2))
                if Instruments_6DoF:
                    Hz_angle_from_aprox = a(horizontal_angle_6DoF(Aproximates,
													   point, inst),a.T_RAD,True).angle
#                    Hz_angle_from_aprox = a(horizontal_angle_from_Coords(
#                                    Aproximates[point],Aproximates[inst]) + \
#                                    X_vector[Ori_inst_i],a.T_RAD,True).angle
                else:
                    Hz_angle_from_aprox = a(horizontal_angle_from_Coords(
                                    Aproximates[point],Aproximates[inst]) - \
                                    X_vector[Ori_inst_i],a.T_RAD,True).angle
                        
                LX0_subv_Hz = np.append(LX0_subv_Hz, Hz_angle_from_aprox)
                L_subv_Hz_HR.append(('Hz', inst, point))
                
                if Instruments_6DoF:
                    A_matrix[Hz_offset+counter,Point_i:Point_i+3] = dX_Hz,\
                                                                dY_Hz, dZ_Hz
                    A_matrix[Hz_offset+counter,instrument_i:instrument_i+3] = \
                                                       -dX_Hz, -dY_Hz, -dZ_Hz
#                    A_matrix[Hz_offset+counter,Ori_inst_i-2:Ori_inst_i+1] = \
#                                                    dRx_Hz, dRy_Hz, dRz_Hz
                    A_matrix[Hz_offset+counter,Ori_inst_i-2:Ori_inst_i+1] = -dRx_Hz, -dRy_Hz, -dRz_Hz  # Markus
                                                    
                    A_matrixHR[(Hz_offset+counter,Point_i)] = ['Hz/dX', inst, 
                                                                   point]
                    A_matrixHR[(Hz_offset+counter,Point_i+1)] = ['Hz/dY', inst, 
                                                                   point]
                    A_matrixHR[(Hz_offset+counter,Point_i+2)] = ['Hz/dZ', inst, 
                                                                   point]
                    A_matrixHR[(Hz_offset+counter,instrument_i)] = ['Hz/dX',
                                                               point,inst]
                    A_matrixHR[(Hz_offset+counter,instrument_i+1)] = ['Hz/dY',
                                                               point,inst]
                    A_matrixHR[(Hz_offset+counter,instrument_i+2)] = ['Hz/dZ',
                                                               point,inst]
                    A_matrixHR[(Hz_offset+counter,Ori_inst_i-2)] = ['Hz/dRx',
                                                               point, inst]
                    A_matrixHR[(Hz_offset+counter,Ori_inst_i-1)] = ['Hz/dRy',
                                                                    point,inst]
                    A_matrixHR[(Hz_offset+counter,Ori_inst_i)] = ['Hz/dRz',
                                                                    point,inst]
                else:
                    A_matrix[Hz_offset+counter,Point_i:Point_i+3] = Hz_dX,\
                                                                Hz_dY, Hz_dZ
                    A_matrix[Hz_offset+counter,Ori_inst_i] = Hz_dO
                    A_matrix[Hz_offset+counter,instrument_i:instrument_i+3] = \
                                                       -Hz_dX, -Hz_dY, -Hz_dZ
                    A_matrixHR[(Hz_offset+counter,Point_i)] = ['Hz/dX', inst, 
                                                                   point]
                    A_matrixHR[(Hz_offset+counter,Point_i+1)] = ['Hz/dY', inst, 
                                                                   point]
                    A_matrixHR[(Hz_offset+counter,Point_i+2)] = ['Hz/dZ', inst, 
                                                                   point]
                    A_matrixHR[(Hz_offset+counter,instrument_i)] = ['Hz/dX',
                                                               point,inst]
                    A_matrixHR[(Hz_offset+counter,instrument_i+1)] = ['Hz/dY',
                                                               point,inst]
                    A_matrixHR[(Hz_offset+counter,instrument_i+2)] = ['Hz/dZ',
                                                               point,inst]
                
            try:
                Pol_measurements[inst][point]['V']
            except KeyError:
                pass#print(inst, point + " V")
            else:
                V = Pol_measurements[inst][point]['V']
                StDev_V = Pol_measurements[inst][point]['StDev_V']
                L_subv_V = np.append(L_subv_V, a(V,a.T_GON,True).angle)
                P_subv_V = np.append(P_subv_V, pow(cg.Sigma_0,2) / pow(
                                                  gon2rad(StDev_V/1000),2))
                Q_subv_V = np.append(Q_subv_V, 1/pow(cg.Sigma_0,2) * pow(
                                                  gon2rad(StDev_V/1000),2))
                if Instruments_6DoF:
                    LX0_subv_V = np.append(LX0_subv_V, vertical_angle_6DoF(
                                           Aproximates,point,inst))
                else:
                    LX0_subv_V = np.append(LX0_subv_V, vertical_angle_from_Coords(
                                        Aproximates[point],Aproximates[inst]))
                L_subv_V_HR.append(('V', inst, point))

                if Instruments_6DoF:
                    A_matrix[V_offset+counter,Point_i:Point_i+3] = dX_V, dY_V,\
                                                                    dZ_V        
                    A_matrix[V_offset+counter,instrument_i:instrument_i+3] = \
                                                          -dX_V, -dY_V, -dZ_V
#                    A_matrix[V_offset+counter,Ori_inst_i-2:Ori_inst_i+1] = \
#                                                          -dRx_V, -dRy_V, -dRz_V 
                    A_matrix[V_offset+counter,Ori_inst_i-2:Ori_inst_i+1] = -dRx_V, -dRy_V, -dRz_V # Markus 
                    A_matrixHR[(V_offset+counter,Point_i)] = ['V/dX', inst, 
                                                                        point]
                    A_matrixHR[(V_offset+counter,Point_i+1)] = ['V/dY', inst, 
                                                                        point]
                    A_matrixHR[(V_offset+counter,Point_i+2)] = ['V/dZ', inst,
                                                                        point]
                    A_matrixHR[(V_offset+counter,instrument_i)] = ['V/dX',
                                                                   point,inst]
                    A_matrixHR[(V_offset+counter,instrument_i+1)] = ['V/dY',
                                                                   point,inst]
                    A_matrixHR[(V_offset+counter,instrument_i+2)] = ['V/dZ',
                                                                   point,inst]
                    A_matrixHR[(V_offset+counter,Ori_inst_i-2)] = ['V/dRx',
                                                               point, inst]
                    A_matrixHR[(V_offset+counter,Ori_inst_i-1)] = ['V/dRy',
                                                                    point,inst]
                    A_matrixHR[(V_offset+counter,Ori_inst_i)] = ['V/dRz',
                                                                    point,inst]
                else:
                    A_matrix[V_offset+counter,Point_i:Point_i+3] = V_dX, V_dY,\
                                                                   V_dZ         
                    A_matrix[V_offset+counter,instrument_i:instrument_i+3] = \
                                                          -V_dX, -V_dY, -V_dZ
                    A_matrixHR[(V_offset+counter,Point_i)] = ['V/dX', inst, 
                                                                        point]
                    A_matrixHR[(V_offset+counter,Point_i+1)] = ['V/dY', inst, 
                                                                        point]
                    A_matrixHR[(V_offset+counter,Point_i+2)] = ['V/dZ', inst,
                                                                        point]
                    A_matrixHR[(V_offset+counter,instrument_i)] = ['V/dX',
                                                                   point,inst]
                    A_matrixHR[(V_offset+counter,instrument_i+1)] = ['V/dY',
                                                                   point,inst]
                    A_matrixHR[(V_offset+counter,instrument_i+2)] = ['V/dZ',
                                                                   point,inst]
            
            counter += 1
    L_vector = np.concatenate((L_vector,L_subv_Hz))
    L_vector = np.concatenate((L_vector,L_subv_V))
    L_vector = np.concatenate((L_vector,L_subv_Sd))
    LX0_vector = np.concatenate((LX0_vector,LX0_subv_Hz))
    LX0_vector = np.concatenate((LX0_vector,LX0_subv_V))
    LX0_vector = np.concatenate((LX0_vector,LX0_subv_Sd))
    L_vectorHR = L_vectorHR + L_subv_Hz_HR + L_subv_V_HR + L_subv_Sd_HR
    P_vector = np.concatenate((P_vector,P_subv_Hz))
    P_vector = np.concatenate((P_vector,P_subv_V))
    P_vector = np.concatenate((P_vector,P_subv_Sd))
    Q_vector = np.concatenate((Q_vector,Q_subv_Hz))
    Q_vector = np.concatenate((Q_vector,Q_subv_V))
    Q_vector = np.concatenate((Q_vector,Q_subv_Sd))


# =============================================================================    
# Constraints - filling A, L, LX0, P and HR versions
# =============================================================================
    if cg.LSM_incl_Cons:
	    L_lenght = len(L_vector)
	    for index, const in enumerate(Combinations_for_constraints):
	        A_row_index = L_lenght + index
	        PointFrom = Combinations_for_constraints[index][0]
	        PointTo = Combinations_for_constraints[index][1]
	        Sd = slope_distance(Nominal_coords[PointFrom],Nominal_coords[PointTo])
	        Sd_aprox = slope_distance(Aproximates[PointFrom],Aproximates[PointTo])
	        PointFrom_i = 3*unknowns.index(PointFrom) 
	        PointTo_i = 3*unknowns.index(PointTo)
	        if Instruments_6DoF:
	            dX_Sd, dY_Sd, dZ_Sd = Par_6DoF_noRot(PointTo, PointFrom, Aproximates,
																				     epsilon)
	            A_matrix[A_row_index,PointTo_i:PointTo_i+3] = dX_Sd, dY_Sd, dZ_Sd
	            # for point "From" the partial derivatives change sign
	            A_matrix[A_row_index,PointFrom_i:PointFrom_i+3] = -dX_Sd, -dY_Sd, \
	                                                              -dZ_Sd
	        else:
	            dX,dY,dZ = ParD_Sd(Nominal_coords[PointTo],Nominal_coords[PointFrom])
	            A_matrix[A_row_index,PointTo_i:PointTo_i+3] = dX,dY,dZ
	            # for point "From" the partial derivatives change sign
	            A_matrix[A_row_index,PointFrom_i:PointFrom_i+3] = -dX,-dY,-dZ
	        # Documenting in human readible format what are the A and L elements
	        L_vector = np.append(L_vector,Sd)
	        LX0_vector = np.append(LX0_vector,Sd_aprox)
	        P_vector = np.append(P_vector,pow(cg.Sigma_0,2)/pow(cg.Constraint_StDev,2))
	        L_vectorHR.append(('constraint', PointFrom, PointTo))
	        A_matrixHR[(A_row_index,PointFrom_i)] = ['dX', 'distance constraint', 
	                   PointFrom, PointTo]
	        A_matrixHR[(A_row_index,PointFrom_i+1)] = ['dY', 'distance constraint', 
	                   PointFrom, PointTo]
	        A_matrixHR[(A_row_index,PointFrom_i+2)] = ['dZ', 'distance constraint', 
	                   PointFrom, PointTo]
	        A_matrixHR[(A_row_index,PointTo_i)] = ['dX', 'distance constraint', 
	                   PointTo, PointFrom]
	        A_matrixHR[(A_row_index,PointTo_i+1)] = ['dY', 'distance constraint', 
	                   PointTo, PointFrom]
	        A_matrixHR[(A_row_index,PointTo_i+2)] = ['dZ', 'distance constraint', 
	                   PointTo, PointFrom]
	    try:
	        del dX, dY, dZ
	    except UnboundLocalError:
	        pass
	    try:
	        del dX_Sd, dY_Sd, dZ_Sd
	    except UnboundLocalError:
	        pass
    
    P_matrix = np.diagflat(P_vector)
    Q_matrix = np.diagflat(Q_vector)
        
    LSM_can_be_done = (len(P_vector) ==  len(L_vector)) and (
                    len(L_vector) == A_matrix.shape[0]) and (
                                         A_matrix.shape[0] == len(LX0_vector))

    if not LSM_can_be_done:
        sys.exit("Error during filling P, L, A, and LX0, sizes are not same")
  
    if len(X_vector) != A_matrix.shape[1]:
        sys.exit("Error during filling A, X, sizes don't match.")
        LSM_can_be_done = False
    return LSM_can_be_done,A_matrix,L_vector,P_matrix,Q_matrix,\
           LX0_vector,A_matrixHR,L_vectorHR

def create_constraints(Aproximates):
    Combinations_for_constraints = []
    for magnet in cg.Names_of_magnets:
        points_on_magnet = [key for key in Aproximates.keys() if magnet in key]
        magnet_combinations = [(PointA, PointB) for index, PointA in enumerate(
            points_on_magnet) for PointB in points_on_magnet[index + 1:]]
        Combinations_for_constraints.extend(magnet_combinations)
    count_constraints = len(Combinations_for_constraints)
    return Combinations_for_constraints, count_constraints


def LSM(Epoch_num, Nominal_coords, Aproximates, measured_distances_in_lines,
				 sorted_measured_points_in_lines,instruments, count_instruments,
				 Pol_measurements,count_Pol_measurements, count_IFM, unknowns, 
          count_unknowns, IFM_StDev, Instruments_6DoF, Trans_par, epsilon):
	
    Combinations_for_constraints,count_constraints =\
														create_constraints(Aproximates)

    G_matrix = filling_G(count_unknowns, unknowns, 
									Aproximates,count_instruments, Instruments_6DoF)

    X_vector, X_vectorHR = filling_X(Aproximates, unknowns, count_unknowns, 
                                     count_instruments, Instruments_6DoF)

    LSM_can_be_done,A_matrix,L_vector,P_matrix,Q_matrix,LX0_vector,A_matrixHR,\
    L_vectorHR  = Filling_A_L_P_LX0(Nominal_coords,Aproximates,Trans_par,
                          Combinations_for_constraints,
                          measured_distances_in_lines,
                          sorted_measured_points_in_lines,
                          instruments, count_instruments,
                          Pol_measurements,
                          unknowns,count_unknowns,
                          X_vector, X_vectorHR, IFM_StDev,Instruments_6DoF,
                          epsilon
                          )

    
    metric = cg.LSM_Threshold + 1
    counter = 0
    print('\n Processing Epoch_%s:' %(str(Epoch_num)))
    while (metric > cg.LSM_Threshold) and (counter <= cg.LSM_Max_iterations):
#        print('\n Iteration', counter)
        l = np.ndarray(L_vector.size)
        for i, lelement in enumerate(L_vector):
            if L_vectorHR[i][0] == "Hz":
                l[i] = (a(LX0_vector[i] - L_vector[i],a.T_RAD,True).angle)
            else:
                l[i] = LX0_vector[i] - L_vector[i]
        del lelement, i
#        print('l',max(abs(l)), np.argmax(abs(l)))
        N = A_matrix.transpose() @ P_matrix @ A_matrix
#        print("Determinant: ",linalg.det(A_matrix.transpose() @ A_matrix))
        if Instruments_6DoF:
            O = np.zeros([6,6])
        else:
            O = np.zeros([4,4])
        N_extended = np.block([[N,G_matrix],[np.transpose(G_matrix),O]])
#        print('Rank N_extended:', np.linalg.matrix_rank(N_extended))
#        try:
#            print("Log-Determinant of G-extended N: ",
#    			   np.linalg.slogdet(N_extended))
#        except:
#            print("Log-Determinant of G-extended N gives overflow ")
#        print("Condition number N: ", np.linalg.cond(N_extended))
        n = A_matrix.transpose() @ P_matrix @ l
        if Instruments_6DoF:
            N_inv = np.linalg.inv(N_extended)[:-6,:-6]
        else:
            N_inv = np.linalg.inv(N_extended)[:-4,:-4]
        dx = N_inv @ n
#        print('dx',max(abs(dx)), np.argmax(abs(dx)))
#        print("X")
#        pretty_print((X_vector*pow(10,9))[-6:])
#        print("dx")
#        pretty_print((dx*pow(10,9))[-6:])
#        pretty_print((X_vector*pow(10,9))[-7:])
        dx_1 = -dx[:-6]
        dx_2 = dx[-6:]
        dx = np.concatenate((dx_1,dx_2))
        X_vector += dx
        v = A_matrix @ dx - l
#        print('v max ', v[np.argmax(abs(v))], np.argmax(abs(v)))
#        print('[vv]  ', sum(v*v))
#        print('[v]  ', sum(v))

        P_inv = inv(P_matrix)
#        Qllest = A_matrix @ N_inv @ A_matrix.transpose()
#        Qvv = P_inv - Qllest
#        Qvv = Q_matrix - Qllest
        Qvv = P_inv - A_matrix @ N_inv @ A_matrix.transpose()
        dof = int(round(np.trace(Qvv@P_matrix)))
#        print("dof from Qvv: ",dof)
        s02 = (v @ P_matrix @ v)/dof
        LSM_results = Aproximates.copy()
        Qxx = N_inv
        Cov_matrix = s02 * Qxx
        w = []
        vcount=0
        for i,vi in enumerate(v):
            if Qvv[i,i] < 0:
                vcount += 1
#                print(f"{i}, vi {vi:8.2}, Qll {Q_matrix[i,i]:7.2}, Ql^l^ {Qllest[i][i]:8.2}, {Qvv[i,i]:8.2}, {L_vectorHR[i]}")
#                pass
            else:
#                if (vi/(m.sqrt(s02*Qvv[i,i]))>4):
#                    print(f"{vi/(m.sqrt(s02*Qvv[i,i])):16.12}, {L_vectorHR[i]}")
                    pass
#        print("n/o diag elem <0: ",vcount)
        if cg.LSM_incl_IFM:
            i = count_IFM
            s02_IFM = pow(cg.Sigma_0,-2)*(v[0:i-1] @ P_matrix[0:i-1,0:i-1] @ v[0:i-1])/\
										round(np.trace(Qvv[0:i-1,0:i-1] @ P_matrix[0:i-1,0:i-1]))
        else:
            i = 0
            s02_IFM = np.nan

        p = count_Pol_measurements

   
        s02_Hz = pow(cg.Sigma_0,-2)*(v[i:i+p-1] @ P_matrix[i:i+p-1,i:i+p-1] @ \
						   v[i:i+p-1])/ round(np.trace(Qvv[i:i+p-1,i:i+p-1] @ \
						   P_matrix[i:i+p-1,i:i+p-1]))			
        s02_V = pow(cg.Sigma_0,-2)*(v[i+p:i+2*p-1] @ \
				           P_matrix[i+p:i+2*p-1,i+p:i+2*p-1] @ v[i+p:i+2*p-1])/\
				           round(np.trace(Qvv[i+p:i+2*p-1,i+p:i+2*p-1] @ \
							    P_matrix[i+p:i+2*p-1,i+p:i+2*p-1]))
        s02_Sd = pow(cg.Sigma_0,-2)*(v[i+2*p:i+3*p-1] @ \
					         P_matrix[i+2*p:i+3*p-1,i+2*p:i+3*p-1] @ v[i+2*p:i+3*p-1])\
									/ round(np.trace(Qvv[i+2*p:i+3*p-1,i+2*p:i+3*p-1]\
											   @ P_matrix[i+2*p:i+3*p-1,i+2*p:i+3*p-1]))

        if cg.LSM_incl_Cons:
            cc = count_constraints
            s02_con = pow(cg.Sigma_0,-2)*(v[-cc:] @ P_matrix[-cc:,-cc:] @ v[-cc:])/\
										round(np.trace(Qvv[-cc:,-cc:] @ P_matrix[-cc:,-cc:]))
        else:
            s02_con = np.nan

        metric = max(abs(dx))
        counter += 1

        Aproximates = filling_Aproximates(unknowns, X_vector, instruments, 
                                          Instruments_6DoF)
        LSM_can_be_done,A_matrix,L_vector,P_matrix,Q_matrix,LX0_vector,A_matrixHR, \
        L_vectorHR = Filling_A_L_P_LX0(Nominal_coords,Aproximates,Trans_par,
                                      Combinations_for_constraints,
                                      measured_distances_in_lines,
                                      sorted_measured_points_in_lines,
                                      instruments, count_instruments,
                                      Pol_measurements,unknowns,count_unknowns,
                                      X_vector, X_vectorHR, IFM_StDev, 
                                      Instruments_6DoF, epsilon
                                      )

        

    print (f"s02    {s02:8.3}")
    print (f"s0     {m.sqrt(s02):8.3}")
    print (f"sig0   {cg.Sigma_0:8.3}")
    if cg.LSM_incl_IFM:
        print (f"s02ifm {s02_IFM:8.3}")
    print (f"s02Hz  {s02_Hz:8.3}")
    print (f"s02V   {s02_V:8.3}")
    print (f"s02Sd  {s02_Sd:8.3}")
    try:
        del cc
    except UnboundLocalError:
        pass
    del metric, counter, i, p
    return P_matrix, LSM_results, Qxx, Qvv, Cov_matrix, s02, dof, w, s02_IFM, s02_Hz, s02_V, s02_Sd, s02_con