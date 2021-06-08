# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:23:20 2021

@author: jbarker
"""

""" SET BACKEND A magical trick from 
https://stackoverflow.com/questions/47356726/
fix-matplotlib-not-installed-as-framework-error-w-out-changing-matplotlib-con"""

#import matplotlib as mpl
#mpl.use('TkAgg')
from collections import namedtuple
import numpy as np
#import matplotlib.pyplot as plt
#import scipy as sp
import math
#import sys
import functions as fc

def X_Rotation(alpha):
    Rxc = math.cos(alpha)
    Rxs = math.sin(alpha)
    Rx = np.array([[1, 0, 0], [0, Rxc, Rxs], [0, -Rxs, Rxc]])
    return Rx

def Y_Rotation(beta):
    Ryc = math.cos(beta)
    Rys = math.sin(beta)
    Ry = np.array([[Ryc, 0, -Rys], [0, 1, 0], [Rys, 0, Ryc]])
    return Ry

def Z_Rotation(gamma):
    Rzc = math.cos(gamma)
    Rzs = math.sin(gamma)
    Rz = np.array([[Rzc, Rzs, 0], [-Rzs, Rzc, 0], [0, 0, 1]])
    return Rz

def dX_Rotation(alpha):
    Rxc = math.cos(alpha)
    Rxs = math.sin(alpha)
    Rdx = np.array([[0, 0, 0], [0, -Rxs, Rxc], [0, -Rxc, -Rxs]])
    return Rdx

def dY_Rotation(beta):
    Ryc = math.cos(beta)
    Rys = math.sin(beta)
    Rdy = np.array([[-Rys, 0, -Ryc], [0, 0, 0], [Ryc, 0, -Rys]])
    return Rdy

def dZ_Rotation(gamma):
    Rzc = math.cos(gamma)
    Rzs = math.sin(gamma)
    Rdz = np.array([[-Rzs, Rzc, 0], [-Rzc, -Rzs, 0], [0, 0, 0]])
    return Rdz
 
def Rotation_matrix(angle_tuple):
    R = np.array(X_Rotation(angle_tuple[0])*Y_Rotation(
                                    angle_tuple[1])*Z_Rotation(angle_tuple[2]))
    return R

def Rotation_matrix(alpha, beta, gamma):
    R = np.array(X_Rotation(alpha)*Y_Rotation(beta)*Z_Rotation(gamma))
    return R

def Helmert_transformation(x, From):
    """3D Helmert transformation with known transformation Key 
    From is a dictionary of points
    (Rotation matrix parameters, Translation vector and scale in tuple)"""
    T = np.array(x[0:3])
    q = float(x[3])
    R = Rotation_matrix(x[-3:])
    From_transformed = {}
    for point in From:
        From_transformed[point] = T + q * R.dot(np.array(From[point]))
    return From_transformed

def Helmert_aproximate_parameters(From,To):
    identicals = list(set(To.keys()) & set(From.keys()))
    if len(identicals) > 3:
        #MAKE BETTER CHOICE ON POINTS WHEN YOU HAVE TIME, IF YOU WANT... PLEASE
        point1_To = np.array(To[identicals[0]])
        point2_To = np.array(To[identicals[-1]])
        point3_To = np.array(To[identicals[len(identicals)//2]])
        point1_To_original = point1_To.copy().transpose()
        #// truncating division = rounds result to lower int
        point1_From = np.array(From[identicals[0]])
        point2_From = np.array(From[identicals[-1]])
        point3_From = np.array(From[identicals[len(identicals)//2]])
        point1_From_original = point1_From.copy().transpose()
        # Translation of the points to have origin in point 1:
        point2_To = point2_To - point1_To
        point3_To = point3_To - point1_To
        point1_To = (0,0,0)
        point2_From = point2_From - point1_From
        point3_From = point3_From - point1_From
        point1_From = (0,0,0)
        # Calculating the aproximate parameters 
        # (based on Angle between planes paper)
        # First rotation angle calculations
        psi_To = math.atan2(point2_To[1],point2_To[0])
        psi_From = math.atan2(point2_From[1],point2_From[0])
        # Applying the calculated angle to point 2 and 3
        FirstZrotation2_To = (Z_Rotation(psi_To)).dot(point2_To).transpose()
        FirstZrotation3_To = (Z_Rotation(psi_To)).dot(point3_To.transpose())
        FirstZrotation2_From = (Z_Rotation(psi_From)).dot(
                                                         point2_To.transpose())
        FirstZrotation3_From = (Z_Rotation(psi_From)).dot(
                                                         point3_To.transpose())
        # Calculating second rotation angle from the first rotated coordinates
        fi_To = math.atan2(FirstZrotation2_To[1],FirstZrotation2_To[2])
        fi_From = math.atan2(FirstZrotation2_From[1],FirstZrotation2_From[2])
        # Applying the calculated angle to point 3
        SecondYrotation3_To = (Y_Rotation(fi_To)).dot(
                                                FirstZrotation3_To.transpose())
        SecondYrotation3_From = (Y_Rotation(fi_From)).dot(
                                              FirstZrotation3_From.transpose())
        # Calculating third rotation angle
        theta_To = math.atan2(SecondYrotation3_To[1],SecondYrotation3_To[0])
        theta_From = math.atan2(SecondYrotation3_From[1],
                                SecondYrotation3_From[0])
        # Using all three angles, the full rotation matrix for From To is made
        R_To = Z_Rotation(theta_To).dot(Y_Rotation(fi_To)).dot(Z_Rotation(
                                                                       psi_To))
        R_From = Z_Rotation(theta_From).dot(Y_Rotation(fi_From)).dot(
                                    Z_Rotation(psi_From))
        # The translation vector is calculated by rotating the original 
        # point1_From to the "To" coordinate frame. By substracting the rotated 
        # point1_From from point1_To we get the translation vector.
        R0 = R_To.transpose().dot(R_From)
        Translation = tuple(point1_To_original - R0.dot(point1_From_original))
        # Euler rotation angles
        alpha = math.atan2(R0[1,2],R0[2,2])
        beta = - math.atan(R0[0,0])
        gamma = math.atan2(R0[0,1],R0[0,0])
        R_angles = (alpha, beta, gamma)
        x0 = Translation + (1.0,) + R_angles
    else:
        print('Not enough identical points for transformation calculations.')
    return R0, x0# RX + T

To = {'Point1': (3970673.003, 1018563.740, 4870369.178),
      'Point2': (3970667.574, 1018565.195, 4870373.010),
      'Point3': (3970659.461, 1018571.269, 4870377.881),
      'Point4': (3970654.604, 1018577.517, 4870380.020),
      'Point5': (3970650.090, 1018580.577, 4870382.774),
      'Point6': (3970646.096, 1018581.683, 4870385.620)
      }
    
From = {'Point1': (744970.551, 1040944.109, 224.592),
        'Point2': (744966.969, 1040938.331, 224.390),
        'Point3': (744958.051, 1040931.492, 224.057),
        'Point4': (744950.344, 1040928.731, 223.676),
        'Point5': (744945.677, 1040924.795, 223.472),
        'Point6': (744943.006, 1040920.538, 223.352)
      }

""" Testing of aproximate parameters
From_transed = {}
for point in From:
    From_transed[point] = Translation + R.dot(np.array(From[point]).transpose())
    distance12 = fc.slope_distance(From['Point1'],From[point])
    distance121 = fc.slope_distance(tuple(From_transed['Point1']),tuple(From_transed[point]))
    print(distance12, distance121)
print(From_transed)"""

def Build_TFrom(x,From,identicals):
    """x are the transform parameters T,q,R in a tuple"""
    TFrom = np.array([]) #Transformed "From" coords
    for i in range(len(identicals)):
        PointID = identicals[i] #string
        From_ith = np.asarray(From[PointID]) # array
        R = Rotation_matrix(x[4],x[5],x[6])
        T = np.array([x[0],x[1],x[2]])
        TFrom_ith = T + x[3] * R.dot(From_ith)
        TFrom = np.append(TFrom,TFrom_ith)
    return TFrom

def Build_A(x,From,identicals):
    # x (TX,TY,TZ,q,alpha,beta,gamma)
    equation_count = 3 * len(identicals)
    RX = X_Rotation(x[4])
    RY = Y_Rotation(x[5])
    RZ = Z_Rotation(x[6])
    dRX = dX_Rotation(x[4])
    dRY = dY_Rotation(x[5])
    dRZ = dZ_Rotation(x[6])
    R = Rotation_matrix(x[4],x[5],x[6])
    A = np.zeros((equation_count,7))
    for i in range(len(identicals)):
        iii = 3*i
        PointID = identicals[i] #string
        From_ith = np.asarray(From[PointID]) # array
        A[iii,0] = 1
        A[iii+1,1] = 1
        A[iii+2,2] = 1
        D_alpha = dRX.dot(RY).dot(RZ).dot(From_ith)
        D_beta = RX.dot(dRY).dot(RZ).dot(From_ith)
        D_gamma = RX.dot(RY).dot(dRZ).dot(From_ith)
        D_q = R.dot(From_ith)
        A[iii,3] = D_q[0]
        A[iii+1,3] = D_q[1]
        A[iii+2,3] = D_q[2]
        A[iii,4] = D_alpha[0]
        A[iii+1,4] = D_alpha[1]
        A[iii+2,4] = D_alpha[2]
        A[iii,5] = D_beta[0]
        A[iii+1,5] = D_beta[1]
        A[iii+2,5] = D_beta[2]
        A[iii,6] = D_gamma[0]
        A[iii+1,6] = D_gamma[1]
        A[iii+2,6] = D_gamma[2]
    return A

def Helmert_LSM(From,To):
    R0, x0 = Helmert_aproximate_parameters(From,To)
    print(x0)
    identicals = list(set(To.keys()) & set(From.keys()))
    equation_count = 3*len(identicals)
    dx = np.zeros(7)
    x = np.array(x0)
#    To_ith = To[PointID] #tuple
    vI = np.empty((equation_count))
    vII = np.empty((equation_count))
    To_array = np.empty((equation_count)) # Transposed "To" coordinates
    Transformed_From = 1
    To_array = np.empty(0)
    From_array = np.empty(0)
    for i in range(len(identicals)):
        PointID = identicals[i] #string
        To_ith = To[PointID] #tuple
        From_ith = From[PointID] #tuple
        To_array = np.concatenate((To_array,np.array(To_ith)),axis = None)
        From_array = np.concatenate((From_array,np.array(From_ith)),
                                                                   axis = None)
    threshold = 0.000001 #basic unit
    metric = threshold + 1
    while metric > threshold:
        dx = np.empty(7)
        x += dx
        vI = Build_A(x,From,identicals).dot(dx) + vII
        metric = 0.000000000001
    return x,Transformed_From

T,R,q,Transformed_From = Helmert_LSM(From,To)