# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:23:20 2021

@author: jbarker
"""

import numpy as np
import math
from angle import Angle as a

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
    R = X_Rotation(angle_tuple[0]) @ Y_Rotation(angle_tuple[1]) \
                                                @ Z_Rotation(angle_tuple[2])
    return R

def Transformation(x, From):
    """3D Helmert transformation with known transformation Key
    From is a dictionary of points
    (Rotation matrix parameters, Translation vector and scale in tuple)"""
    T = np.array(x[0:3])
    q = float(x[3])
    R = Rotation_matrix(x[-3:])
    From_transformed = {}
    for point in From:
        From_transformed[point] = tuple(T + q * R @ np.array(From[point][:3]))
    return From_transformed

def Helmert_aproximate_parameters(From,To):
    identicals = list(set(To.keys()) & set(From.keys()))
    if len(identicals) > 3:
        #MAKE BETTER CHOICE ON POINTS WHEN YOU HAVE TIME, IF YOU WANT... PLEASE
        point1_To = np.array(To[identicals[0]][:3])
        point2_To = np.array(To[identicals[-1]][:3])
        point3_To = np.array(To[identicals[len(identicals)//2]][:3])
        point1_To_original = point1_To.copy().transpose()
        #// truncating division = rounds result to lower int
        point1_From = np.array(From[identicals[0]][:3])
        point2_From = np.array(From[identicals[-1]][:3])
        point3_From = np.array(From[identicals[len(identicals)//2]][:3])
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
        FirstZrotation2_To = Z_Rotation(psi_To) @ point2_To.transpose()
        FirstZrotation3_To = Z_Rotation(psi_To) @point3_To.transpose()
        FirstZrotation2_From = Z_Rotation(psi_From) @ point2_To.transpose()
        FirstZrotation3_From = Z_Rotation(psi_From) @ point3_To.transpose()
        # Calculating second rotation angle from the first rotated coordinates
        fi_To = math.atan2(FirstZrotation2_To[1],FirstZrotation2_To[2])
        fi_From = math.atan2(FirstZrotation2_From[1],FirstZrotation2_From[2])
        # Applying the calculated angle to point 3
        SecondYrotation3_To = Y_Rotation(fi_To)@ FirstZrotation3_To.transpose()
        SecondYrotation3_From = Y_Rotation(
                                   fi_From) @ FirstZrotation3_From.transpose()
        # Calculating third rotation angle
        theta_To = math.atan2(SecondYrotation3_To[1],SecondYrotation3_To[0])
        theta_From = math.atan2(SecondYrotation3_From[1],
                                SecondYrotation3_From[0])
        # Using all three angles, the full rotation matrix for From To is made
        R_To = Z_Rotation(theta_To) @ Y_Rotation(fi_To) @ Z_Rotation(psi_To)
        R_From = Z_Rotation(theta_From) @ Y_Rotation(fi_From) @ Z_Rotation(
                                                                      psi_From)
        # The translation vector is calculated by rotating the original
        # point1_From to the "To" coordinate frame. By substracting the rotated
        # point1_From from point1_To we get the translation vector.
        R0 = R_To.transpose() @ R_From
        Translation = tuple(point1_To_original - R0 @ point1_From_original)
        # Euler rotation angles
        alpha = math.atan2(R0[1,2],R0[2,2])
        beta = - math.atan(R0[0,0])
        gamma = math.atan2(R0[0,1],R0[0,0])
        R_angles = (alpha, beta, gamma)
        x0 = Translation + (1.0,) + R_angles
    else:
        print('Not enough identical points for transformation calculations.')
    return R0, x0# RX + T

def Build_TFrom(x,From,identicals):
    """x are the transform parameters T,q,R in a tuple"""
    TFrom = np.array([]) #Transformed "From" coords
    for i in range(len(identicals)):
        PointID = identicals[i] #string
        From_ith = np.asarray(From[PointID][:3]) # array
        R = Rotation_matrix(x[4:])
        T = np.array([x[0:3]])
        TFrom_ith = T + x[3] * R @ From_ith
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
    R = Rotation_matrix(x[4:7])
    A = np.zeros((equation_count,7))
    for i in range(len(identicals)):
        iii = 3*i
        PointID = identicals[i] #string
        From_ith = np.asarray(From[PointID][:3]) # array
        A[iii,0] = 1
        A[iii+1,1] = 1
        A[iii+2,2] = 1
        D_alpha = dRX @ RY @ RZ @ From_ith
        D_beta = RX @ dRY @ RZ @ From_ith
        D_gamma = RX @ RY @ dRZ @ From_ith
        D_q = R @ From_ith
        A[iii:iii+3,3] = D_q
        A[iii:iii+3,4] = D_alpha
        A[iii:iii+3,5] = D_beta
        A[iii:iii+3,6] = D_gamma
    return A

def Helmert_transform(From,To):
    R0, x0 = Helmert_aproximate_parameters(From,To)
#    print(x0)
    identicals = list(set(To.keys()) & set(From.keys()))
    equation_count = 3*len(identicals)
    dx = np.zeros(7)
    x = np.array(x0)
    vI = np.empty((equation_count))
    vII = np.empty((equation_count))
    To_array = np.empty(0)
    TFrom = Build_TFrom(x,From,identicals)
    A = Build_A(x,From,identicals)
    for i in range(len(identicals)):
        PointID = identicals[i] #string
        To_ith = To[PointID][:3] #tuple
        To_array = np.concatenate((To_array,np.array(To_ith)),axis = None)
    threshold = 0.000000001 #fraction of basic unit
    metric = threshold + 1
    counter = 0
    while (metric > threshold) and (counter < 1000):
        l_prime = TFrom - To_array
        dx = -np.linalg.inv(A.transpose() @ A) @ A.transpose() @ l_prime
        x += dx
        vI = A @ dx + l_prime
        TFrom = Build_TFrom(x,From,identicals)
        vII = TFrom - To_array
        v = vI-vII
        metric = max(abs(v))
        counter += 1
        A = Build_A(x,From,identicals)
    if counter == 100:
        print("Too many iterations")
#    Transformed_From = Transformation(x,From)
    Trans_par = np.array([x[0], x[1], x[2], x[3],
                 a(x[4],a.T_RAD, True).angle,
                 a(x[5],a.T_RAD, True).angle,
                 a(x[6],a.T_RAD, True).angle])
    return Trans_par

# =============================================================================
# Testing data [m]
# =============================================================================
"""
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
"""

Nominal = {'Point1': (0.00, 0.00, 0.00),
      'Point2': (1000.00, 0.00, 0.00),
      'Point3': (0.00, 1000.00, 0.00),
      'Point4': (0.00, 0.00, 1000.00)
      }

X_45 = {'Point1': (0.00, 0.00, 0.00),
        'Point2': (1000.00, 0.00, 0.00),
        'Point3': (0.00, 760.41, -649.45),
        'Point4': (0.00, 649.45, 760.41)
      }

Y_45 = {'Point1': (0.00, 0.00, 0.00),
        'Point2': (760.41, 0.00, -649.45),
        'Point3': (0.00, 1000.00, 0.00),
        'Point4': (649.45, 0.00, 760.41)
      }

Z_45 = {'Point1': (0.00, 0.00, 0.00),
        'Point2': (760.41, -649.45, 0.00),
        'Point3': (649.45, 760.41, 0.00),
        'Point4': (0.00, 0.00, 1000.00)
      }

LotR = {'Point1': (0.00, 0.00, 0.00),
        'Point2': (572.06, -572.06, -587.79),
        'Point3': (441.35, 818.73, -367.29),
        'Point4': (691.35, -49.31, 720.84)
      }


RoTr = {'Point1': (100.00, -200.00, 300.00),
        'Point2': (672.06, -772.06, -287.79),
        'Point3': (541.35, 618.73, -67.29),
        'Point4': (791.35, -249.31, 1020.84)
      }


parX_45 = Helmert_transform(X_45,Nominal)
NewX_45 = Transformation(parX_45,X_45)

#parY_45 = Helmert_transform(Y_45,Nominal)
#NewY_45 = Transformation(parY_45,Y_45)
#
#parZ_45 = Helmert_transform(Z_45,Nominal)
#NewZ_45 = Transformation(parZ_45,Z_45)

#parLotR = Helmert_transform(LotR,Nominal)
#NewLotR = Transformation(parLotR,LotR)
#
#parRoTr = Helmert_transform(RoTr,Nominal)
#NewRoTr = Transformation(parRoTr,RoTr)
# =============================================================================
# Testing data - Results
# =============================================================================
'''
T:   2.98027817e+06   1.37948616e+06   5.59696765e+06 [m]
q:   1.00016673
Rot: 0.204758250  -0.670278274   1.11095571 [rad]

{'Point1': (3970673.003288134, 1018563.7398536116, 4870369.178586198),
'Point2': (3970667.5735479398, 1018565.1948317708, 4870373.0091473022),
'Point3': (3970659.461105776, 1018571.2697034873, 4870377.8815089483),
'Point4': (3970654.6043181177, 1018577.5169809447, 4870380.0196550041),
'Point5': (3970650.0896756388, 1018580.5766302142, 4870382.7734877616),
'Point6': (3970646.09606406, 1018581.6829998301, 4870385.6206144663)}'''
