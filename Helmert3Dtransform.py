# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:23:20 2021

@author: jbarker
"""

import numpy as np
import math

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

def Rotation_matrix(angles):
    R = X_Rotation(angles[0]) @ Y_Rotation(angles[1]) \
                                                @ Z_Rotation(angles[2])
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
        point2_From = point2_From[:3] - point1_From[:3]
        point3_From = point3_From[:3] - point1_From[:3]
        point1_From = (0,0,0)
        # Calculating the aproximate parameters 
        # (based on Angle between planes paper)
        # First rotation angle calculations
        psi_To = math.atan2(point2_To[1],point2_To[0])
        psi_From = math.atan2(point2_From[1],point2_From[0])
        # Applying the calculated angle to point 2 and 3
        FirstZrotation2_To = Z_Rotation(psi_To) @ point2_To.transpose()
        FirstZrotation3_To = Z_Rotation(psi_To) @ point3_To.transpose()
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

def Build_TFrom_dict(x,From,identicals):
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

def Build_TFrom(x,From,identicals):
    """x are the transform parameters T,q,R in a tuple"""
    TFrom = np.array([]) #Transformed "From" coords
    for i in range(len(identicals)):
        iii = 3*i
        From_ith = From[i][:3] # array
        R = Rotation_matrix(x[4:])
        T = np.array([x[0:3]])
        TFrom_ith = T + x[3] * R @ From_ith
        TFrom = np.append(TFrom,TFrom_ith)
#    print(TFrom)
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
#        PointID = identicals[i] #string
        From_ith = np.asarray(From[i][:3]) # array
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
    identicals = list(set(To.keys()) & set(From.keys()))
    equation_count = 3*len(identicals)
    dx = np.zeros(7)
    x = np.array(x0)
    vI = np.empty((equation_count))
    vII = np.empty((equation_count))
    To_array = np.empty(0)
    TFrom = Build_TFrom_dict(x,From,identicals)
    for i in range(len(identicals)):
        PointID = identicals[i] #string
        To_ith = To[PointID][:3] #tuple
        To_array = np.concatenate((To_array,np.array(To_ith)),axis = None)
    threshold = 0.000001 #fraction of basic unit
    metric = threshold + 1
    counter = 0
    while (metric > threshold) and (counter < 10):
        A = Build_A(x,np.reshape(TFrom,(len(identicals),3)),identicals)
        l_prime = TFrom - To_array
        dx = -np.linalg.inv(np.transpose(A) @ A) @ np.transpose(A) @ l_prime
        x += dx
        vI = A @ dx + l_prime # First corrections
        TFrom = Build_TFrom(x,np.reshape(TFrom,(len(identicals),3)),identicals)
        vII = TFrom - To_array # Second corrections
        print(vI,vII)
        metric = max(abs(vII-vI)) #max(abs(vII-vI)) # 
        print(metric)
        counter += 1
    if counter > 9:
        print("Too many iterations")
    Transformed_From = Transformation(x,From)
    return x,Transformed_From

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
To = {
"Point0" : (0.0, 0.0, 0.0),
"Point1" : (1.0, 0.0, 0.0),
"Point2" : (0.0, 1.0, 0.0),
"Point3" : (0.0, 0.0, 1.0),
"Point4" : (1.0, 2.0, 3.0),
"Point5" : (4.0, 5.0, 6.0),
"Point6" : (7.0, 8.0, 9.0),
"Point7" : (7.01, 8.01, 9.01),
"Point8" : (-34.0, -83.0, -91.0),
}

From = {
"Point0" : (-8.780027463096536, -2.720601476546373, -24.573464270144104),
"Point1" : (-8.429707044341491, -3.6567499586382843, -24.46398378079296),
"Point2" : (-7.897140498768689, -2.410695613428955, -24.870129716134883),
"Point3" : (-8.524688128230633, -2.5362355448093794, -23.646537459579648),
"Point4" : (-5.959344316325122, -2.4650414723865803, -22.273224988441825),
"Point5" : (-1.5184015094523113, -3.6723193600871027, -20.036236357176854),
"Point6" : (2.882002148199907, -4.8793938866627755, -17.820044628943787),
"Point7" : (2.8869739874161393, -4.898918508351885, -17.80203952689806),
"Point8" : (-115.61731391068102, -15.249107085839404, -88.78855881835221),
}
# =============================================================================
# Testing data - Results
# ============================================================================
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

#To = {'Girder_17': (3580.033, 319.23, -450.0), 'Girder_15': (3200.305, 319.177, -450.0), 'Girder_11': (2266.398, 319.081, -450.0), 'Girder_7': (1606.962, 319.299, -450.0), 'Girder_3': (702.931, 319.295, -450.0), 'Girder_1': (323.147, 318.65, -450.0), 'Girder_2': (323.225, -319.282, -450.0), 'Girder_4': (702.676, -318.731, -450.0), 'Girder_8': (1607.732, -319.047, -450.0), 'Girder_12': (2267.506, -319.204, -450.0), 'Girder_16': (3200.582, -318.769, -450.0), 'Girder_18': (3580.677, -318.683, -450.0), 'Girder_13': (2744.555, 319.123, -259.682), 'Girder_5': (1169.555, 319.067, -259.725), 'Girder_6': (1169.555, -318.987, -259.516), 'Girder_14': (2744.555, -318.886, -259.257), 'Girder_10': (1932.286, -318.933, -445.804), 'Girder_9': (1931.985, 319.004, -446.523), 'PQK36_1': (3580.129, 319.168, -259.561), 'PQK36_2': (3580.392, 319.0, 260.176), 'PQK36_3': (3580.851, 190.734, 381.651), 'PQK36_4': (3580.31, -189.781, 381.618), 'PQK36_5': (3579.981, -319.023, 260.177), 'PQK36_6': (3580.49, -318.774, -258.959), 'PQK36_7': (3200.33, 319.143, -259.602), 'PQK36_8': (3200.397, 319.051, 260.128), 'PQK36_9': (3200.797, 190.148, 381.562), 'PQK36_10': (3200.68, -190.38, 381.225), 'PQK36_11': (3199.917, -319.228, 259.852), 'PQK36_12': (3200.404, -318.892, -259.514), 'PQK62_1': (703.196, 319.314, -260.024), 'PQK62_2': (703.923, 319.366, 260.716), 'PQK62_3': (703.392, 189.911, 381.514), 'PQK62_4': (702.991, -189.47, 381.848), 'PQK62_5': (703.82, -318.919, 260.03), 'PQK62_6': (702.983, -318.782, -259.336), 'PQK62_7': (323.31, 318.828, -259.505), 'PQK62_8': (323.754, 319.315, 260.722), 'PQK62_9': (322.645, 190.05, 381.429), 'PQK62_10': (322.986, -190.123, 381.216), 'PQK62_11': (323.408, -319.094, 260.357), 'PQK62_12': (323.274, -319.232, -259.953), 'PQL6_1': (2266.658, 319.01, -260.071), 'PQL6_2': (2267.37, 318.814, 260.553), 'PQL6_3': (2267.308, 189.988, 381.281), 'PQL6_4': (2267.461, -189.958, 381.464), 'PQL6_5': (2267.16, -318.638, 260.34), 'PQL6_6': (2267.413, -319.052, -259.248), 'PQL6_7': (1607.165, 319.106, -259.459), 'PQL6_8': (1607.718, 318.579, 260.443), 'PQL6_9': (1607.669, 189.952, 381.423), 'PQL6_10': (1607.666, -190.164, 381.566), 'PQL6_11': (1607.678, -318.455, 259.748), 'PQL6_12': (1607.717, -318.888, -259.318)}
#From = {'PQK36_4': (-1864.7391985337786, -1412.5845681076676, 337.5129251994295, 0.008484731470596056, 0.006427402985937218, 0.0015358679488451294), 'PQL6_6': (-563.6822583192709, -1636.9194273603023, -301.1856214205907, 0.0034206517670883666, 0.009933339443147556, 0.0018277972469305753), 'PQK62_11': (1308.9304639165696, -2155.815215165814, 221.20123913470763, 0.00557323317712973, 0.009179096480500828, 0.0009421056061386079), 'PQK62_5': (942.252947958758, -2054.6074218479876, 220.31446634117754, 0.004456264337266509, 0.009716910937035155, 0.0010421649035059612), 'PQK62_6': (944.0385231639474, -2054.147867656249, -299.0300855126724, 0.004447064969496911, 0.009676346240253202, 0.0014087881127118277), 'PQK62_12': (1310.1064657966838, -2154.92028065426, -299.11910328778254, 0.005561624388256889, 0.009147955442434755, 0.0012700024219439107), 'PQK62_10': (1274.7635284591167, -2280.4113882071556, 341.7819795813937, 0.005222980013401145, 0.00934327459497035, 0.001400528200078131), 'PQK62_4': (908.3305069086985, -2179.7601403242315, 341.860375875407, 0.004095392729722672, 0.0098277752489944, 0.0015414840951922001), 'PQK62_3': (806.830912644317, -2545.2781093678836, 340.7679868916535, 0.0032384119172258713, 0.010215862609070306, 0.0013679155418259586), 'Girder_5': (334.5118195376464, -2547.217877625032, -304.03077270736935, 0.0013951196466667477, 0.010622091958219249, 0.0012680327926292093), 'Girder_6': (503.3121340146122, -1930.035362377214, -304.746965806248, 0.0026697994678817004, 0.010237506953488414, 0.0016166062301026345), 'PQL6_10': (36.58919024011306, -1937.9792688962732, 340.28757472149886, 0.0001999117287525675, 0.010532217906334013, 0.0018494546552070045), 'PQL6_11': (71.00914002153476, -1814.131411508509, 218.73670545097795, 0.00041499645057480137, 0.010589549100618103, 0.0012769840459721886), 'PQL6_12': (72.07039018287901, -1812.8993074333082, -300.3212318963134, 0.000418870850903713, 0.010524223661527032, 0.0017435390102903954), 'Girder_10': (-238.6981197284123, -1726.0549631030801, -496.83050346103494, 0.0014060054268077657, 0.010166021003826074, 0.0029262696837823146), 'PQL6_4': (-599.375734960091, -1762.3231107437275, 339.2348623602857, 0.003384419925984741, 0.009950930336624687, 0.0019155919016397517), 'PQL6_5': (-564.5588090963123, -1638.1859087703535, 218.3813577277245, 0.003446658061159245, 0.010001052381830794, 0.0013333594038629074), 'Girder_14': (-1012.5429340518019, -1513.040647466852, -304.6923997474439, 0.005856660572607021, 0.008751565925990528, 0.0017624836477635015), 'Girder_13': (-1185.3143477882425, -2129.7260215327265, -307.6879265322543, 0.005195649787168174, 0.009335279614276566, 0.001348877585422272), 'PQK36_12': (-1462.966608342805, -1388.3937374205786, -302.78760842258276, 0.007679522774037704, 0.007288071794799409, 0.0015895515779077473), 'PQK36_10': (-1498.7070188344346, -1513.1937054501027, 337.6672274985987, 0.0074529925263003465, 0.007525033489584116, 0.0016793320246100253), 'Girder_TP': (-1502.2728809673715, -2013.3261783342223, -504.2373155014483, 0.006323252752681873, 0.008474316184378654, 0.002122505283184329), 'PQK36_5': (-1829.7463318971477, -1287.92193683412, 216.33985656689072, 0.008739478316604151, 0.006151563631458317, 0.0010335368257714105), 'PQK36_6': (-1829.3170688530968, -1287.21589756316, -302.8024710258846, 0.008702893151010147, 0.006123890554077427, 0.0014407290796265384), 'PQK36_9': (-1600.2496449549192, -1879.9153371503346, 337.20649843174596, 0.00691926129491641, 0.008128486026675614, 0.0014582005298869336), 'Girder_18': (-1829.3398830901938, -1283.5020424633437, -489.94690957692603, 0.008590413304673802, 0.006027227780747151, 0.002300842557097709), 'Girder_12': (-563.2926950899798, -1635.8549201831534, -491.1499923606186, 0.0033418674818390544, 0.00970494977779732, 0.0029138845219463276), 'Girder_8': (71.84568961356774, -1813.5629650638004, -489.42774964436546, 0.00040874513531349334, 0.010305488511270723, 0.0027812217845583883), 'Girder_4': (944.1641588237248, -2056.1063581024837, -488.9871486370865, 0.004383704073475598, 0.009546302505903763, 0.002270420029525201), 'Girder_2': (1308.93177535996, -2158.246795570055, -488.08931509879767, 0.005491349168622281, 0.009054425289329917, 0.0020477854216283606)}
x,Points = Helmert_transform(From,To)
#print(x,Points)
