# -*- coding: utf-8 -*-
"""
Created on Fri May 21 11:21:32 2021

@author: jbarker
"""

#import scipy as sp
import numpy as np
#import sys
#import string
import math as m
import config as cg
import functions as fc
import MainCode as MC

#from operator import itemgetter
#from collections import namedtuple

# Are there two epochs to calculate?
if len(cg.Which_epochs)>1:
    Two_epochs = True
else:
    Two_epochs = False
    
Nominal_coords = MC.Nominal_coords   
LoS_measurements = MC.LoS_measurements
Pol_measurements = MC.Pol_measurements
Pol_measurements_cart = MC.Pol_measurements_cart
StDevs_IFM_measurements = MC.StDevs_IFM_measurements
measured_distances_in_lines = MC.measured_distances_in_lines
sorted_measured_points_in_lines = MC.sorted_measured_points_in_lines

if Two_epochs:
    LoS_measurements_E1 = MC.LoS_measurements_E1
    Pol_measurements_E1 = MC.Pol_measurements_E1
    Pol_measurements_cart_E1 = MC.Pol_measurements_cart_E1
    StDevs_IFM_measurements_E1 = MC.StDevs_IFM_measurements_E1
    measured_distances_in_lines_E1 = MC.measured_distances_in_lines_E1
    sorted_measured_points_in_lines_E1 = MC.sorted_measured_points_in_lines_E1


    