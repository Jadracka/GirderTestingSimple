# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 15:27:51 2021

@author: jbarker
"""

"""just a small code to create test measurements"""

import config as cg
import NewCode as nc
import random as rnd
import functions as fc

out_file = open("Testing_measurements_Epoch1.txt","w+")

lines_of_measurements = ()
for line in cg.Lines_of_sight:
    distance_offset = round(rnd.uniform(1500, 3000), 5)
    H_rand = rnd.uniform(0,400)
    Z_rand = rnd.uniform(45,135)
    for i in range (1,len(cg.Lines_of_sight[line])):
        distance_error = rnd.gauss(0,0.001)
        H_error = rnd.gauss(0,0.00015)
        Z_error = rnd.gauss(0,0.00015)
        distance_nom = fc.slope_distance(nc.Nominal_coords[
                                               cg.Lines_of_sight[line][i-1]],
                                         nc.Nominal_coords[
                                                 cg.Lines_of_sight[line][i]
                                         ])
        H = H_rand + H_error
        Z = Z_rand + Z_error
        dist = distance_offset + distance_nom + distance_error
        H_res = round(H, 5)
        Z_res = round(Z, 5)
        dist_res = round(dist, 5)
        out_file.write(f"{line} {cg.Lines_of_sight[line][i]} {dist_res} {H_res} {Z_res}\n")
del i, distance_offset, line, H_error, H_rand, H, H_res, Z, Z_error,  Z_rand, Z_res, dist, dist_res, distance_error, distance_nom

out_file.close()