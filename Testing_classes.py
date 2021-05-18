# -*- coding: utf-8 -*-
"""
Created on Wed May  5 13:07:50 2021

@author: jbarker
"""
import functions as fc

class Point3DCart(dict):
    def __init__(self, PointID, X, Y, Z):
        self.PointID = PointID
        self.X = X
        self.Y = Y
        self.Z = Z
    def display(self):
        print(f'PointID = {self.PointID}, X = {self.X}, '
              f'Y = {self.Y}, Z = {self.Z}')
    def trans2Pol(self):
        return(self.PointID,fc.cart2polal3Dgon(self.X,self.Y,self.Z))
    
Point = Point3DCart('PQL6_12',138.13,12393.31,13.13)

class Point3DPol(dict):
    def __init__(self, PointID, Sd, Hz, Z):
        self.PointID = PointID
        self.Sd = Sd
        self.Hz = Hz
        self.Z = Z
    def display(self):
        print(f'PointID = {self.PointID}, Sd = {self.Sd}, '
              f'Hz = {self.Hz}, Z = {self.Z}')
    def trans2Cart(self):
        return(self.PointID,fc.polar2cart3Dgon(self.Sd,self.Hz,self.Z))