# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 09:45:14 2021

@author: schloe
"""

import math as m


class Angle:
#     minimum = 0
#     maximum = 400
    minimum = 0
    maximum = 2*m.pi
#    minimum = -m.pi
#    maximum = +m.pi

    def __init__(self, angle=0):
        if isinstance(angle, (float,int)):
            Angle.intervall = self.maximum - self.minimum
            Angle.half_intervall = self.intervall / 2
            self.angle = float(self.normalize(angle))
        else:
            raise TypeError("Type of angle has to be either float or int")

    def __str__(self):
        return (format(self.angle))

    def normalize(self, angle):
        if (angle < self.minimum) or (angle >= self.maximum):
            angle = ((angle-self.minimum)%self.intervall)+self.minimum
        return(angle)

    def __neg__(self):
        return(-self.angle)

    def __add__(self, add):
        if isinstance(add, Angle):
            angle = self.angle + add.angle
        elif isinstance(add, (float, int)):
            angle = self.angle + add
        else:
            return(None)
        angle = self.normalize(angle)
        return(angle)

    def __radd__(self,add):
        return(self.__add__(add))

    def __sub__(self, sub):
        if isinstance(sub, Angle):
            angle = self.angle - sub.angle
        elif isinstance(sub, (float, int)):
            angle = self.angle - sub
        else:
            return(None)
        angle = self.normalize(angle)
        return(angle)

    def __rsub__(self, sub):
        if isinstance(sub, (float, int)):
            angle = sub - self.angle
        else:
            return(None)
        angle = self.normalize(angle)
        return(angle)

    def __mul__(self, mul):
        if isinstance(mul, (float, int)):
            angle = self.angle * mul
        else:
            return(None)
        angle = self.normalize(angle)
        return(angle)

    def __rmul__(self, mul):
        return(self.__mul__(mul))

    def __truediv__(self, div):
        if isinstance(div, (float, int)):
            angle = self.angle / div
        else:
            return(None)
        angle = self.normalize(angle)
        return(angle)

    def __rtruediv__(self, div):
        return(None)

    def sin(arg):
        return(m.sin(arg.angle/Angle.half_intervall*m.pi))

    def cos(arg):
        return(m.cos(arg.angle/Angle.half_intervall*m.pi))

    def tan(arg):
        return(m.tan(arg.angle/Angle.half_intervall*m.pi))

    def asin(arg):
        angle = m.asin(arg)/m.pi*Angle.half_intervall
        angle = Angle.normalize(Angle,angle)
        return(angle)

    def acos(arg):
        angle = m.acos(arg)/m.pi*Angle.half_intervall
        angle = Angle.normalize(Angle,angle)
        return(angle)

    def atan(arg):
        angle = m.atan(arg)/m.pi*Angle.half_intervall
        angle = Angle.normalize(Angle,angle)
        return(angle)

    def atan2(arg1,arg2):
        angle = m.atan2(arg1,arg2)/m.pi*Angle.half_intervall
        angle = Angle.normalize(Angle,angle)
        return(angle)


#==============================================================================
# Demonstration data
#==============================================================================
#w1=Angle(500)
#w2=Angle(200)
#w3=Angle(50)
#w4=Angle()
#w4=w4+800
#
#print('w4      =',w4)
#print('w1      =',w1)
#print('w2      =',w2)
#print('w1+w2   =',w1+w2)
#print('w1+30   =',w1+30)
#print('w1+30.0 =',w1+30.0)
#print('w1+\'30\' =',w1+'30')
#print('30+w1   =',30+w1)
#print('30.0+w1 =',30.0+w1)
#print('\'30\'+w1 =','30'+w1)
#print('w2-w1   =',w2-w1)
#print('w1-30   =',w1-30)
#print('w1-30.0 =',w1-30.0)
#print('w1-\'30\' =',w1-'30')
#print('30-w1   =',30-w1)
#print('30.0-w1 =',30.0-w1)
#print('\'30\'-w1 =','30'-w1)
#print('w1*w2   =',w1*w2)
#print('w1*3    =',w1*3)
#print('w1*3.0  =',w1*3.0)
#print('3*w1    =',3*w1)
#print('-3.0*w1 =',-3.0*w1)
#print('w1/3    =',w1/3)
#print('3/w1    =',3/w1)
#print('w1/0    = division by zero error')
#print('sin(w1) =',Angle.sin(w3))
#print('atan2(3,2)=',Angle.atan2(3,2))

