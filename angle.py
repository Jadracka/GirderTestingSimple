# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 09:45:14 2021

@author: schloe
"""

import math as m
import numpy

class Angle():
    T_RAD=0
    T_GON=1
    T_DEG=2
    T_SELF_DEFINED=3

    minimum = 0
    maximum = m.tau
    interval = m.tau         ## maximum - minimum
    half_interval = m.pi     ## intervall / 2
    
    __INTERNAL_MINIMUM = -m.pi
    __INTERNAL_MAXIMUM = +m.pi
    __INTERNAL_INTERVAL = m.tau      ## __int_maximum - __int_minimum
    __INTERNAL_HALF_INTERVAL = m.pi  ## __int_intervall / 2

    angle=0
    a_type=T_RAD
    a_symmetric=False
    
    def __init__(self, angle=0, type=T_RAD, symmetric=False, minimum=0, maximum=m.tau):
        if not isinstance(angle, (float,numpy.longdouble,int)):
            raise TypeError("Type of angle has to be float or int")
        if not isinstance(symmetric, (bool)):
            raise TypeError("Type of symmetric has to be Boolean")
        if not (isinstance(minimum, (float, int)) and isinstance(maximum, (float,int))):
            raise TypeError("Type of minimum and maximum have to be float or int")
        if not (type in [self.T_RAD, self.T_GON, self.T_DEG, self.T_SELF_DEFINED]):
            raise TypeError("Type of type has to be Angle.t_rad, .t_gon, .t_deg or .t_self")
        
        if type==self.T_SELF_DEFINED:
            self.minimum = minimum
            self.maximum = maximum
        else:
            self.minimum = 0
            if type==self.T_RAD:
                self.maximum = m.tau
            elif type==self.T_GON:
                self.maximum = 400
            elif type==self.T_DEG:
                self.maximum = 360
            else:
                raise TypeError("How did you get here? This should never happen.")
        
        self.interval = self.maximum - self.minimum
        self.half_interval = self.interval / 2
        if ((symmetric==True) and (type!=self.T_SELF_DEFINED)):
            self.minimum -= self.half_interval
            self.maximum -= self.half_interval

        self.angle = self.__ext2int(self.__normalize(angle))
        self.a_type = type
        self.a_symmetric = symmetric

    def __str__(self):
        return (format(self.__normalize(self.__int2ext(self.angle))))
        
    def __internal_normalize(self, angle):
        if (angle < self.__INTERNAL_MINIMUM) or (angle >= self.__INTERNAL_MAXIMUM):
            angle = ((angle-self.__INTERNAL_MINIMUM)%(self.__INTERNAL_INTERVAL))+self.__INTERNAL_MINIMUM
        return(angle)

    def __normalize(self, angle):
        if (angle < self.minimum) or (angle >= self.maximum):
            angle = ((angle-self.minimum)%(self.interval))+self.minimum
        return(angle)
        
    def __ext2int(self,value):
        return(value/self.half_interval*self.__INTERNAL_HALF_INTERVAL)

    def __int2ext(self,value):
        return(self.__normalize(value*self.half_interval/self.__INTERNAL_HALF_INTERVAL))
        
    def get_intervall(self):
        return(self.intervall)
        
    def get_half_intervall(self):
        return(self.half_intervall)
        
    def __neg__(self):
        return(-self.__int2ext(self.angle))
        
    def __invert__(self):
        return(self.__int2ext(self.angle))

    def __add__(self, add):
        if isinstance(add, Angle):
            angle = self.angle + add.angle
        elif isinstance(add, (float, int)):
            angle = self.angle + self.__ext2int(add)
        else:
            return(None)
        return(self.__int2ext(angle))
        
    def __radd__(self,add):
        return(self.__add__(add))
        
    def __sub__(self, sub):
        if isinstance(sub, Angle):
            angle = self.angle - sub.angle
        elif isinstance(sub, (float, int)):
            angle = self.angle - self.__ext2int(sub)
        else:
            return(None)
        return(self.__int2ext(angle))
        
    def __rsub__(self, sub):
        if isinstance(sub, (float, int)):
            angle = self.__ext2int(sub) - self.angle
        else:
            return(None)
        return(self.__int2ext(angle))
        
    def __mul__(self, mul):
        if isinstance(mul, (float, int)):
            angle = self.angle * mul
        else:
            return(None)
        return(self.__int2ext(angle))

    def __rmul__(self, mul):
        return(self.__mul__(mul))
        
    def __truediv__(self, div):
        if isinstance(div, (float, int)):
            angle = self.angle / div
        else:
            return(None)
        return(self.__int2ext(angle))
        
    def __rtruediv__(self, div):
        return(None)
        
    def comp_error(self):
        raise TypeError("angles can only be compared to angles")
        return(None)

    def __lt__(self, other):
        if isinstance(other, Angle):
            if other.angle == self.angle:
                return(False)
            if (self.__internal_normalize(self.angle - other.angle)) < 0:
                return(True)
            else:
                return(False)
        else:
            self.comp_error()

    def __le__(self, other):
        if isinstance(other, Angle):
            if other.angle == self.angle:
                return(True)
            if (self.__internal_normalize(self.angle - other.angle)) < 0:
                return(True)
            else:
                return(False)
        else:
            self.comp_error()

    def __gt__(self, other):
        if isinstance(other, Angle):
            if other.angle == self.angle:
                return(False)
            if (self.__internal_normalize(self.angle - other.angle)) > 0:
                return(True)
            else:
                return(False)
        else:
            self.comp_error()

    def __ge__(self, other):
        if isinstance(other, Angle):
            if other.angle == self.angle:
                return(True)
            if (self.__internal_normalize(self.angle - other.angle)) > 0:
                return(True)
            else:
                return(False)
        else:
            self.comp_error()

    def __eq__(self, other):
        if isinstance(other, Angle):
            if self.__internal_normalize(other.angle - self.angle) == 0:
                return(True)
            else:
                return(False)
        else:
            self.comp_error()

    def __ne__(self, other):
        if isinstance(other, Angle):
            if other.angle != self.angle:
                return(True)
            else:
                return(False)
        else:
            self.comp_error()

    @staticmethod
    def is_similar(arg1, arg2, delta):
        if abs(arg1-arg2)<abs(delta):
            return(True)
        else:
            return(False)

    # all trigonometric functions are declared as staitic methods, so they
    # can be called as Angle.sin(x). x can be either of Angle-Type or any type
    # the math library provides (float, int). It's not possible to call these
    # methods as x.sin(), when x is of Angle type.

    @staticmethod
    def sin(x):
        if isinstance(x, Angle):
            return(m.sin(x.__ext2int(x)))
        else:
            return(m.sin(x))
    
    @staticmethod        
    def cos(x):
        if isinstance(x, Angle):
            return(m.cos(x.__ext2int(x)))
        else:
            return(m.cos(x))

    @staticmethod
    def tan(x):
        if isinstance(x, Angle):
            return(m.tan(x.__ext2int(x)))
        else:
            return(m.tan(x))

    def asin(self, arg, type=T_RAD, symmetric=False, minimum=0, maximum=m.tau):
        if not isinstance(arg, (float, int)):
            raise TypeError("Type of arg has to be float or int")
        self.__init__(0,type,symmetric,minimum,maximum)
        self.angle = m.asin(arg)
        return(self.__int2ext(self.angle))

    def acos(self, arg, type=T_RAD, symmetric=False, minimum=0, maximum=m.tau):
        if not isinstance(arg, (float, int)):
            raise TypeError("Type of arg has to be float or int")
        self.__init__(0,type,symmetric,minimum,maximum)
        self.angle = m.acos(arg)
        return(self.__int2ext(self.angle))

    def atan(self, arg, type=T_RAD, symmetric=False, minimum=0, maximum=m.tau):
        if not isinstance(arg, (float, int)):
            raise TypeError("Type of arg has to be float or int")
        self.__init__(0,type,symmetric,minimum,maximum)
        self.angle = m.atan(arg)
        return(self.__int2ext(self.angle))
        
    def atan2(self, arg1, arg2, type=T_RAD, symmetric=False, minimum=0, maximum=m.tau):
        if not (isinstance(arg1, (float, int)) and isinstance(arg2, (float, int))):
            raise TypeError("Type of arg1 and arg2 has to be float or int")
        self.__init__(0,type,symmetric,minimum,maximum)
        self.angle = m.atan2(arg1,arg2)
        return(self.__int2ext(self.angle))
      
        
"""        
w1=Angle(500,Angle.T_GON)
w2=Angle(200,Angle.T_GON)
w3=Angle(50,Angle.T_GON)
w4=Angle(0,Angle.T_DEG)
w4=w4+720

print('w4      =',w4)
print('w1      =',w1)
print('-w1      =',-w1)
print('~w1     =',~w1)
print('w2      =',w2)
print('w1+w2   =',w1+w2)
print('w1+30   =',w1+30)
print('w1+30.0 =',w1+30.0)
print('w1+\'30\' =',w1+'30')
print('30+w1   =',30+w1)
print('30.0+w1 =',30.0+w1)
print('\'30\'+w1 =','30'+w1)
print('w2-w1   =',w2-w1)
print('w1-30   =',w1-30)
print('w1-30.0 =',w1-30.0)
print('w1-\'30\' =',w1-'30')
print('30-w1   =',30-w1)
print('30.0-w1 =',30.0-w1)
print('\'30\'-w1 =','30'-w1)
print('w1*w2   =',w1*w2)
print('w1*3    =',w1*3)
print('w1*3.0  =',w1*3.0)
print('3*w1    =',3*w1)
print('-3.0*w1 =',-3.0*w1)
print('w1/3    =',w1/3)
print('3/w1    =',3/w1)
print('w1/0    = division by zero error')
print('sin(w1) =',Angle.sin(w1))
print('atan2(-3,2)=',w3.atan2(-3,2,Angle.T_GON))
print('20<19   =',Angle(20,Angle.T_GON)<Angle(19,Angle.T_GON))
print('19<20   =',Angle(19,Angle.T_GON)<Angle(20,Angle.T_GON))
print('19<19   =',Angle(19,Angle.T_GON)<Angle(19,Angle.T_GON))
print('19<=19  =',Angle(19,Angle.T_GON)<=Angle(19,Angle.T_GON))
print('50<251  =',Angle(50,Angle.T_GON)<Angle(251,Angle.T_GON))
print('50<250  =',Angle(50,Angle.T_GON)<Angle(250,Angle.T_GON),' ... undefined!')
print('20>19   =',Angle(20,Angle.T_GON)>Angle(19,Angle.T_GON))
print('19>20   =',Angle(19,Angle.T_GON)>Angle(20,Angle.T_GON))
print('19>19   =',Angle(19,Angle.T_GON)>Angle(19,Angle.T_GON))
print('19>=19  =',Angle(19,Angle.T_GON)>=Angle(19,Angle.T_GON))
print('50>251  =',Angle(50,Angle.T_GON)>Angle(251,Angle.T_GON))
print('50>250  =',Angle(50,Angle.T_GON)>Angle(250,Angle.T_GON),' ... undefined!')
print('50==450 =',Angle(50,Angle.T_GON)==Angle(450,Angle.T_GON))
print('50!=450 =',Angle(50,Angle.T_GON)!=Angle(450,Angle.T_GON))
print('50!=451 =',Angle(50,Angle.T_GON)!=Angle(451,Angle.T_GON))
print('w1!=w2  =',w1!=w2)
w3=Angle(50,Angle.T_GON)
print('w2==4*w3=',w2==Angle(w3*4,Angle.T_GON))
print('0.001 is_similar -0.001 by 0.1',Angle.is_similar(-0.001,+0.001,0.1))
print('w1 is_similar w2 by 0.1',Angle.is_similar(w1,w2,0.1))
w1=Angle(150,Angle.T_GON)
w2=Angle(135,Angle.T_DEG)
print(w1)
print(w2)
print('w1 is_similar w2 by 0.1',Angle.is_similar(w1,w2,0.1))
"""