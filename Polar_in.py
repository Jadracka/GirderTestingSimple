# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 15:15:24 2021

@author: jbarker
"""

import re
import config as cg
import math as m
from functions import StDev_sys_ppm

"""IT IS NOT WORKING!!!!!"""


#cg.Epochs_dictionary['Pol'][0]
#cg.Print_2F_checks = True
def Polar_2F_meas_read_in(Meas_filename, 
                          Sd_StDev=(0.010,0), Hz_StDev=0.15, V_StDev=0.15):
    """Standard deviations for standard measurements with Leica AT960/AT930)"""
    Meas_file = open(Meas_filename,'r')
    max_diff_Hz = (Hz_StDev * m.sqrt(2))/1000
    max_index_err = (V_StDev * m.sqrt(2))/1000
    Measurements = {}
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
        words = re.split(';+|,+|\t+| +',row.strip())
        if words[0] not in Measurements.keys():
            Measurements[words[0]] = {}
            Measurements[words[0]][words[1]] =                         \
              (float(words[4]), float(words[2]), float(words[3]))
        else:
            if words[1] not in Measurements[words[0]].keys():
                Measurements[words[0]][words[1]] =                         \
                    (float(words[4]), float(words[2]), float(words[3]),1)
            else:
                Meas1 = Measurements[words[0]][words[1]][:3]
                Meas2 = (float(words[4]), float(words[2]), float(words[3]))
                avg_Sd = (Meas1[0]+Meas2[0])/2
                max_diff_Sd = StDev_sys_ppm(avg_Sd,Sd_StDev)
                diff_Sd = abs(Meas1[0]+Meas2[0])
                if max_diff_Sd > diff_Sd and cg.Print_2F_checks:
                    print("Point: %s, measured by %s fails 2Face check in"
                          "distance. Maximum difference is %1.3f mm and "
                          "measured difference is %1.3f mm."
                          %(words[1], words[0], max_diff_Sd, diff_Sd))
                if Meas2[2] >= 200 and Meas1[2] < 200:
                    diff_Hz = (Meas1[1] - (Meas2[1]-200))%400.0
                    index_err = (Meas1[2] + Meas2[2])%400.0
                                        
    del words, row
    Meas_file.close()
    return Measurements

Polar_2F_meas_read_in(cg.Epochs_dictionary['Pol'][0])