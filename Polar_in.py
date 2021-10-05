# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 15:15:24 2021

@author: jbarker
"""

import re
import config as cg
import math as m
import numpy as np
from functions import StDev_sys_ppm




Measurements = Polar_2F_meas_read_in(cg.Epochs_dictionary['Pol'][1])
