# -*- coding: utf-8 -*-
"""
Created on Tue Apr 13 11:12:15 2021

@author: jbarker
"""

#=============================================================================
# Overflow code
#=============================================================================
import config as cg
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import NewCode as nc

if cg.Print_FIDs:
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.axis('equal')
    xdata = []
    ydata = []
    zdata = []
    PointID = []
    for key in nc.Nominal_coords:
        PointID.append(key)
        xdata.append(nc.Nominal_coords[key][0])
        ydata.append(nc.Nominal_coords[key][1])
        zdata.append(nc.Nominal_coords[key][2])
#    ax.scatter3D(xdata, ydata, zdata);
    for i in range(len(PointID)): #plot each point + it's index as text above
        ax.scatter(xdata[i],ydata[i],zdata[i],color='b') 
        ax.text(xdata[i],ydata[i],zdata[i],  '%s' % (str(PointID[i])), size=9, zorder=1, color='k')
    plt.show()