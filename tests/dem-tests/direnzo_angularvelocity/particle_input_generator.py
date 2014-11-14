#!/usr/bin/env python

from math import sin,cos,pi

radius = 0.250
density = 4.0

angles =[6.00246,
         10.98401,
         21.19311,
         31.09471,
         39.02829,
         50.2214,
         60.246]

for ii in range(0,len(angles)):
    angle = angles[ii]
    uu =  390*sin(angle*pi/180)
    vv = -390*cos(angle*pi/180)
    ww = 0.0
    print "%g %g %g %g %g %g %g %g" % (1.0,1.0,ii+1, radius, density, uu, vv, ww )
