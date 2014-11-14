#!/usr/bin/env python

from math import sin,cos,pi

radius = 0.250
density = 4.0

angles = [2.075099,
          3.913043,
          9.664032,
          20.09881,
          30.23715,
          39.3083,
          49.56522,
          59.88142,]

for ii in range(0,len(angles)):
    angle = angles[ii]
    uu =  390*sin(angle*pi/180)
    vv = -390*cos(angle*pi/180)
    ww = 0.0
    print "%g %g %g %g %g %g %g %g" % (1.0,1.0,ii+1, radius, density, uu, vv, ww )
