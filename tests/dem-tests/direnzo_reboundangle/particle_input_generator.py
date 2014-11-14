#!/usr/bin/env python

from math import sin,cos,pi

radius = 0.250
density = 4.0

angles = [0.239726,
          2.294521,
          4.292237,
          9.942922,
          20.55936,
          30.37671,
          39.62329,
          49.78311,
          60.05708,]

for ii in range(0,len(angles)):
    angle = angles[ii]
    uu =  390*sin(angle*pi/180)
    vv = -390*cos(angle*pi/180)
    ww = 0.0
    print "%g %g %g %g %g %g %g %g" % (1.0,1.0,ii+1, radius, density, uu, vv, ww )
