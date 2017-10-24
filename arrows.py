# Peter Vaillancourt
# Generates tangent vectors for elliptical orbits in a plane

from __future__ import division
from pylab import *
from numpy import *
import scipy.optimize as so
from mpl_toolkits.mplot3d import Axes3D
from functions import *

random.seed()

G=1
M=1
a = 1 # arbitrary
e = 0.9 # high enough eccentricity that it matters
dtheta=0.001 # "sensitivity" of the derivative

#Initialize
X1=empty(3)
X2=empty(3)
V=empty(3)
Tangent=empty(3)
p=empty(len(theta))
p2=empty(len(theta))
NormT=empty(3)

def magnitude(V):
        mag = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
        return mag

def visviva(r,a):
        v = sqrt( G*M*((2/r) - (1/a)) )
        return v

def tangent(e,a,EA,r_spot):
        r1 = r_spot # r at spot
        r2 = a*(1.-(e*cos(EA+dtheta))) # r at an increment dtheta in front of spot
        x,y = polar2(r2,EA)
        dr = (a*e*(1-e*e)*sin(EA))/((1+e*cos(EA))*(1+e*cos(EA))) #derivative of r(TA)
        X1 = array([r1*cos(EA),r1*sin(EA),0]) #position vector
        X2 = array([r2*cos(EA+dtheta),r2*sin(EA+dtheta),0]) #position vector delta theta
        Tangent=X2-X1 #tangent vectors calculated my way.
        NormT=Tangent/magnitude(Tangent)
        return NormT,x,y

def orbit():
        #Pre-written stuff for orbital velocity
        e,a = pick_ea()
        #e = 0.9 # for science!
        #a = 1 # for science!
        P = sqrt( (4*pi*pi*a*a*a)/(G*M) ) # Period of the orbit
        R = random.random()
        t = R*P # T = 0
        EA = eccanom(e,a,P,t)
        TA = PutInRange(trueanom(EA,e),0,2*pi)
        r = (a*(1-e*e))/(1+e*cos(theta))
        x_rad,y_rad = polar(r)
        v = visviva(r,a) 
        r_spot = a*(1.-(e*cos(EA)))
        v_spot = visviva(r_spot,a) #vis-viva velocity at spot
        x_spot,y_spot = polar2(r_spot,TA) # spot on the orbit
        # Tangent Vector
        Tangent,vector1,vector2 = tangent(e,a,EA,r_spot)
        radial = v_spot*Tangent
        test = magnitude(radial)
        plot(x_rad,y_rad,x_spot,y_spot,'p')
        arrow(x_spot,y_spot,vector1,vector2)
        savefig('Diagram.png')
        return P,t,v,radial

number = 5
expansion = 1.
LOS = array([0.,1.,0.])
redshift = expansion*LOS
Radial = empty(number)
for i in range(0,len(Radial)):
        P,t,v,Rad = orbit()
        Rad = Rad + redshift
        Radial[i] = inner(Rad,LOS)

#savefig('histogramY.png')

