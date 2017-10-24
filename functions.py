# Peter Vaillancourt 
# A collection of functions taylored to the model 

from __future__ import division
from pylab import *
from numpy import *
import math as m
import scipy.optimize as so
from mpl_toolkits.mplot3d import Axes3D

random.seed()

G = 6.67e-11 # N*(m/kg)^2
theta = linspace(0,2*pi,1000) #radians
M =  2e42  #kg, arbitrary - approximate mass of our galaxy
arange = 100. # Upper limit on semi-major axis, arbitrary
dtheta = 0.001 #"sensitivity" of the derivative
expansion = 1.
LOS = array([0.,1.,0.])
redshift = expansion*LOS

def pick_ea():
        e = random.random() # Random eccentricity
        a = random.uniform(1,arange) # Random semi-major axis
        return e,a

def polar(r):
        x = r*cos(theta)
        y = r*sin(theta)
        return (x,y)

def randphi():
        phi = random.uniform(0,2*pi) 
        return phi 

def funk(e,E,P,t):
        T = 0
        return E-e*sin(E)-2*pi*(t-T)/P

def dfunk(e,E):
        return 1-e*cos(E)

def eccanom(e,a,P,t):
        EA_0 = pi # initial guess
        TA = 0. # intialize true eccentric anomaly
        EA = so.newton(lambda EA_0: funk(e,EA_0,P,t),EA_0,fprime = lambda EA_0:
        dfunk(e,EA_0),tol=1e-09, maxiter=500000) 

        return EA #TA

def trueanom(EA,e):
        return 2*m.atan2(sqrt(1+e)*tan(EA/2),sqrt(1-e))

#instead, calculate TA, then make sure TA is in the range 0 to 2pi
def PutInRange(V,ll,uu):
# V is value to put into range.
#ll lower limit, eg 0
#uu upper limit, eg 2*pi
  if ll <= V <= uu:
    return V
  elif V > uu:
    while V > uu:
     V=V-uu
    return V 
  elif V < ll:
    while V < ll:
      V=V+uu
    return V
  else:
    return 0	  
  return TA


def polar2(r,angle):
        x = r*cos(angle)
        y = r*sin(angle)
        return (x,y)

def magnitude(V):
        mag = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
        return mag

def visviva(r,a):
        v = sqrt( G*M*((2/r) - (1/a)) )
        return v

def tangent(e,a,EA,r_spot):
        r1 = r_spot # r at spot
        r2 = a*(1.-(e*cos(EA+dtheta))) # r at an increment dtheta in front of spot
        dr = (a*e*(1-e*e)*sin(EA))/((1+e*cos(EA))*(1+e*cos(EA))) #derivative of r(TA)
        X1 = array([r1*cos(EA),r1*sin(EA),0]) #position vector
        X2 = array([r2*cos(EA+dtheta),r2*sin(EA+dtheta),0]) #position vector delta theta
        #V = array([ (dr*cos(EA) - r1*sin(EA)), (dr*sin(EA) + r1*cos(EA)), 0]) #radial velocity vector
        Tangent=X2-X1 #tangent vectors calculated my way.
                #p[i]=inner(X1[i,:],Tangent[i,:]) #test with inner product 	
                #p2[i]=inner(X1[i,:],V[i,:])
        NormT=Tangent/magnitude(Tangent)
#        return Tangent,p,p2
        return NormT


