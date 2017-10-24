# Peter Vaillancourt
# Creates orbits in 3D, places spot, finds velocities, rotates orbit, rotates spot, 

from __future__ import division
from pylab import *
from numpy import *
import scipy.optimize as so
from mpl_toolkits.mplot3d import Axes3D
from functions import *

number = 50 

#Initialize
X1=empty(3)
X2=empty(3)
V=empty(3)
Tangent=empty(3)
NormT=empty(3)
X = empty([number,len(theta)])
Y = empty([number,len(theta)])
Z = empty([number,len(theta)])
NewX = empty([number,len(theta)])
NewY = empty([number,len(theta)])
NewZ = empty([number,len(theta)])

def orbit():
        e,a = pick_ea()
        P = sqrt( (4*pi*pi*a*a*a)/(G*M) ) # Period of the orbit
        t = random.random()*P # T = 0
        return e,a,P,t

# rotations in 3D
def rotations(galaxy):
        e,a,P,t = galaxy
        r = (a*(1-e*e))/(1-e*cos(theta))
        x,y = polar(r)
        z = zeros(len(x))
        r = (a*(1-e*e))/(1-e*cos(theta))
        qx,qy = polar(r)
        qz = zeros(len(qx))
        new = [qx, qy, qz]
        phi = randphi()
        Rx = array([[1, 0, 0],
                   [0, cos(phi), -sin(phi)],
                   [0, sin(phi),  cos(phi)]])
        Ry = array([[cos(phi), 0, sin(phi)],
                   [0, 1, 0],
                   [-sin(phi), 0, cos(phi)]])
        Rz = array([[cos(phi), -sin(phi), 0],
                   [sin(phi), cos(phi), 0],
                   [0, 0, 1]])
        new = dot(Rx,new) + dot(Ry,new) + dot(Rz,new)
        return x,y,z,new #P,t

#Velocities
def velocities(galaxy):
        e,a,P,t = galaxy
        r = (a*(1-e*e))/(1-e*cos(theta))
        EA = eccanom(e,a,P,t)
        TA = PutInRange(trueanom(EA,e),0,2*pi)
        r = (a*(1-e*e))/(1+e*cos(theta))
        v = visviva(r,a) 
        r_spot = a*(1.-(e*cos(EA)))
        v_spot = visviva(r_spot,a) #vis-viva velocity at spot
        x_spot,y_spot = polar2(r_spot,TA) # spot on the orbit
        # Tangent Vector
        Tangent = tangent(e,a,EA,r_spot)
        radial = v_spot*Tangent
        test = magnitude(radial)
        return radial

def setup(number):
        Galaxy = zeros(shape=(number,4))
        for q in range(0,number):
                Galaxy[q,:] = orbit()
        return Galaxy

Galaxy = setup(number)
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
for j in range(0,number):
        x,y,z,new = rotations(Galaxy[j])
        X[j,:] = x
        Y[j,:] = y
        Z[j,:] = z
        NewX[j,:] = new[0,:]
        NewY[j,:] = new[1,:]
        NewZ[j,:] = new[2,:]
        ax.plot_wireframe(x,y,z)
        axis('equal')
#savefig('original.png')
        
fig = figure(2)
ax = fig.add_subplot(111, projection='3d')
for k in range(0,number):
        ax.plot_wireframe(NewX[k,:],NewY[k,:],NewZ[k,:])
axis('equal')
title('Rotated orbit in 3D')
#savefig('rotated.png')

Radial = empty(number)
for i in range(0,number):
        Rad = velocities(Galaxy[i])
        Rad = Rad + redshift
        Radial[i] = inner(Rad,LOS)

fig = figure(3)
hist(Radial)
show()
#savefig('histogramY.png')
