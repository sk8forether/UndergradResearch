# Peter Vaillancourt
# Creating orbits from e and a
# Rotates orbits in 3D 

from __future__ import division
from pylab import *
from numpy import *
import scipy.optimize as so
from mpl_toolkits.mplot3d import Axes3D
from functions import *

def orbit():
        e,a = pick_ea()
        P = sqrt( (4*pi*pi*a*a*a)/M ) # Period of the orbit
        t = random.random()*P # T = 0
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
        alpha = random.random()
        beta = random.random()
        gamma = random.random()
        new = alpha*dot(Rx,new) + beta*dot(Ry,new) + gamma*dot(Rz,new)
        return x,y,z,new #P,t

number = 50 
X = empty([number,len(theta)])
Y = empty([number,len(theta)])
Z = empty([number,len(theta)])
NewX = empty([number,len(theta)])
NewY = empty([number,len(theta)])
NewZ = empty([number,len(theta)])
fig = figure(1)
ax = fig.add_subplot(111, projection='3d')
for j in range(0,number):
        x,y,z,new = orbit()
        X[j,:] = x
        Y[j,:] = y
        Z[j,:] = z
        NewX[j,:] = new[0,:]
        NewY[j,:] = new[1,:]
        NewZ[j,:] = new[2,:]
        ax.plot_wireframe(x,y,z)
        axis('equal')
savefig('original.png')
        
fig = figure(2)
ax = fig.add_subplot(111, projection='3d')
for k in range(0,number):
        ax.plot_wireframe(NewX[k,:],NewY[k,:],NewZ[k,:])
axis('equal')
title('Rotated orbit in 3D')
savefig('rotated.png')




