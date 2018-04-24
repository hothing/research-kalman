###
Universidade Federal de Campina Grande
eROBOTICA UFCG
Author: Thiago de Freitas Oliveira Araújo
Antonio Marcus Nogueira Lima
Sample trajectory Kalman Filter Sample
using a P2dx(Mobile Robots Inc.)Robot model
2008
###
import numpy
import pylab
import time
import math
from numpy import *
from math import *
from scipy import *
 
# C-style declaration for didactic purposes
T = 50e-3
N = 1000
 
xi = 1
yi = 1
 
vd = 0.1
ve = 0.1222
B = 0.12
 
Q = 0.001*eye(2)
R = 50000*eye(2)
P = 1e6*eye(3)
 
xvR = 0
yvR = 0
pvR = 0
 
xvE = 0.0
yvE = 0.0
pvE = 0
 
u10 = []
u20 = []
 
xv0 = []
yv0 = []
pv0 = []
 
ri0 = []
oi0 = []
 
xE0 = []
yE0 = []
pE0 = []
rE0 = []
oE0 = []
 
# Robot Trajectory Definition
 
for i in range(N):
 
    if (i >= 1 and i < N/5):
        vd = 0.1
        ve = 0.1020
 
    elif (i >= N/5 and i < 2*N/5):
        vd = 0.2
        ve = 0.1960
 
    elif (i >= 2*N/5 and i < 3*N/5):
        vd = 0.1
        ve = 0.1020
 
    elif (i >= 3*N/5 and i <= 4*N/5):
        vd = 0.2
        ve = 0.1960
 
    elif (i >= 4*N/5 and i <= 5*N/5):
        vd = 0.1
        ve = 0.1020
 
    soma = (vd + ve)/2
    dif = vd - ve
    xvR = xvR + T*soma*math.cos(pvR)
    yvR = yvR + T*soma*math.sin(pvR)
    pvR = pvR + T*dif/B
 
    wr = rand(1,1)/50
    wo = rand(1,1)/50
    riR = math.sqrt(pow((xvR-xi),2) + pow((yvR-yi),2)) + wr
    oiR = math.atan((yvR-yi)/(xvR-xi)) + wo
    yR = bmat([[riR],[oiR]])
    xvE = xvE + T*soma*math.cos(pvE)
    yvE = yvR + T*soma*math.sin(pvE)
    pvE = pvE + T*dif/B
 
# Kalman Filter Calculus Steps
    F = mat([[1,0,(-T*(0.5*ve + 0.5*vd)*math.sin(pvE))],[0,1,(T*(0.5*ve + 0.5*vd)*math.cos(pvE))],[0,0,1]])
    G = mat([[0.5*T*math.cos(pvE),0.5*T*math.cos(pvE)],[0.5*T*math.sin(pvE),0.5*T*math.sin(pvE)],[T*1/(B),-T*1/(B)]])
    d0 = pow(xi,2) - 2*xi*xvE + pow(yi,2) - 2*yi*yvE + pow(xvE,2) + pow(yvE,2)
    d1 = math.sqrt(d0)
    H = mat([[-(xi-xvE)/d1,-(yi-yvE)/d1,0],[(yi-yvE)/d0,-(xi-xvE)/d0,-1]])
    I = mat([[1,0],[0,1]])
 
    riE = math.sqrt(pow((xvE-xi),2) + pow((yvE-yi),2))
    oiE = atan((yvE-yi)/(xvE-xi))
    yE = mat([[riE],[oiE]])
 
    P = F*P*transpose(F) + G*Q*transpose(G)
    e = yR - yE
    S = H*P*transpose(H) + I*R*transpose(I)
    W = P*transpose(H)*numpy.linalg.inv(S)
 
    xE = mat([[xvE],[yvE],[pvE]])
    print xE
    print 'WWW', W
    xE = xE + W*e
    xvE = xE[0,0]
    yvE = xE[1,0]
    pvE = xE[2,0]
 
    P = P - W*S*transpose(W)
 
    u10.append(vd)
    u20.append(ve)
 
    xv0.append(xvR)
    yv0.append(yvR)
    pv0.append(pvR)
 
    ri0.append(riR)
    oi0.append(oiR)
 
    xE0.append(xvE)
    yE0.append(yvE)
    pE0.append(pvE)
    rE0.append(riE)
    oE0.append(oiE)
 
pylab.xlabel("Xv")
pylab.ylabel("Yv")
 
pylab.plot(xv0, yv0, 'r-',label='Analysis of the Real Value')
pylab.plot(xE0, yE0, 'g-',label='Analysis of the Estimated Value')
pylab.xlabel('x Real')
pylab.ylabel('y Real')
pylab.legend()
pylab.grid(True)
 
pylab.show()