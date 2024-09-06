# -*- coding: utf-8 -*-
from pylab import *

def h(u):
    return 1.e-6*(exp(u/0.026) - 1)

def dh(u):
    return exp(u/0.026)/26000


def Ue(t):
    return 0.4*sin(200*pi*t)

def dUe(t):
    return 80*pi*cos(200*pi*t)

Ub = 6

R0 = 1000
R1 = 9000
R2 = 9000
R3 = 9000
R4 = 9000
R5 = 9000

C1 = 1.e-6
C2 = 2.e-6
C3 = 3.e-6

itermax = 100
y0 = array([-Ub*R1/(R1+R2),Ub*R1/(R1+R2),Ub,0])
z0 = array([0,Ub])
t0 = 0.0
tf = 0.2
tol = 1.e-12
N = 1000

def f(y,z):
  f1 = (Ue(y[3]) - z[0])/(R0*C1)
  f2 = h(z[0] - y[0] - y[1])/C2 - y[1]/(C2*R3)
  f3 = (z[1] - y[2])/(C3*R5)
  f4 = 1
  return array([f1,f2,f3,f4])

def g(y,z):
  hval = h(z[0] - y[0] - y[1])
  g1 = (Ue(y[3]) - z[0])/R0 + Ub/R2 + (y[0] - z[0])*(1/R1 + 1/R2) - hval/100
  g2 = (Ub - z[1])/R4 - 0.99*hval + (y[2] - z[1])/R5
  return array([g1,g2])

def fy(y,z):
  dhval = dh(z[0] - y[0] - y[1])
  return array([[0, 0, 0, dUe(y[3])/(R0*C1)],
                [-dhval/C2, -dhval/C2 - 1/(C2*R3), 0, 0]
                , [0, 0, -1/(C3*R5), 0], 
                [0,0,0,0]])

def fz(y,z):
  return array([[-1/(R0*C1), 0],
                [dh(z[0] - y[0] - y[1])/C2, 0], 
                [0, 1/(C3*R5)],
                [0,0]])

def gy(y,z):
  dhval = dh(z[0] - y[0] - y[1])
  return array([[1/R1 + 1/R2 + dhval/100, dhval/100, 0, dUe(y[3])/R0],
                [0.99*dhval, 0.99*dhval, 1/R5, 0]])

def gz(y,z):
  dhval = dh(z[0] - y[0] - y[1])
  return array([[-1/R0 - 1/R1 - 1/R2 - dhval/100, 0],
                [-0.99*dhval, -1/R4 - 1/R5]])