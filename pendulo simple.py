# -*- coding: utf-8 -*-

from pylab import *

itermax = 100 #iteraciones maximas
y0 = array([1,0,0,0]) # condiciones iniciales para y = [pos_x, pos_y, velocidad_x, velocidad_y]
z0 = 0.0 # tension
t0 = 0.0 # tiempo inicial
tf = 5.0 # tiempo final 
tol = 1.e-12 # tolerancia para metodo de newton
N = 400 # numero de particiones del intervalo
m = 1 # peso pendulo
l = 1 # longitud pendulo 
gra = -9.81 #gravedad

def f(y,z):
  f1 = y[2]
  f2 = y[3]
  f3 = -z/(l*m)*y[0] 
  f4 = -z/(l*m)*y[1] + gra
  return array([f1,f2,f3,f4])

def g(y,z):
  g1 = m*(y[2]*y[2] + y[3]*y[3]) - z*l + gra*m*y[1]
  return g1

def fy(y,z):
  return array([[0,0,1,0], [0,0,0,1], [-z/(l*m), 0, 0, 0], [0, -z/(l*m), 0, 0]])

def fz(y,z):
  return array([0,0,-y[0]/(l*m), -y[1]/(l*m)])

def gy(y,z):
  return array([0, gra*m, 2*m*y[2], 2*m*y[3]])

def gz(y,z):
  return -l

