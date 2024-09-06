# -*- coding: utf-8 -*-
from pylab import *

def metodo_epsilon_encajado(f,fy,fz,g,gy,gz,x0,xf,N,y0,z0,tol,itermax,A,b):    
    # Paso de malla y condiciones iniciales
    h = (xf-x0)/float(N)
    y = zeros(N+1)
    z = zeros(N+1)
    t = linspace(x0,xf,N+1)
    y[0] = y0
    z[0] = z0
    s = int(len(b))
    
    identidad = eye(s)
    ceros = zeros((s,s))
    invA = linalg.inv(A)
    for k in range(0,N):        
        # Bucle para calcular los estados intermedios del metodo RK
        cont = 0
        d = 1.0 + tol
        yk_iter = y[k]*ones(s)
        zk_iter = z[k]*ones(s)
        uk_iter = block([yk_iter,zk_iter])
        
        while(d >= tol) and (cont < itermax):
            # Creamos los bloques de la Jacobiana
            A1 = identidad -h*(A*fy(yk_iter,zk_iter))
            A2 = -h*(A*fz(yk_iter,zk_iter))
            A3 = diag(gy(yk_iter,zk_iter))
            A4 = diag(gz(yk_iter,zk_iter))
            
            # Matriz Jacobiana de la funcion H a la que le buscamos un cero
            DH = block([[A1,A2],[A3,A4]])
            
            # Calculamos el valor de H
            H = block([yk_iter - y[k]*ones(s)- h*A@f(yk_iter,zk_iter), g(yk_iter,zk_iter)])
            uk_new = uk_iter - solve(DH, H)
            
            d = max(abs(uk_iter - uk_new))
            
            uk_iter = uk_new
            yk_iter = uk_iter[:s]
            zk_iter = uk_iter[-s:]

            cont += 1
        if cont == itermax:
            print('El método no va bien: numero maximo de iteraciones alcanzado')
        else:
            # Calculamos la siguiente compontente de y. El segundo sumando es la funcion incremento del metodo
            y[k+1] = y[k] + h*sum(b*f(yk_iter,zk_iter))
            z[k+1] = (1 - b@(invA@ones(s)))*z[k] + b@(invA@zk_iter)
    return (t,y,z)
