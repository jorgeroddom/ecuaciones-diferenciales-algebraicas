from pylab import *

def state_space_form(f,fy,fz,g,gy,gz,t0,tf,N,y0,z0,tol,itermax,A,b):
  # Creamos el paso de malla
  h = (tf - t0)/float(N)
  
  # Dimension de y y de z
  n = len(y0)
  m = 1

  # Definimos las matrices de la solucion 
  # numerica e incluimos las condiciones iniciales
  y = zeros([n,N+1])
  z = zeros(N+1)
  y[:,0] = y0
  z[0] = z0

  # Creamos el vector de tiempos
  t = linspace(t0,tf,N+1)

  # Calculamos el numero de etapas del metodo RK
  s = len(b)

  # Vamos utilizar mucho lo siguiente
  ns = n*s
  ms = m*s
  I = identity(s)
  E = -h*A
  
  def unir_matrices(A):
      s1, = shape(A)
      n1,m1 = shape(A[0]) 
      M = zeros((n1,m1), dtype='O')
      for i in range(n1):
          for j in range(m1):
              M[i,j] = zeros(s1)
              for k in range(s1):
                  M[i,j][k] = A[k][i,j]
      return M
  
  def unir_vector(A):
      s1, = shape(A)
      n1, = shape(A[0]) 
      M = zeros(n1, dtype='O')
      for i in range(n1):
              M[i] = zeros(s1)
              for k in range(s1):
                  M[i][k] = A[k][i]
      return M
  
  for k in range(0,N):
    cont = 0
    d = 1.0 + tol
    # Creamos los y_k^(1),...,y_k^(s),
    # z_k^(1),...,z_k^(s) iniciales
    yk_iter = y[:,k]*ones((s,n))
    zk_iter = z[k]*ones(s)

    # Creamos un unico vector con todas las variables
    uk_iter = block([yk_iter.flatten(order='F'), zk_iter.flatten(order='F')])
    
    # Inciamos el metodo de Newton
    while (d >= tol) and (cont < itermax):
      F = zeros((s,n))
      G = zeros(s)
      Fy = zeros(s, dtype='O')
      Fz = zeros(s, dtype='O')
      Gy = zeros(s, dtype='O')
      Gz = zeros(s)
      # Evaluamos f, g, fy, fz, gy, gz       
      for i in range(0,s):
          F[i] = f(yk_iter[i],zk_iter[i])
          G[i] = g(yk_iter[i],zk_iter[i])
          Fy[i] = fy(yk_iter[i],zk_iter[i])
          Fz[i] = fz(yk_iter[i],zk_iter[i])
          Gy[i] = gy(yk_iter[i],zk_iter[i])
          Gz[i] = gz(yk_iter[i],zk_iter[i])
      
      # Evaluacion de cada etapa
      Fy = unir_matrices(Fy)
      Fz = unir_vector(Fz)
      Gy = unir_vector(Gy)
      
      # Creamos y rellenamos la matriz Jacobiana 
      DH = zeros((ns + ms, ns + ms))       
      for i in range(0,n):
        for j in range(0,n):
          if i == j:
            DH[i*s:(i+1)*s,j*s:(j+1)*s] = I + E*Fy[i,j]
          else:
            DH[i*s:(i+1)*s,j*s:(j+1)*s] = E*Fy[i,j]
        for j in range(0,m):
          DH[j*s + s*n:(j+1)*s + s*n , i*s :(i+1)*s] = diag(Gy[j])
          
      for j in range(0,m):
        for i in range(0,n):
          DH[i*s:(i+1)*s, j*s + s*n:(j+1)*s + s*n] = E*Fz[i]
        for i in range(0,m):
          DH[j*s + s*n:(j+1)*s + s*n,i*s + s*n:(i+1)*s + s*n] = diag(Gz)
          
      # Creamos H
      H = zeros(ns+ms)
      for i in range(n):
          H[i*s : (i+1)*s] = yk_iter[:,i] - y[i,k]*ones(s) + E@F[:,i]
          
      H[-ms:] = G
      
      uk_new = uk_iter - solve(DH,H)
      d = max(abs(uk_iter - uk_new))
      uk_iter = uk_new
      yk_iter = uk_iter[:ns].reshape((n,s)).transpose()
      zk_iter = uk_iter[-ms:]

      cont += 1

    if cont == itermax:
        print('El metodo no va bien (y_{k+1})')
    else:          
        # Calculamos y_{k+1}
        for i in range(0,n):
            y[i,k+1] = y[i,k] + h*b@F[:,i]
            
        cont = 0
        d = 1.0 + tol
        zk_iter = z[k]
        # Método de Newton para z_{k+1}
        while(d >= tol) and (cont < itermax):  
            znew =  zk_iter - g(y[:,k+1],zk_iter)/gz(y[:,k+1],zk_iter) 
            d = abs(znew - zk_iter)
            zk_iter = znew
            cont += 1
        if cont == itermax:
            print('El metodo no va bien (z_{k+1})')
        else:
            # Calculamos z_{k+1}
            z[k+1] = zk_iter
             
  return (t,y,z)

