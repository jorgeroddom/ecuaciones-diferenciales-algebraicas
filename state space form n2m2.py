from pylab import *

def state_space_form(f,fy,fz,g,gy,gz,t0,tf,N,y0,z0,tol,itermax,A,b):
  # Creamos el paso de malla
  h = (tf - t0)/float(N)
  
  # Dimension de y y de z
  n = len(y0)
  m = len(z0)

  # Definimos las matrices de la solucion 
  # numerica e incluimos las condiciones iniciales
  y = zeros([n,N+1])
  z = zeros([m,N+1])
  y[:,0] = y0
  z[:,0] = z0

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
  
  for k in range(0,N):
    cont = 0
    d = 1.0 + tol
    # Creamos los y_k^(1),...,y_k^(s),
    # z_k^(1),...,z_k^(s) iniciales
    yk_iter = y[:,k]*ones((s,n))
    zk_iter = z[:,k]*ones((s,m))

    # Creamos un unico vector con todas las variables
    uk_iter = block([yk_iter.flatten(order='F'), zk_iter.flatten(order='F')])
    
    # Inciamos el metodo de Newton
    while (d >= tol) and (cont < itermax):
      F = zeros((s,n))
      G = zeros((s,m))
      Fy = zeros(s, dtype='O')
      Fz = zeros(s, dtype='O')
      Gy = zeros(s, dtype='O')
      Gz = zeros(s, dtype='O')
      # Evaluamos f, g, fy, fz, gy, gz       
      for i in range(0,s):
          F[i] = f(yk_iter[i],zk_iter[i])
          G[i] = g(yk_iter[i],zk_iter[i])
          Fy[i] = fy(yk_iter[i],zk_iter[i])
          Fz[i] = fz(yk_iter[i],zk_iter[i])
          Gy[i] = gy(yk_iter[i],zk_iter[i])
          Gz[i] = gz(yk_iter[i],zk_iter[i])
      
      # Evaluacion de cada etapa
      Fy2 = unir_matrices(Fy)
      Fz2 = unir_matrices(Fz)
      Gz2 = unir_matrices(Gz)
      Gy2 = unir_matrices(Gy)
      
      
      # Creamos y rellenamos la matriz Jacobiana 
      DH = zeros((ns + ms, ns + ms))       
      for i in range(0,n):
        for j in range(0,n):
          if i == j:
            DH[i*s:(i+1)*s,j*s:(j+1)*s] = I + E*Fy2[i,j]
          else:
            DH[i*s:(i+1)*s,j*s:(j+1)*s] = E*Fy2[i,j]
        for j in range(0,m):
          DH[j*s + s*n:(j+1)*s + s*n , i*s :(i+1)*s] = diag(Gy2[j,i])
          
      for j in range(0,m):
        for i in range(0,n):
          DH[i*s:(i+1)*s, j*s + s*n:(j+1)*s + s*n] = E*Fz2[i,j]
        for i in range(0,m):
          DH[j*s + s*n:(j+1)*s + s*n,i*s + s*n:(i+1)*s + s*n] = diag(Gz2[j,i])
          
      # Creamos H
      H = zeros(ns+ms)
      for i in range(n):
          H[i*s : (i+1)*s] = yk_iter[:,i] - y[i,k]*ones(s) + E@F[:,i]
          
      H[-ms:] = G.flatten(order='F')
      
      uk_new = uk_iter - solve(DH,H)
      d = max(abs(uk_iter - uk_new))
      uk_iter = uk_new
      yk_iter = uk_iter[:ns].reshape((n,s)).transpose()
      zk_iter = uk_iter[-ms:].reshape((m,s)).transpose()

      cont += 1

    if cont == itermax:
        print('El metodo no va bien (y_{k+1})')
    else:          
        # Calculamos y_{k+1}
        for i in range(0,n):
            y[i,k+1] = y[i,k] + h*b@F[:,i]
            
        cont = 0
        d = 1.0 + tol
        zk_iter = z[:,k]
        # Método de Newton para z_{k+1}
        while(d >= tol) and (cont < itermax):  
            znew = zk_iter -  solve(gz(y[:,k+1],zk_iter), g(y[:,k+1],zk_iter))
            d = max(abs(znew - zk_iter))
            zk_iter = copy(znew)
            cont += 1
        if cont == itermax:
            print('El metodo no va bien (z_{k+1})')
        else:
            # Calculamos z_{k+1}
            z[:,k+1] = zk_iter
             
  return (t,y,z)

