from numpy import *
import Gnuplot

def condicaoContorno(f,L,nx):
 f[0] = 0.0
 f[nx-1] = 1.0
 return f

def eulerExplicito(nx,dt,k,dx,f):
 for i in range(1,nx-1): # loop sobre a malha
  a = +1
  b = -2
  c = +1
  cte=(dt*k)/(dx*dx)
  varX = ( a*f[i-1] + b*f[i] + c*f[i+1] )
  fnew[i]=f[i] + cte * ( varX )
 return fnew

def residuo(f,fnew,nx):
 residuo = 0.0
 for i in range(0,nx):
  residuo = residuo + abs( f[i]-fnew[i] )
 return residuo

def atualizaF(nx,fnew):
 for i in range(1,nx-1):
  f[i]=fnew[i];
 return f

# Execucao do programa

nx  = 6     # numero de pontos em x
L   = 1.0   # comprimento total
dx  = L/nx  # intervalo dx
dt  = 0.01  # intervalo de tempo
k   = 1
erro = 5
EPS = 1.0E-5 # precisao

f    = zeros( (nx) )
fnew = zeros( (nx) )

gp = Gnuplot.Gnuplot(persist = 1) 
f = condicaoContorno(f,L,nx)
for i in range(1,100):
 fnew = eulerExplicito(nx,dt,k,dx,f)
 erro = residuo(f,fnew,nx)
 f    = atualizaF(nx,fnew)
 gp.plot( f )

#f.tofile("heat1d.txt", sep='\n', format = "%e") 


