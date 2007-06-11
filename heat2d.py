from numpy import zeros,double,reshape

def malha( L,nx,ny,X,Y ):
 count=0;
 for i in range( 0,nx ):
  for j in range( 0,ny ):
   X[count] = i*L/(nx-1);
   Y[count] = j*L/(ny-1);
   count=count+1;
 return 0

def condicaoContorno(f,L,nx,ny):
 for j in range(1,nx):
  f[0][j] = j*L/(nx-1)
  fnew[0][j] = j*L/(nx-1)

 for i in range(1,ny):
  f[i][0] = i*L/(ny-1)
  fnew[i][0] = i*L/(ny-1)

 for j in range(1,nx):
  f[ny-1][j] = f[0][nx-1] + j*L/(nx-1) * j*L/(nx-1)
  fnew[ny-1][j] = f[0][nx-1] + j*L/(nx-1) * j*L/(nx-1)
 
 for i in range(1,ny):
  f[i][nx-1] = f[ny-1][0] + i*L/(ny-1) * i*L/(ny-1)
  fnew[i][nx-1] = f[ny-1][0] + i*L/(ny-1) * i*L/(ny-1)
 return f

def eulerExplicito(nx,ny,dt,k,dx,dy,f):
 for i in range( 1,nx-1 ): # loop sobre a malha
  for j in range( 1,ny-1 ):
   a = +1
   b = -2
   c = +1
   cteX=(dt*k)/(dx*dx)
   cteY=(dt*k)/(dy*dy)
   varX = ( a*f[i-1][j] + b*f[i][j] + c*f[i+1][j] )
   varY = ( a*f[i][j-1] + b*f[i][j] + c*f[i][j+1] )
   fnew[i][j]=f[i][j] + cteX * varX + cteY * varY 
 return fnew

def eulerImplicito(nx,ny,dt,k,dx,dy,f):
 K = zeros( (nx,ny),double )
 for i in range( 1,nx-1 ):
  for j in range( 1,ny-1 ):
   cteX=(dt*k)/(dx*dx)
   cteY=(dt*k)/(dy*dy)
   a = -cteX
   b = (1+2*(cteX+cteY))
   c = -cteY
   K[i][j] = a*(fnew[i-1][j]+fnew[i+1][j]) + b*fnew[i][j] + c*(fnew[i][j-1]+fnew[i][j+1])
   print K
 return K

def residuo(f,fnew,nx,ny):
 residuo = 0.0
 for i in range(0,nx):
  for j in range(0,ny):
   residuo = residuo + abs( f[i][j]-fnew[i][j] )
 return residuo

def atualizaF( nx,ny,fnew ):
 for i in range( 1,nx-1 ):
  for j in range( 1,ny-1 ):
   f[i][j]=fnew[i][j];
 return f

def vtkOut(X,Y,f,nx,ny):
 vtkFile = open( 'field.vtk', 'w' )
 vtkFile.write( "# vtk DataFile Version 1.0\n" )
 vtkFile.write( "Simulacao Numerica Heat 2D\n" )
 vtkFile.write( "ASCII\n" )
 vtkFile.write( "DATASET STRUCTURED_GRID\n" )
 vtkFile.write( "DIMENSIONS %d %d %d\n" %(nx,ny,1) )
 vtkFile.write( "POINTS %d float\n" %X.size )

 fOut = zeros( (nx*ny),float )
 fOut = reshape( f,(nx*ny) ) 
 for i in range( 0,nx*ny ):
  vtkFile.write( "%f %f %f\n" %(X[i],Y[i],X[i]*0) )
 
 vtkFile.write( "POINT_DATA %d\n" %X.size )
 vtkFile.write( "SCALARS heat float\n" )
 vtkFile.write( "LOOKUP_TABLE default\n" )
 for i in range( 0,nx*ny ):
  vtkFile.write( "%f\n" %fOut[i] )

 return 0


# Execucao do programa

nx  = 3     # numero de pontos em x
ny  = 3     # numero de pontos em y
L   = 1.0   # comprimento total
dx  = L/nx  # intervalo dx
dy  = L/ny  # intervalo dy
dt  = 0.002  # intervalo de tempo
k   = 1
EPS = 1.0E-5 # precisao

f    = zeros( (nx,ny),double );
fnew = zeros( (nx,ny),double );
X    = zeros( (nx*ny),double );
Y    = zeros( (nx*ny),double );

malha(L,nx,ny,X,Y)
f = condicaoContorno(f,L,nx,ny)
for i in range(1,1000):
 fnew = eulerExplicito(nx,ny,dt,k,dx,dy,f)
 erro = residuo(f,fnew,nx,ny)
 if( erro < EPS ):
  break
 f    = atualizaF(nx,ny,fnew)

print( "Convergiu para uma solucao com erro < 10e-05" )
print( "erro = %f" %erro )
print( "iteracao = %d" %i )

vtkOut(X,Y,f,nx,ny)

