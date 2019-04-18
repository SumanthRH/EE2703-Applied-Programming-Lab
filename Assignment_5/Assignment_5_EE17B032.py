import numpy as np 
import pylab
import mpl_toolkits.mplot3d.axes3d as p3
import argparse
import matplotlib.pyplot as plt
# using argparse library for accepting the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--Nx',default=25,type=int,help='Size along the x axis')
parser.add_argument('--Ny',default=25,type=int,help='Size along the y axis')
parser.add_argument('--radius',default=8,type=float,help='Radius of central lead')
parser.add_argument('--Niter',default=1000,type=int,help='Number of iterations to perform')
args = parser.parse_args()

Nx,Ny,radius,Niter = args.Nx,args.Ny,args.radius,args.Niter

def phi_update(phi,oldphi):
    ''' Function to update phi using a finite differentation approximation'''
    phi[1:-1,1:-1] = 0.25*(oldphi[1:-1,0:-2]+oldphi[1:-1,2:]+oldphi[0:-2,1:-1]+oldphi[2:,1:-1]) 
    return phi

def boundary_conditions(phi,inds):
    ''' Function to enforce appropriate boundary conditions''' 
    phi[1:-1,0] = phi[1:-1,1]
    phi[0,1:-1] = phi[1,1:-1]
    phi[1:-1,-1] = phi[1:-1,-2]
    phi[-1,1:-1] = 0
    phi[inds] = 1.0
    return phi

def fit_exp(x,A,B):
    ''' Function to obtain an exponential '''
    return A*np.exp(B*x)

def error_fit(x,y):
    ''' Evaluates the parameters of an exponent. x is the data vector, y is the vector of true values ''' 
    logy=np.log(y)
    xvec=np.zeros((len(x),2))
    xvec[:,0]=x
    xvec[:,1]=1
    B,logA=np.linalg.lstsq(xvec, np.transpose(logy))[0]
    return (np.exp(logA),B)

def max_error(A,B,N):
    ''' Finds an upper bound for the error ''' 
    return -A*(np.exp(B*(N+0.5)))/B

phi = np.zeros((Ny,Nx))
x1,y1 = np.linspace(-(Ny-1)/2,(Ny-1)/2,Ny),np.linspace(-(Nx-1)/2,(Nx-1)/2,Nx)
Y,X = np.meshgrid(y1,x1)
inds = np.where((X**2 + Y**2) <= ((radius*Nx)/100)**2)
phi[inds] = 1.0
error = np.zeros(Niter)

# contour of phi 
pylab.figure(0)
pylab.contourf(Y,X,phi,cmap=pylab.cm.get_cmap("autumn"))
pylab.xlabel(r'x$\rightarrow$',fontsize=15)
pylab.ylabel(r'y$\rightarrow$',fontsize=15)
pylab.title('Potential Configuration')
pylab.show()

# Evaluating the potential 
for i in range(Niter):
    oldphi = phi.copy()
    phi = phi_update(phi,oldphi)
    phi = boundary_conditions(phi,inds)
    error[i] = (abs(phi-oldphi)).max()

# Fitting an exponential to the error data
A,B = error_fit(range(Niter),error) # fit1
A_500,B_500 = error_fit(range(Niter)[500:],error[500:]) # fit2

# Evolution of the error function 
pylab.figure(1)
pylab.plot(range(Niter),error,'-r',markersize=3,label='original')
pylab.xlabel(r'Niter$\rightarrow$')
pylab.ylabel(r'Error$\rightarrow$')
pylab.title('Plot of Error vs number of iterations')
pylab.show()

# Error, fit1 and fit2 in a semilog plot 
pylab.figure(2)
pylab.semilogy(range(Niter)[::50],error[::50],'ro',label='original')
pylab.semilogy(range(Niter)[::50],fit_exp(range(Niter)[::50],A,B),'go',label='fit1')
pylab.semilogy(range(Niter)[::50],fit_exp(range(Niter)[::50],A_500,B_500),'bo',label='fit2')
pylab.legend(loc='upper right')
pylab.xlabel(r'Niter$\rightarrow$')
pylab.ylabel(r'Error$\rightarrow$')
pylab.title('Semilog plot of Error vs number of iterations')
pylab.show()


pylab.figure(3)
pylab.loglog(range(Niter)[::50],error[::50],'ro',markersize=3)
pylab.xlabel(r'Niter$\rightarrow$')
pylab.ylabel(r'Error$\rightarrow$')
pylab.title('Loglog plot of Error vs number of iterations')
pylab.show()

# Upper bound on the error
pylab.figure(4)
pylab.semilogy(range(Niter)[::50],max_error(A,B,np.arange(0,Niter,50)),'ro',markersize=3)
pylab.xlabel(r'Niter$\rightarrow$')
pylab.ylabel(r'Error$\rightarrow$')
pylab.title('Semilog plot of Cumulative Error vs number of iterations')
pylab.show()

# 3-D plot of potential
fig1= pylab.figure(5)
ax=p3.Axes3D(fig1) 
pylab.title('The 3-D surface plot of the potential')
surf = ax.plot_surface(Y, X, phi.T, rstride=1, cstride=1,cmap=pylab.cm.jet)
print(X.shape,Y.shape,phi.shape)
pylab.xlabel(r'x$\rightarrow$',fontsize=15)
pylab.ylabel(r'y$\rightarrow$',fontsize=15)
ax.set_zlabel(r'$\phi\rightarrow$',fontsize=15)
pylab.show()

# Contour plot of the potential
pylab.figure(6)
pylab.contour(Y,X[::-1],phi)
pylab.title("Contour plot of the potential")
pylab.plot(inds[1]-(Nx-1)/2,inds[0]-(Ny-1)/2,'ro')
pylab.ylabel(r"y$\rightarrow$")
pylab.xlabel(r'x$\rightarrow$')
pylab.show()

# Current evaluation
Jx = np.zeros((Ny,Nx))
Jy = np.zeros((Ny,Nx))
Jx[:,1:-1] = 0.5*(phi[:,0:-2]-phi[:,2:])
Jy[1:-1,:] = 0.5*(phi[2:, :]-phi[0:-2,:])
fig,ax = plt.subplots()
fig = ax.quiver(Y[2:-1],X[::-1][2:-1],Jx[2:-1],Jy[2:-1],scale=5)
pylab.plot(inds[1]-(Nx-1)/2,inds[0]-(Ny-1)/2,'ro')
plt.title("Current Density")
plt.xlabel(r"x$\rightarrow$")
plt.ylabel(r'y$\rightarrow$')
plt.show()