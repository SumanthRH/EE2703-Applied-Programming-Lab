import numpy as np 
from pylab import *



# Breaking the square region into 101x101 mesh
x = np.linspace(-20,20,103)
x = x[1:-1] # drop last two points as the interval is unbounded
y = x.copy()
X,Y = meshgrid(x,y)

# Defining the wire section wise
wire1 = np.linspace(-30,-10,61)
x = np.array([0])
wire1 = wire1[:-1]
Wire1,x = meshgrid(wire1,x)

wire2 = np.linspace(10,30,61)
wire2 = wire2[1:]
Wire2,x = meshgrid(wire2,x)

theta = np.linspace(0,pi,81)
circ = c_[10*sin(theta),np.zeros_like(theta),10*cos(theta)] # circ is a 3-d array of points spaced equally on the semicircle
# plotting the wire
plot(x[0],Wire1[0],'r')
plot(x[0],Wire2[0],'b')
plot(circ[:,0],circ[:,2],'g')
xlim([-20,20])
ylim([-30,30])
show()

# r_ is an array of coordinates of points on the wire
r_ = np.zeros((201,3))
zero  = np.zeros_like(wire1)
r_[:60] = stack((zero,zero,wire1),axis=-1)
r_[60:141] = circ
r_[141:201] = stack((zero,zero,wire2),axis=-1)
r_ = r_[:-1]

# obtaining an array of length vectors for each section
dl = np.zeros((200,3))
dl[:60] = 20/60.0*np.array([0,0,1]) 
dl[60:140] = 10*pi/81.0*stack([cos(theta[::-1])[:-1],zeros_like(theta)[:-1],-sin(theta[::-1])[:-1]],axis=-1)
dl[140:] = 20/60.0*np.array([0,0,1])

r = np.zeros((len(x),len(y),3))
r = stack((X,Y,np.zeros_like(X)),axis=-1)

def calc(k,r=r,r_=r_,dl=dl):
    ''' calculates B for kth index in r_'''
    R = r - r_[k]
    #print(R.shape)
    #print(dl[k].shape)
    R3 = sum(square(R),axis=2)**1.5
    #print(R3.shape)
    cross_p = cross(dl[k],R)
    #print(cross_p.shape)
    return np.divide(cross_p,stack((R3,R3,R3),axis=-1))

B = sum(np.array([calc(k) for k in range(len(dl))]),axis=0) # Can't avoid the for loop as the function defined evaluates magnetic field section-wise
print(B.shape) 
modB = sum(square(B),axis=-1)**0.5
logB = log10(modB)
B_ = np.divide(B,stack((modB,modB,modB),axis=-1))
print(logB.shape)
# Contour plot of magnitude in log scale
figure(0)
contour(X,Y,logB)
show()
# quiver plot of nromalized magnetic field
figure(1)
quiver(X,Y,B_[0],B_[1])
show()