import numpy as np 
import scipy.signal as sp 
import pylab 
def build_H(omega,alpha):
    poly1 = np.poly1d([1,omega])
    poly2 = np.poly1d([1,(2*omega),((omega*omega)+(alpha*alpha))])
    return poly1,poly2

def display(i,x,y,xlabel='t',ylabel='x'):
    ''' Function to plot graphs '''
    pylab.figure(i)
    pylab.plot(x,y,'-r')
    pylab.xlabel(r'{}'.format(xlabel),fontsize=15)
    pylab.ylabel(r'{}'.format(ylabel),fontsize=15)
    pylab.show()

def build_f(omega,alpha,t):
    x = np.cos(alpha*t)
    x2 = np.multiply(np.exp(-omega*t),np.heaviside(t,0.5))
    return np.multiply(x,x2)

# Spring system
num,den = build_H(0.5,1.5)
den = np.polymul([1,0,2.25],den)
H1 = sp.lti(num,den)
t = np.linspace(0,50,1000)
sol1 = sp.impulse(H1,T=t)
display(0,sol1[0],sol1[1],'time','x')

# Spring system with different decay
num,den = build_H(0.05,1.5)
den = np.polymul([1,0,2.25],den)
H2 = sp.lti(num,den)
sol2 = sp.impulse(H2,T=t)
display(0,sol2[0],sol2[1],'time','x')

# LTI response over different frequencies of applied force
i=0
for alpha in np.arange(1.4,1.6,0.05):
    H = sp.lti([1],[1,0,2.25])
    t1 = np.linspace(0,100,1000)
    f,x,_ = sp.lsim(H,build_f(0.05,alpha,t1),t1)
    i+=1
    pylab.subplot(3,2,i)
    pylab.plot(t1,x,'-r')

# Coupled spring system
t = np.linspace(0,20,1000)
H_x = sp.lti(np.poly1d([1,0,2]),np.poly1d([1,0,3,0]))
sol_x = sp.impulse(H_x,T=t)
H_y = sp.lti(np.poly1d([2]),np.poly1d([1,0,3,0]))
sol_y = sp.impulse(H_y,T=t)
pylab.figure(5)
pylab.plot(sol_x[0],sol_x[1])
pylab.plot(sol_y[0],sol_y[1])
pylab.xlabel('time')
pylab.ylabel('x/y')
pylab.show()
# Two port network
H = sp.lti(np.poly1d([1000000]),np.poly1d([0.000001,100,1000000]))
w,S,phi=H.bode()
pylab.subplot(2,1,1)
pylab.semilogx(w,S)
pylab.ylabel(r'$|H(s)|$')
pylab.subplot(2,1,2)
pylab.semilogx(w,phi)
pylab.ylabel(r'$\angle(H(s))$')
pylab.show()

t = np.linspace(0,30*0.000001,1000)
vi = np.multiply(np.cos(1000*t)-np.cos(1000000*t),np.heaviside(t,0.5))
_,y,svec = sp.lsim(H,vi,t)
display(5,t,y,'t',r'$v_{o}(t)$')
t = np.linspace(0,10*0.001,100000)
vi = np.multiply(np.cos(1000*t)-np.cos(1000000*t),np.heaviside(t,0.5))
_,y2,svec = sp.lsim(H,vi,t)
display(6,t,y2,'t',r'$v_{o}(t)$')