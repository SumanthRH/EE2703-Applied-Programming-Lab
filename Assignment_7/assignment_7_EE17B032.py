import numpy as np 
from sympy import *
import scipy.signal as sp 
init_session
import pylab
def display(i,x,y,xlabel='t',ylabel='x'):
    ''' Function to plot graphs '''
    pylab.figure(i)
    pylab.plot(x,y,'-r',label=r'$V_{o}$')
    pylab.xlabel(xlabel,fontsize=15)
    pylab.ylabel(ylabel,fontsize=15)
    pylab.legend(loc ='upper right')
    pylab.grid(True)
    pylab.show()

def sympy_to_lti(xpr, s=symbols('s')):
    """ Convert Sympy transfer function polynomial to Scipy LTI """
    num, den = simplify(xpr).as_numer_denom()  # returns the expressions
    p_num_den = poly(num, s), poly(den, s)
    c_num_den = [expand(p).all_coeffs() for p in p_num_den]  # returns the coefficients
    l_num, l_den = [lambdify((), c)() for c in c_num_den]  # convert to floats
    return sp.lti(l_num, l_den)


def lowpass(R1,R2,C1,C2,G,Vi):
    ''' Lowpass filter '''
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
    b=Matrix([0,0,0,-Vi/R1])
    V = A.inv()*b
    return (A,b,V)

# Obtaining the step response
s = symbols('s')
A,b,V=lowpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
H = sympy_to_lti(Vo)
t = np.linspace(0,0.001,1000)
Vo = sp.step(H,T=t)
display(0,Vo[0],Vo[1],r't$\rightarrow$',r'$V_{o}\rightarrow$')

# Obtaining the response for mixed frequency sinusoid
t = np.linspace(0,0.01,100000)
Vi = np.multiply((np.sin(2000*np.pi*t)+np.cos(2000000*np.pi*t)),np.heaviside(t,0.5))
Vo = sp.lsim(H,Vi,T=t)
pylab.figure(1)
pylab.plot(Vo[0],Vi,label=r'$V_{in}$')
display(1,Vo[0],Vo[1],r't$\rightarrow$',r'$V\rightarrow$')

def highpass(R1,R3,C1,C2,G,Vi):
    ''' High pass filter '''
    s=symbols('s')
    A=Matrix([[0,0,1,-1/G],[-1/(1+1/(s*R3*C2)),1,0,0],[0,-G,G,1],[-s*C1-s*C2-1/R1,s*C2,0,1/R1]])
    b=Matrix([0,0,0,-Vi*s*C1])
    V = A.inv()*b
    return (A,b,V)

# Magnitude response of the high pass filter
A,b,V = highpass(10000,10000,1e-9,1e-9,1.586,1)
Vo = V[3]
H = sympy_to_lti(Vo)
w=pylab.logspace(0,8,801)
ss=1j*w
hf=lambdify(s,Vo,'numpy')
v=hf(ss)
pylab.loglog(w,abs(v),lw=2)
pylab.xlabel(r'$w\rightarrow$')
pylab.ylabel(r'$|H(jw)|\rightarrow$')
pylab.grid(True)
pylab.show()

# response for damped sinusoids
t = np.linspace(0,10,1000)
Vi = np.multiply(np.multiply(np.exp(-0.5*t),np.sin(2*np.pi*t)),np.heaviside(t,0.5))
Vo = sp.lsim(H,Vi,T=t)
pylab.figure(2)
pylab.plot(Vo[0],Vi,label=r'$V_{in}$')
display(2,Vo[0],Vo[1],r't$\rightarrow$',r'$V\rightarrow$')

t = np.linspace(0,0.0001,10000)
Vi = np.multiply(np.multiply(np.exp(-0.5*t),np.sin(2*np.pi*200000*t)),np.heaviside(t,0.5))
Vo = sp.lsim(H,Vi,T=t)
pylab.figure(2)
pylab.plot(Vo[0],Vi,label=r'$V_{in}$')
display(2,Vo[0],Vo[1],r't$\rightarrow$',r'$V\rightarrow$')

# Step response
t = np.linspace(0,0.001,1000)
Vo = sp.step(H,T=t)
display(0,Vo[0],Vo[1],r't$\rightarrow$',r'$V_{o}\rightarrow$')









