
from pylab import *
import scipy.special as sp
N=101 # no of data points
k=9
# no of sets of data with varying noise
# generate the data points and add noise
t=linspace(0,10,N)
# t vector
y=1.05*sp.jn(2,t)-0.105*t # f(t) vector
Y=meshgrid(y,ones(k),indexing="ij")[0] # make k copies
scl=logspace(-1,-3,k) # noise stdev
n=dot(randn(N,k),diag(scl)) # generate k vectors
yy=Y+n
# add noise to signal
# shadow plot
plot(t,yy)
xlabel(r"$t$",size=20)
ylabel(r"$f(t)+n$",size=20)
title(r"Plot of the data to be fitted")
grid(True)
savetxt("fitting4.dat",c_[t,yy]) # write out matrix to file
show()