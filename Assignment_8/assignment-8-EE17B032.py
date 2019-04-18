from pylab import *
title_dict = {'cos': r"Spectrum of $cos^3(t)$",'sin': r"Spectrum of $sin^3(t)$",'fm': r"Spectrum of $cos(20t+5cos(t))$",
                'gauss': r"Spectrum of $\exp(-t^2/2)$"}
func_dict = {'cos' : lambda x : cos(x)**3,'sin' : lambda x : sin(x)**3,'fm' : lambda x : cos(20*x+5*cos(x)),
               'gauss' : lambda x : exp(-x**2/2) }
def find_dft(func,N=512,r=4*pi,phase_lim=1e-3,xlimit=40,w_lim=40,plot_=True):
    '''Function to find and plot the dft for the inputs given by the dictionaries '''
    t = linspace(-r,r,N+1)
    t = t[:-1]
    # Finding the dft
    y= func_dict[func](t)
    Y = fftshift(fft(y))/N
    w = linspace(-w_lim,w_lim,N+1)
    w = w[:-1]
    if func == 'gauss' : 
        Y = fftshift(abs(fft(y)))/N
        # Normalizing for the case of gaussian
        Y = Y*sqrt(2*pi)/max(Y)
        Y_ = exp(-w**2/2)*sqrt(2*pi)
        print("max error is {}".format(abs(Y-Y_).max()))
    if plot_ == True :
        figure()
        subplot(2,1,1)
        plot(w,abs(Y),lw=2)
        xlim([-xlimit,xlimit])
        ylabel(r"$|y|$", size = 16)
        ttl = title_dict[func]
        title(ttl)
        grid(True)
        subplot(2,1,2)
        ii = where(abs(Y)>phase_lim)
        plot(w[ii], angle(Y[ii]), 'go', lw=2)
        xlim([-xlimit,xlimit])
        ylabel(r"Phase of $Y$", size=16)
        xlabel(r"$k$", size=16)
        grid(True)
        show()
    return Y,w

Y,w = find_dft('cos',xlimit=15)
Y,w = find_dft('sin',xlimit=15)
Y,w = find_dft('fm',xlimit=40)
Y1,w = find_dft('gauss',N=512,r=8*pi,plot_=True,w_lim=32,xlimit=10)
