from pylab import *

import mpl_toolkits.mplot3d.axes3d as p3
get_title = {'cos3' : r"Spectrum of $\cos^3(0.86t)$",'cos':r"Spectrum of $\cos{t}$",'chirp':r"Spectrum of $\cos\left(16\left(1.5+\frac{t}{2\pi}\right)t\right)$"}
get_func = {'cos3':lambda x: cos(0.86*x)**3,'cos':cos,'chirp' : lambda t : cos(16*t*(1.5+t/(2*pi)))}
def find_fft_np(func,N=512,t_lim_1=-pi,t_lim_2=pi,windowing=True,x_limit=8,plot_=True):
    '''Function to find the fft of a non periodic function. Returns the fft and the frequency array as : (fft,freqs)
    
    args :
        func :
            function key in the get_func dictionary
        N : 
            number of samples
        t_lim_1,t_lim_2 : 
            range in time domain
        x_limit : 
            frequency limit for the  plot
        plot_ : 
            Boolean to specify plotting of the magnitude and phase of the fft
        windowing :
            Boolean to specify usage of Hamming window
    '''    
    t=linspace(t_lim_1,t_lim_2,N+1);t=t[:-1]
    dt=t[1]-t[0];fmax=1/dt
    wnd = 1
    if windowing == True:
        n=arange(N)
        wnd=fftshift(0.54+0.46*cos(2*pi*n/(N-1)))
    y=get_func[func](t)
    y=y*wnd
    y[0]=0 
    y=fftshift(y)
    Y=fftshift(fft(y))/float(N)
    w=linspace(-pi*fmax,pi*fmax,N+1);w=w[:-1]
    if plot_ == True:
        figure()
        subplot(2,1,1)
        plot(w,abs(Y),'-bo',lw=2)
        xlim([-x_limit,x_limit])
        ylabel(r"$|Y|$",size=16)
        title(get_title[func])
        grid(True)
        subplot(2,1,2)
        plot(w,angle(Y),'ro',lw=2)
        xlim([-x_limit,x_limit])
        ylabel(r"Phase of $Y$",size=16)
        show()
    return Y,w


Y,w = find_fft_np('cos3',t_lim_1 = -4*pi,t_lim_2=4*pi,windowing=False)
Y ,_ = find_fft_np('cos3',t_lim_1 = -4*pi,t_lim_2=4*pi,windowing=True)

t_vec = linspace(-pi,pi,129);t_vec = t_vec[:-1]

def get_cos(t,w,delta):
    return cos((w*t)+delta)

cos1 = get_cos(t_vec,0.5,pi)

def estimate_w(actual_func,t_vec,x_limit=8,pow_=1.7):
    N = len(t_vec)
    dt = t_vec[1]-t_vec[0];fmax = 1/dt
    w = linspace(-pi*fmax,pi*fmax,N+1);w=w[:-1]
    y = actual_func
    n=arange(N)
    wnd=fftshift(0.54+0.46*cos(2*pi*n/(N-1)))
    y = y*wnd
    y[0]=0
    Y = fftshift(fft(fftshift(y)))/float(N)
    delta = angle(Y[::-1][argmax(abs(Y[::-1]))])
    w0 = sum(abs(Y**pow_*w))/sum(abs(Y)**pow_)
    w_inds = where(abs(Y)==abs(Y).max())
    figure()
    subplot(2,1,1)
    plot(w,abs(Y),'-bo',lw=2)
    xlim([-x_limit,x_limit])
    ylabel(r"$|Y|$",size=16)
    title(r"Spectrum of the function")
    grid(True)
    subplot(2,1,2)
    plot(w,angle(Y),'ro',lw=2)
    xlim([-x_limit,x_limit])
    ylabel(r"Phase of $Y$",size=16)
    show()
    return w0,delta

print(" Estimates for w and delta (without noise) are :",estimate_w(cos1,t_vec,pow_=1.7))
cos2 = cos1 +  0.1*randn(128)
print(" Estimates for w and delta (with noise) are :",estimate_w(cos2,t_vec,pow_=2.4))

find_fft_np('chirp',t_lim_1=-pi,t_lim_2=pi,N=1024,plot_=True,x_limit=50,windowing=False)
Y =[]
t = linspace(-pi,pi,1024)
figure()
plot(t,get_func['chirp'](t),'b');xlabel(r'$t\rightarrow$',size=15);ylabel(r'Chirped signal$\rightarrow$',size=15);title(r'Plot of $\cos\left(16\left(1.5+\frac{t}{2\pi}\right)t\right)$')
show()

for i in range(16):
    t_lim_1 = -pi + (2*pi)*i/16
    t_lim_2 = -pi + (2*pi)*(i+1)/16
    Y.append(find_fft_np('chirp',t_lim_1=t_lim_1,t_lim_2=t_lim_2,N=64,plot_=False,windowing=False)[0])

Y = np.array(Y)
t1 = linspace(-pi,pi,16)
t =  linspace(-pi,pi,1025);t=t[:-1]
dt = t[1]-t[0];fmax=1/dt
w = linspace(-pi*fmax,pi*fmax,65);w=w[:-1]
ax=p3.Axes3D(figure())
Y1 = Y.copy()
inds = where(abs(w)>150)
Y1[:,inds] = NaN
t1,w = meshgrid(t1,w) 
surf = ax.plot_surface(t1,w,abs(Y1).T,rstride=1,cstride=1,cmap=cm.jet)
ylabel(r'$\omega\rightarrow$',size=16)
xlabel(r'$t\rightarrow$',size=16)
ax.set_ylim([-150,150])
ax.set_zlabel(r'$|Y|$')
show()