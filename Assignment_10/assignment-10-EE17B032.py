from pylab import *
import argparse 
from scipy import signal
import csv
import time
parser = argparse.ArgumentParser()
parser.add_argument('--file1',type=str,default='./h.csv',help="Path to the first csv file")
parser.add_argument('--file2',type=str,default='./x1.csv',help="Path to the Zadoff-Chu sequence")
args = parser.parse_args()

h = genfromtxt(args.file1,delimiter=',')
w, H = signal.freqz(h)


figure();subplot(2,1,1);plot(w,20*log10(abs(H)),'b');ylabel(r'$|Y|$')
title('Magnitude and Phase response of the filter')
subplot(2,1,2);plot(w,angle(H),'b');ylabel(r'Phase of $Y$');xlabel(r'$\omega$')

show()

n = linspace(1,2**10,2**10)
x = cos(0.2*pi*n) + cos(0.85*pi*n)
figure()
plot(n,x,'r');xlabel(r'$n\rightarrow$',size=15);ylabel(r'$x\rightarrow$',size=15);title(r'Plot of $\cos\left(0.2\pi t\right)+\cos\left(0.85\pi t\right)$')
xlim([1,100]);show()

# Linear convolution
y =  convolve(x,h)
figure()
plot(range(len(n)+len(h)-1),y,'r');xlabel(r'$n\rightarrow$',size=15);ylabel(r'$y\rightarrow$',size=15);title(r'Plot of $y=x*h$')
xlim([1,100]);show()

# Using the DFTs
x_ = np.concatenate((x,zeros(len(h)-1)))
y1= np.fft.ifft(np.fft.fft(x_) * np.fft.fft( np.concatenate( (h,np.zeros(len(x_)-len(h))) )))
figure()
plot(range(len(y1)),y1,'r');xlabel(r'$n\rightarrow$',size=15);ylabel(r'$y\rightarrow$',size=15);title(r'Plot of $y=x*h$ obtained using DFT ')
xlim([1,100]);show()


def circular_conv(x,h):
    P = len(h)
    n_ = int(ceil(log2(P)))
    h_ = np.concatenate((h,np.zeros(int(2**n_)-P)))
    P = len(h_)
    n1 = int(ceil(len(x)/2**n_))
    x_ = np.concatenate((x,np.zeros(n1*(int(2**n_))-len(x))))
    y = np.zeros(len(x_)+len(h_)-1)
    for i in range(n1):
        temp = np.concatenate((x_[i*P:(i+1)*P],np.zeros(P-1)))
        y[i*P:(i+1)*P+P-1] += np.fft.ifft(np.fft.fft(temp) * np.fft.fft( np.concatenate( (h_,np.zeros(len(temp)-len(h_))) ))).real
    return y


y2 = circular_conv(x,h)
figure()
plot(range(len(y2)),y2,'r');xlabel(r'$n\rightarrow$',size=15);ylabel(r'$y\rightarrow$',size=15);title(r'Plot of $y=x*h$ obained from cicular convolution')
xlim([0,100]);show()


lines = []
with open(args.file2,'r') as file2:
    csvreader = csv.reader(file2)
    for row in csvreader:
        lines.append(row)
lines2 = []
for line in lines:
    line = list(line[0])
    try :
        line[line.index('i')]='j'
        lines2.append(line)
    except ValueError:
        lines2.append(line)
        continue
x = [complex(''.join(line)) for line in lines2]
X = np.fft.fft(x)
x2 = np.roll(x,5)
cor = np.fft.ifftshift(np.correlate(x2,x,'full'))
print(len(cor))
figure()
stem(linspace(0,len(cor)-1,len(cor)),abs(cor),'b');xlabel(r'$t\rightarrow$',size=15);ylabel(r'Correlation$\rightarrow$',size=15);title(r'Plot of the auto-correlation')
xlim([1,20]);show()
