import scipy
import scipy.optimize
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.gridspec as mplgs
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import pickle
from scipy.special import jn,yn
import scipy.special

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('font',size=24)#14)#size=14

c = ['b','r','g','m','k']
colors = {'g':'#228B22','r':'#A02C2C','b':'#0E525A','m':'#800080','k':'k'}

E = 511.

def energy(R1, R2, m1=2, m2=3):
    return E*(m2*R1/(m1*R2)-1)

def plot(Ebin=100, rr=[.44,.92], pts=3e3, m1=2, m2=3):
    rtemp = scipy.linspace(rr[0],rr[1],pts)
    Ein = (scipy.mgrid[0:(E*(float(m2)/float(m1)-1)):Ebin])[1:]
    r1,r2 = scipy.meshgrid(rtemp,rtemp)
    lim = float(m1)/float(m1+1)*rr[1]

    idx = scipy.searchsorted(rtemp,[lim])[0]

    print(rtemp[idx])


    loc = [[.63,.8],
           [.74,.8]]


    temp = plt.contour(r1[:,idx:],r2[:,idx:],energy(r1,r2,m1=m1,m2=m2)[:,idx:],levels=Ein,linewidths=2.,colors=colors['b'])
    plt.clabel(temp,inline=1, fontsize=16,fmt = '%1.f',manual=loc)

    Ein = (scipy.mgrid[0:(E*(float(m2+1)/float(m1)-1)):Ebin])[1:]


    loc = [[.63,.88],
           [.72,.88],
           [.805,.88],
           [.895,.88]]

    temp = plt.contour(r1[:,idx:],r2[:,idx:],energy(r1,r2,m1=m1,m2=m2+1)[:,idx:],levels=Ein,linewidths=2.,colors=colors['g'])
    plt.clabel(temp,inline=1, fontsize=16, fmt ='%1.f',manual=loc)

    plt.fill_between(rtemp,rr[0]*scipy.ones((pts,)),rtemp,alpha=.4,facecolor='k')
    #plt.fill_between(rtemp[0:idx],rtemp[0:idx],rr[1]*scipy.ones((idx,)),alpha=.4,facecolor=colors['r'])


    plt.xlim([lim,rr[1]])
    plt.ylabel(r'$R_{emitted}$')
    plt.xlabel(r'$R_{observed}$')
    plt.show()
