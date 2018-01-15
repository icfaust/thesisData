import scipy
import scipy.optimize
import scipy.integrate
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import ticker
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

mes = 9.1094e-31
mi = 3.3435e-27
e0 = 8.8541e-12
ec = 1.6021e-19


def nu(n, T):
    return .795e-15*n*pow(T,-1.5)

def P(n, T, B, omega):
    me = mes*(1-1j*nu(n,T)/omega)
    wpe = n*pow(ec,2)/(e0*me)
    wpi = n*pow(ec,2)/(e0*mi)
    return 1 -wpe*pow(omega,-2) - wpi*pow(omega,-2)

def S(n, T, B, omega):
    me = mes*(1-1j*nu(n,T)/omega)
    wce = ec*B/me
    wci = ec*B/mi
    wpe = n*pow(ec,2)/(e0*me)
    wpi = n*pow(ec/omega,2)/(e0*mi)
    return 1 - wpe/(pow(omega,2)-pow(wce,2)) - wpi/(pow(omega,2)-pow(wci,2))

def D(n, T, B, omega):
    me = mes*(1-1j*nu(n,T)/omega)
    wce = ec*B/me
    wci = ec*B/mi
    wpe = n*pow(ec,2)/(e0*me)
    wpi = n*pow(ec/omega,2)/(e0*mi)
    return - wce*wpe/(pow(omega,2)-pow(wce,2))/omega + wci*wpi/(pow(omega,2)-pow(wci,2))/omega

def As(n, T, B, omega, npar):
    return S(n, T, B, omega)

def Bs(n, T, B, omega, npar):
    return (pow(npar,2) - S(n, T, B, omega))*(P(n,T,B,omega) + S(n,T,B,omega)) + pow(D(n,T,B,omega),2)

def Cs(n, T, B, omega, npar):
    return P(n,T,B,omega)*(pow(pow(npar,2) - S(n,T,B,omega),2) - pow(D(n,T,B,omega),2))

def plot(nmin=-2,nmax=1,tmin=0,tmax=2,num=301,B=5.4,npar=1.8):
    Taxis = scipy.logspace(tmin,tmax,num)*1e-3
    naxis = scipy.logspace(nmin,nmax,num)*1e20
    n,T = scipy.meshgrid(naxis,Taxis)
    om = scipy.pi*2*4.6e9
    #wpe = n*pow(ec,2)/(e0*mes)
    print(nu(n,T))
    Ains = As(n, T, B, om, npar)
    Bins = Bs(n, T, B, om, npar)
    Cins = Cs(n, T, B, om, npar)

    Q = (-Bins + scipy.sqrt(pow(Bins,2) - 4*Ains*Cins))/2/Ains
    output = scipy.sqrt(scipy.sqrt(pow(scipy.real(Q),2)+pow(scipy.imag(Q),2))-scipy.real(Q))/pow(2,.5)
    CS = plt.contour(naxis,Taxis*1e3,om*output/299792458.,locator=ticker.LogLocator(),linewidths=2.,colors=colors['b'])#om*output/299792458.
    plt.clabel(CS, inline=1, fontsize=16,fmt=ticker.LogFormatterMathtext())#

    #CS2 = plt.contour(naxis,Taxis*1e3,pow(2,-1)*abs(scipy.imag(Q)/abs(scipy.sqrt(scipy.real(Q))))*om/299792458.,locator=ticker.LogLocator(),linewidths=2.,colors=colors['g'])#om*output/299792458.
    #plt.clabel(CS2, inline=1, fontsize=16,fmt=ticker.LogFormatterMathtext())#

    wpe = n*pow(ec,2)/(e0*mes)
    wce = ec*B/mes
    wci = ec*B/mi
    output2 = nu(n,T)*pow(wpe,.5)/(2*pow(om,2))*scipy.sqrt(abs((npar**2-1)/(1-wpe/(pow(wce,2)*(npar**2-1)))))
    #CS3 = plt.contour(naxis,Taxis*1e3,output2*om/299792458.,locator=ticker.LogLocator(),linewidths=2.,colors=colors['r'])#om*output/299792458.

    output3 = nu(n,T)*pow(wpe,.5)/(2*pow(om,2))*npar
   # CS2 = plt.contour(naxis,Taxis*1e3,output3*om/299792458.,locator=ticker.LogLocator(),linewidths=2.,colors=colors['g'])#om*output/299792458.
    #plt.clabel(CS2, inline=1, fontsize=16,fmt=ticker.LogFormatterMathtext())#


    #plt.clabel(CS3, inline=1, fontsize=16,fmt=ticker.LogFormatterMathtext())#

    A2 = pow(ec,2)/(e0*mi)/pow(om,2)
    B2 = -2*pow(e0*mes,-.5)*ec/wce*npar
    C2 = npar**2-1

    temp =pow((-B2-pow(B2**2-4*A2*C2,.5))/(2*A2),2)
    plt.fill_between(naxis,Taxis[0]*1e3,Taxis[-1]*1e3,naxis>temp,color='k',alpha=.3,linewidth=2.)

    plt.gca().set_yscale("log") 

    plt.gca().set_xscale("log")

    plt.ylabel(r'Electron Temperature [eV]')
    plt.xlabel(r'$n_e$ [m$^{-3}$]')


    plt.show()

    
