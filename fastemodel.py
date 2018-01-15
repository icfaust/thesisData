import scipy
import scipy.integrate
import matplotlib.pyplot as plt
import eqtools
from matplotlib import rc
import matplotlib.gridspec as mplgs
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import pickle
from scipy.special import iv,kv
from matplotlib.colors import LogNorm

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('font',size=22)#14)#size=14

tokamaks = ['C-Mod','EAST','ARIES-RS','ITER']
zeta = [56.8,5.6,22.6,28.2]

def loss(r,zeta):
    rp = zeta*r
    B2 = -jn(0,zeta)/yn(0,zeta)
    A1 = 1 + B2*(yn(0,rp)/jn(0,rp))

    f1 = (A1*rp*jn(1,rp))+0


def prof(r, zeta, rin= scipy.mgrid[0:1:1e-4]):
    
    rp = zeta*r
    B2 = -iv(0,zeta)/kv(0,zeta)
    temp = kv(0,rp)/iv(0,rp)
    A1 = (1 + B2*temp)

    output = scipy.zeros(rin.shape)
    idx = rin < r
    output[idx] = A1*(iv(0,zeta*rin[idx]))
    output[~idx] = iv(0,zeta*rin[~idx]) + B2*kv(0,zeta*rin[~idx])

    return abs(output),rin


def test2(r, zeta):
    n,r = prof(r, zeta)
    totes = scipy.integrate.trapz(n*r,r)
    dndr = (n[-1]-n[-2])/(r[-1]-r[-2])
    print(1./(1-(zeta**2)*(totes/dndr)))
    


def tester(r, zeta):
    rin = zeta*r
    A1 = kv(0,rin)/iv(0,rin)
    B2 = iv(0,zeta)/kv(0,zeta)

    temp = (iv(1,zeta)-B2*(A1*r*iv(1,rin)+kv(1,zeta)-r*kv(1,rin)))/(iv(1,zeta)-B2*kv(1,zeta))

    return 1./(1-temp)



def contourout(rho=scipy.linspace(.5,1,2e2), zeta=scipy.linspace(0,60,5e2), levels=[-5,0], pts=11):

    levelin = scipy.logspace(levels[0],levels[1],pts)


    j,k = scipy.meshgrid(rho,zeta)


    out = tester(j,k)

    plt.contour(rho,zeta,out,levels=levelin,norm=LogNorm(),linewidths=2)
    temp = plt.colorbar()
    temp.set_label(r'$f_L$',rotation=90)
    plt.xlabel(r'$\rho$ of deposition')
    plt.ylabel(r'$\zeta$')

    plt.show()
