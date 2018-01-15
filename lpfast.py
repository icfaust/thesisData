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

mime = 3670.473

def plot(ft=[-2,-1e3,-1e4],fnmin=-7,fnmax=-2,sigf=0,num=1e3):

    fn = scipy.logspace(-5,-1,num)

    f = ft[0]*(1+fn)/(ft[0]+fn)

    plt.loglog(fn,f-1,label=r'ft = $-10^3$')

    f = ft[1]*(1+fn)/(ft[1]+fn)

    plt.loglog(fn,f-1,label=r'ft = $-10^4$')

    f = ft[2]*(1+fn)/(ft[2]+fn)

    plt.loglog(fn,f-1,label=r'ft = $-10^5$')

    plt.xlabel(r'$n_{ef}/n_{es}$')
    plt.ylabel(r'$f-1$')
    plt.legend()
    plt.show()
