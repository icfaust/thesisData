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
rc('font',size=18)#14)#size=14

c = ['b','r','g','m','k','o']
colors = {'g':'#228B22','r':'#A02C2C','b':'#0E528A','m':'#800080','k':'k','o':' #CF5300'}



def test(taue=35e-3,L=50e-3,pts=200,n=50,Wo=1.):
    output = scipy.zeros((pts,))
    t = scipy.linspace(0,2*L,pts)
    pi = scipy.pi


    A = -Wo/(1 - scipy.exp(-L/taue))
    B = -A - Wo/2

    for i in scipy.arange(n)+1:
        b = i*scipy.pi/L
        output += (A*b/(pow(b,2)+pow(1./taue,2))*(1-pow(-1,i)*scipy.exp(-L/taue)) + 0*2*B/i/pi*(1-pow(-1,i)))*scipy.sin(b*t)

    plt.plot(t,output,t,.1+A*scipy.exp(-t/taue)+B)
    plt.show()

def frac(tpte = 1e1, pts = 1e4):
    t = scipy.linspace(0,tpte,pts)
    

    plt.plot(t[1:],1-scipy.exp(-t[1:]/2),linewidth=2,color='#0E525A')


    plt.xlabel(r' $\tau_{mod}/\tau_E$')
    plt.ylabel(r' $\Delta W/ \Delta W_\infty$')
    plt.ylim([0,1])
    plt.show()


def pltex(tpte = 1e1, pts = 1e4, num = 6):

    ax1 = plt.subplot2grid((3,3), (0,0))
    ax2 = plt.subplot2grid((3,3), (1,0))
    ax3 = plt.subplot2grid((3,3), (2, 0)) 
    
    temp = (ax1,ax2,ax3)
    temp2 = [.3,2.,3.]

    t = scipy.linspace(0,tpte,pts)
    ax4 = plt.subplot2grid((3,3), (0, 1),colspan=3,rowspan=3)

    ax4.plot(t[1:],1-scipy.exp(-t[1:]/2),linewidth=2,color='#0E525A')

    for j in xrange(3):
        t,data1,data2 = intime(temp2[j], num, pts)
        data2[-1] = 0. #kludge
        print(t.shape,data1.shape,data2.shape)
        temp[j].plot(t,data1,color=colors[c[1+j]],linewidth=2.)
        temp[j].plot(t,data2,color='k',linestyle=':',linewidth=2.)
        temp[j].set_ylim([0,1.2])
        temp[j].spines['top'].set_visible(False)
        temp[j].spines['right'].set_visible(False)
        ax4.plot(temp2[j],1-scipy.exp(-temp2[j]/2),linestyle='none',marker='.',markersize=18.,color=colors[c[1+j]])
    
    ax3.set_xlabel(r' t [$\tau_E$]')

    ax4.set_xlabel(r' $\tau_{mod}/\tau_E$')
    ax4.set_ylabel(r' $\Delta W/ \Delta W_\infty$')


    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    plt.tight_layout(w_pad=.5,h_pad=0)
    plt.show()


def intime(tpte, num, pts=1e3, w=1):
    # num is the number of taues to do
    tpte = float(tpte)
    t = scipy.linspace(0,tpte/2,pts)
    temp = int(num/tpte)
    output = scipy.zeros((temp,pts))
    output2 = scipy.zeros(output.shape)

    output[0::2] = w*(1-scipy.exp(-t))#/(1-scipy.exp(-tpte/2))
    output[1::2] = w*(scipy.exp(-t)-scipy.exp(-tpte/2))#/(1-scipy.exp(-tpte/2)) #laziness abides
    output2[0::2] = scipy.ones((pts,))


    return scipy.linspace(0,2*num,output.size),output.flatten(),output2.flatten()
