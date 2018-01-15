import scipy
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy.interpolate

rc('text',usetex=True)
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('font',size=22)#14)#size=14


def plot():
    data = scipy.loadtxt('impact.txt',dtype=float)
    data.T[3][data.T[3] == 0] = scipy.nan
    
    plt.semilogx(data.T[0],1/data.T[1])
    plt.semilogx(data.T[0],1/data.T[2])
    plt.semilogx(data.T[0],1/data.T[3],'.')
    plt.semilogx(data.T[0],1/data.T[4])
    plt.semilogx(data.T[0],1/data.T[5])
    
    plt.xlabel(r'electron energy [eV]')
    plt.ylabel(r'$\lambda_{mfp}n_e$ [$10^{20}$m$^{-2}$]')
    plt.title(r'$e^{-}$ impact ionization mfp for atomic Hydrogen')
    plt.show()

def plot2():
    data = scipy.loadtxt('total_cross_section_hydrogen_yoon_2008.txt',dtype=float)
    
    plt.loglog(data.T[0],data.T[1],linewidth=2.,label=r'H$_2$ total',color='#0E528A')
    
    data = scipy.loadtxt('impact.txt',dtype=float)

    plt.loglog(data.T[0],data.T[4],linewidth=2.,label=r'H ionization',color='#A02C2C')

    output = None
    for i in xrange(9):
        data = scipy.loadtxt('excitation_n'+str(i+2)+'.txt',dtype=float)
        if output is None:
            output = data.T[1]
            x = data.T[0]
        else:
            inter = scipy.interpolate.interp1d(data.T[0],
                                               data.T[1],
                                               fill_value=0,
                                               bounds_error=False)
            output += inter(x)

    plt.loglog(x,output,linewidth=2.,label=r'H excitation',color='#228B22')
    plt.legend()

    plt.ylim(1e-2,5e2)
    plt.xlabel(r'electron energy [eV]')
    plt.ylabel(r'$\sigma$ [$10^{-20}$m$^{2}$]')
    #plt.title(r'$e^{-}$ impact ionization mfp for atomic Hydrogen')
    plt.show()

def plot3():
    output = None
    for i in xrange(9):
        data = scipy.loadtxt('excitation_n'+str(i+2)+'.txt',dtype=float)
        if output is None:
            output = data.T[1]
            x = data.T[0]
        else:
            inter = scipy.interpolate.interp1d(data.T[0],
                                               data.T[1],
                                               fill_value=0,
                                               bounds_error=False)
            output += inter(x)

    plt.loglog(x,output)
    
    plt.xlabel(r'electron energy [eV]')
    plt.ylabel(r'$\lambda_{mfp}n_e$ [$10^{20}$m$^{-2}$]')
    plt.title(r'$e^{-}$ impact ionization mfp for atomic Hydrogen')
    plt.show()
