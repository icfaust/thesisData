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


def zeros(num,a):
    """ returns num values of alpha_n which satisifies 
    J_0(alpha_n*a)Y_0(alpha_n) - J_0(alpha)*Y_0(alpha_n*a) = 0"""
    lim = scipy.special.jn_zeros(0,num+1)/(1-a)
    output = scipy.zeros((num,))
    
    for i in xrange(num):
        output[i] = scipy.optimize.brentq(_func,lim[i],lim[i+1],args=(a,))
        #start = scipy.optimize.newton(_func2,output[i],args=(a,))
    
    return output

def _func(alpha,a):
    return jn(0,alpha*a)*yn(0,alpha) - jn(0,alpha)*yn(0,alpha*a)

def _func2(alpha,a):
    return jn(0,alpha*a)*yn(1,alpha) - jn(1,alpha)*yn(0,alpha*a) + a*(jn(1,alpha*a)*yn(0,alpha) - jn(0,alpha)*yn(1,alpha*a))


def S0(r, alpha, a):
    return yn(0,alpha*r)-(yn(0,alpha*a)/jn(0,alpha*a))*jn(0,alpha*r)

def genB(func, alpha, a):

    S2 = lambda r,alpha,a: r*pow(S0(r,alpha,a),2)
    S = lambda r,alpha,a: r*func(r)*S0(r,alpha,a)
    int1 = scipy.integrate.quad(S,a,1,args=(alpha,a))[0]
    int2 = scipy.integrate.quad(S2,a,1,args=(alpha,a))[0]
    return int1/int2

def coeffs(alpha, a, Tb):
    An = scipy.zeros(alpha.shape)
    Bn = scipy.zeros(alpha.shape)
    
    for i in xrange(len(alpha)):
        Bn[i] = -Tb*genB(scipy.log,alpha[i],a)/scipy.log(a)
        if jn(0,alpha[i]) == 0:
            An[i] = -yn(0,alpha[i]*a)/jn(0,alpha[i]*a)*Bn[i]
        else:
            An[i] = -yn(0,alpha[i])/jn(0,alpha[i])*Bn[i]

    return An,Bn

def ex(num, a, b, Tb):
    #r = scipy.array([.55,.7,.9])
    r = scipy.linspace(a,1,1e3)
    s = ['-',':','-.']
    #t = scipy.linspace(0,1,1e3)
    t = scipy.array([.1,1.,10.])
    toff = scipy.linspace(-1.,0.,1e3)
    r1,t1 = scipy.meshgrid(r,t)
    output = scipy.zeros((len(t),len(r),num))
    alpha = zeros(num, a)
    An,Bn = coeffs(alpha, a, Tb)
    print((jn(0,alpha[0]*r1)).shape,r1.shape,output.shape)
    for i in xrange(num):
        output[:,:,i] = ((An[i]*jn(0,alpha[i]*r1) + Bn[i]*yn(0,alpha[i]*r1))*scipy.exp(-b*pow(alpha[i],2)*t1))
        
    temp  = scipy.sum(output,axis=2) + Tb*scipy.log(r1)/scipy.log(a)
    temp2 = Tb*scipy.log(r1)/scipy.log(a)
    for i in xrange(len(t)):

        #plt.plot(scipy.concatenate((toff,t)),scipy.concatenate((scipy.zeros(toff.shape),temp.T[i])),colors[c[i]],linewidth=2.)
        #plt.plot(t,temp.T[i],colors[c[i]],linewidth=2.)/temp2[:,i]

        #plt.plot(t,temp.T[i]/(Tb*scipy.log(r1)/scipy.log(a))[:,i],colors[c[i]],linewidth=2.)
        
        plt.plot(r,temp[i],'k',linewidth=2.,linestyle=s[i])
    #plt.imshow(temp)
    #plt.colorbar()
    #plt.show()
    #plt.plot(r,output+Tb*scipy.log(r)/scipy.log(a))
    #plt.show()

def ex2(num, a, b, Tb):
    p = [.55,.7,.9,.3,.5]
    s = ['-',':','-.']
    r = scipy.linspace(0,a,1e3)
    #r = scipy.array([.3])
    #t = scipy.linspace(0,1,1e3)
    t = scipy.array([.1,1.,10.])
    toff = scipy.linspace(-1.,0.,1e3)
    r1,t1 = scipy.meshgrid(r,t)
    output = scipy.zeros((len(t),len(r),num))
    temp = scipy.special.jn_zeros(0,num+1)
    alpha = temp/a
    for i in xrange(num):
        output[:,:,i] = (jn(0,alpha[i]*r1)/(alpha[i]*jn(1,alpha[i]*a)))*scipy.exp(-b*pow(alpha[i],2)*t1)
        

    temp  = Tb*(1-2*scipy.sum(output,axis=2)/a)
    print(temp.shape)
    print(-1*scipy.log(a)*(.75*pow(a,2)-.25-.5*pow(a,2)*scipy.log(a)))
    for i in xrange(len(t)):
        #plt.plot(scipy.concatenate((toff,t)),scipy.concatenate((scipy.zeros(toff.shape),temp.T[i])),colors[c[i+3]],linewidth=2.)
        #plt.plot(t,temp.T[i]/Tb,colors[c[i+3]],linewidth=2.)
        
        plt.plot(r,temp[i],'k',linewidth=2.,linestyle=s[i])

    for i in xrange(len(p)):
        plt.axvline(p[i],color=colors[c[i]],linewidth=2.)

    #plt.plot(scipy.concatenate((toff,t)),scipy.concatenate((scipy.zeros(toff.shape),scipy.ones(t.shape))),'-.k',linewidth=2.)
    #plt.gca().get_xaxis().set_visible(False)
    plt.gca().set_ylim([0.,1.1])
    plt.gca().set_yticks(scipy.linspace(0.,1.,5))
    plt.gca().set_xticks([0,1])
    plt.gca().set_yticks([0])
    plt.ylabel(r'$\Delta T$')#\frac{\Delta T}{\Delta T_\infty}$')
    plt.xlabel(r'r')
    plt.show()



def circ(temp = [.55,.7,.9,.3]):
    cent = (0,0)
    #plt.subplot(111,projection='polar')
    for i in xrange(len(temp)):
        plt.gca().add_artist(plt.Circle(cent,temp[i],color=colors[c[i]],fill=False,linewidth=2.))

    plt.gca().add_artist(plt.Circle(cent,.5,color='k',linewidth=2.,fill=False))
    plt.gca().add_artist(plt.Circle(cent,1,color='k',linewidth=1.,fill=False))
    
    #plt.axes([-1.,1.,-1.,1.])
        #plt.gca().add_artist(b)
    plt.xlim([-1.1,1.1])
    plt.ylim([-1.1,1.1])
    plt.gca().set_aspect('equal')
    plt.axis('off')




#a = .5
#x = scipy.linspace(0,1000,1e4)
#y = scipy.linspace(.01,.99,5e2)
#x1,y1 = scipy.meshgrid(x,y)
#plt.plot(x,_func(x,a))
#plt.plot(scipy.special.jn_zeros(0,100)/(1-a)-zeros(100,a),'.r')

#plt.show()
