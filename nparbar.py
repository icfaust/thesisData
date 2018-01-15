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


vals = scipy.sqrt(511/2/scipy.array([.02,.1,1]))/3


a = scipy.array([[1,40]])
plt.figure(figsize=(9, 1.5))
img = plt.imshow(a, cmap="viridis_r")
plt.gca().set_visible(False)
cax = plt.axes([0.1, 0.2, 0.8, 0.6])
temp = plt.colorbar(orientation="horizontal", cax=cax,ticks=vals,extend='max')#[5,10,15,20,25,30,35,40],extend='max')
temp.set_label(r'$n_\parallel$')
plt.show()
