import scipy
import scipy.linalg
import scipy.interpolate
import scipy.io
import TRIPPy
import TRIPPy.plot.mayaplot as mplt
import TRIPPy.BPLY
import mayavi.mlab as mlab
import eqtools
import matplotlib.pyplot as plt

xlim = [.4,1.]
ylim = [-.2,.2]
zlim = [-.5,.5]

x1=201
y1=71
z1=201


k = eqtools.EqdskReader(gfile='/home/ian/python/g1120824019.01400')
pts1=101
xgrid = scipy.mgrid[xlim[0]:xlim[1]:complex(0,x1)]
ygrid = scipy.mgrid[ylim[0]:ylim[1]:complex(0,y1)]
zgrid = scipy.mgrid[zlim[0]:zlim[1]:complex(0,z1)]
x,y,z = scipy.mgrid[xlim[0]:xlim[-1]:complex(0,x1-1),ylim[0]:ylim[-1]:complex(0,y1-1),zlim[0]:zlim[-1]:complex(0,z1-1)]

output = scipy.zeros((x1-1,y1-1,z1-1))


cmod = TRIPPy.Tokamak(k)
#mplt.plotTokamak(cmod,angle=[-scipy.pi/6,scipy.pi/6],section=6)
n=0

elements = TRIPPy.BPLY.BPLY(cmod)
for i in elements[:-1]:
    print(n)
    n = n+1
    beams = TRIPPy.beam.multiBeam(i.split(5,5),elements[-1].split(5,5))
    cmod.trace(beams)
    output += TRIPPy.beam.volWeightBeam3d(beams,xgrid,ygrid,zgrid)

print(output.shape,xgrid.shape,ygrid.shape,zgrid.shape)
you = mlab.pipeline.scalar_field(x,y,z,output*1e9)
#you.spacing = [1e3*(xlim[1]-xlim[0])/(x1-1),1e3*(ylim[1]-ylim[0])/(y1-1),1e3*(zlim[1]-zlim[0])/(z1-1)]

temp = mlab.pipeline.volume(you)
beams = TRIPPy.beam.multiBeam(elements[:-1],elements[-1])
cmod.trace()
mplt.plotLine(beams,colormap='Blues')
    
