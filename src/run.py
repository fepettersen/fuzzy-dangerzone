"""runscript for new test stuff"""
import numpy as np, matplotlib.pyplot as mpl, matplotlib.animation as animation
from combine import Combine

def initial(x,x0=0.5,sigma=0.2):
	return np.exp(-(x-x0)*(x-x0)/(2*sigma*sigma))

N = 32
T = 20
limit = 3

start = 0
stop = 1
x = np.linspace(start,stop,N)
dx = x[1]-x[0]
dy = 0
dt = dx**2
u0 = initial(x)
D = np.ones(N)

im = []
fig = mpl.figure()

stuff = Combine(dt,dx,dy,limit)
stuff.SetInitialCondition(u0,D[limit:])

for i in xrange(T):
	stuff.Solve()
	im.append(mpl.plot(x,stuff.u,'b-'))
	# im.append(mpl.plot(stuff.C,'b-'))
ani = animation.ArtistAnimation(fig,im)
mpl.show()