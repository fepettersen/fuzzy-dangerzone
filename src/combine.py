import numpy as np

class Combine:
	"""
	docstring for Combine:
	This is presently only implemented in 1D
	"""
	def __init__(self,DT,dx,dy,x0=0,x1=1,hc=50,tau=50):
		from diff import Diff
		from walk import Walk

		self.dx = dx
		steplength = np.sqrt(2*DT/tau)
		self.pdesolver = Diff(DT,dx,dy,x0,x1)
		self.limit = 0
		self.length = 8
		self.walksolver = Walk(steplength,1.0/(self.length-1))
		self.a = DT/(2*dx**2)
		self.Hc = hc
		self.tau = tau
		self.C = np.zeros(self.length)
		self.left = self.right = 0

	def SetInitialCondition(self,u0,D):
		dx = 1.0/(self.length+1)
		self.up = u0
		self.pdesolver.Assemble(D,self.a,0)
		self.u = np.zeros(len(self.up))
		self.ConvertToWalkers(u0)
		for i in xrange(self.length):
			for j in xrange(int(self.C[i])):
				pos = i*dx
				self.walksolver.AddWalker(pos)

	def CalculateFluxes(self):
		self.pc = (self.u[self.limit+4]-self.u[self.limit+3])/self.pdesolver.dx
		self.pp = float(self.left - self.right)/(2*self.tau)
		print self.pp
		self.left = self.right = 0

	def Solve(self):
		tau = self.tau
		limit = self.limit
		self.CalculateFluxes()
		
		self.up[limit+1] += self.pp
		self.u[limit:] = self.pdesolver.Solve(self.up[limit:])

		interval = self.pc/tau
		walkers_per_time = 1
		self.ConvertToWalkers(self.u[:self.length])
		for i in xrange(tau):
			if i == interval:
				for j in xrange(walkers_per_time):
					self.walksolver.AddWalker((self.length-1)*self.pdesolver.dx)
			self.walksolver.Solve()
			self.left = self.walksolver.left
			self.right = self.walksolver.right
			# print "-------------------\n nwalers = ",len(self.walksolver.walkers)
		self.ConvertFromWalkers()

		
		self.up = self.u.copy()

	def ConvertToWalkers(self,u):
		for i in xrange(self.length):
			tmp = u[i]*self.Hc
			self.C[i] = int(tmp)

	def ConvertFromWalkers(self):
		dx = 1.0/(self.length)
		self.C = np.zeros(self.length)
		for walker in self.walksolver.walkers:
			indx = int(walker.r/dx)
			# print "r = ",walker.r
			try:
				self.C[indx] += 1
			except IndexError:
				print "IndexError"
				print "index, pos = ",indx,",",walker.r
		self.C[:] /= self.Hc
		# print self.C[:], self.up[:self.length]
		self.up[:self.length] = self.C[:]
		# print self.u[:self.length]
