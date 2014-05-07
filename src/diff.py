import numpy as np

class Diff:
	"""docstring for Diff"""
	def __init__(self,DT,dx,dy,x0=0,x1=1):
		self.preconditioned = False
		self.t = 0
		self.dt = DT
		self.dx = dx
		self.m = int(1.0/dx)-1
		self.n = n = 1
		self.x = np.linspace(x0,x1,self.m)

	def Assemble(self,D,a,b,v=0):
		"""Assemble the matrix A which will be constant as long as dt is constant, 
		and make a LU decomposition of A.
		D is the diffusion \"constant\" given as a np vector, a and b are the 
		prefactors a = dt/(2*dx**2) and b = dt/(2*dy**2)"""

		drift = v
		if not np.shape(D):
			import sys
			print "Cannot assemble from scalar diffusion constant. Use D = np.ones(..)*D. \n exiting \n "
			sys.exit(1)

		try:
			m,n = np.shape(D)
		except ValueError:
			m = len(D)
			n = 1
		self.m = m
		self.n = n
		N = m*n
		self.A = np.zeros((N,N))
		if n==1:
			k=0
			self.A[k,k] = 1+2*a*(D[1]+D[0])
			self.A[k,k+1] = -2*a*(D[0]+D[1])
			k+=1
			for i in xrange(1,m-1):
				self.A[k,k] = 1+a*D[i+1]+2*a*D[i] +a*D[i-1]
				self.A[k,k+1] = -a*(D[i+1]+D[i])+drift
				self.A[k,k-1] = -a*(D[i-1]+D[i])+drift
				k+=1
			self.A[k,k] = 1+2*a*(D[m-2]+D[m-1])
			self.A[k,k-1] = -2*a*(D[m-1]+D[m-2])
			return
		k=0
		for i in xrange(m):
			self.A[i,i+n] = -2*a*(D[0,i]+D[1,i]);
			self.A[N-1-i,N-i-n-1] = -2*a*(D[m-1,m-1-i]+D[m-2,m-1-i]);
			for j in xrange(n):
				if j==0:
					self.A[k,k+1] = -2*b*(D[i,j+1]+D[i,j])
					if i==0:
						self.A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
					elif i==m-1:
						self.A[k,k-1] = 0
						self.A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
					else:
						self.A[k,k-1] = 0
						self.A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) + 2*b*D[i,j+1]
				elif j==n-1:
					self.A[k,k-1] = -2*b*(D[i,j-1]+D[i,j])
					if i==0:
						self.A[k,k+1] = 0
						self.A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + 2*b*D[i,j-1]
					elif i==m-1:
						self.A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) +2*b*D[i,j-1]
					else:
						self.A[k,k+1] = 0
						self.A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) +2*b*D[i,j-1]
				elif j!=0:
					self.A[k,k+1] = -b*(D[i,j+1]+D[i,j])
					self.A[k,k-1] = -b*(D[i,j-1]+D[i,j])
					if i==0:
						self.A[k,k] = 1+2*a*D[i+1,j] + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
					elif i==m-1:
						self.A[k,k] = 1+2*a*D[i-1,j] + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
					else:
						self.A[k,k] = 1+a*(D[i+1,j]+D[i-1,j]) + D[i,j]*(2*a+2*b) + b*(D[i,j+1]+D[i,j-1])
				if (k>(m-1) and k<(N-m)):
					self.A[k,k+m] = -a*(D[i+1,j]+D[i,j]);
					self.A[k,k-m] = -a*(D[i-1,j]+D[i,j]);
				k+=1


	def Precondition(self,A):
		self.preconditioned = True
		n = self.n
		m = self.m
		# n = int(np.sqrt(np.shape(A)[0]))
		N = m*n
		self.H = []
		self.D_ = []
		self.Aa = []
		i=0
		a = np.zeros((n,n))
		b = np.zeros((n,n))
		c = np.zeros((n,n))

		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

		self.D_.append(np.linalg.inv(b))
		self.H.append(-1*np.dot(self.D_[-1],c))
		for i in xrange(1,m-1):
			a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
			b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
			c[:,:] = A[i*n:(i+1)*n,(i+1)*n:(i+2)*n]

			self.Aa.append(a.copy())
			self.D_.append(np.linalg.inv(b+np.dot(a,self.H[-1])))
			self.H.append(-1*np.dot(self.D_[-1],c))
			
		i=m-1

		a[:,:] = A[i*n:(i+1)*n,(i-1)*n:i*n]
		self.Aa.append(a.copy())
		b[:,:] = A[i*n:(i+1)*n,i*n:(i+1)*n]
		self.D_.append(np.linalg.inv(b+np.dot(a,self.H[-1])))

	def modified_Block_Tridiag(self,H,D,a,Up):
		n = self.n
		m = self.m
		N = m*n
		g = []
		k = np.zeros(n)
		x = np.zeros(N)


		i = 0
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[0],k))
		for i in xrange(1,m-1):
			gtmp = g[-1]
			k[:] = Up[i*n:(i+1)*n]
			g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))
		i = m-1
		# print "Up[%d:%d] = "%(i*n,(i+1)*n),Up[i*n],",",Up[(i+1)*n]
		gtmp = g[-1]
		k[:] = Up[i*n:(i+1)*n]
		g.append(np.dot(D[i],(k-np.dot(a[i-1],gtmp))))

		x[i*n:(i+1)*n] = g[i]
		for i in xrange(m-2,-1,-1):
			x[i*n:(i+1)*n] = g[i]+np.dot(H[i],x[(i+1)*n:(i+2)*n])

		return x

	def Block_Tridiag(self,A,Up):
		# Solves Ax=b using an efficient block tridiagonal solver
		H = []
		tmp = []
		g = []
		n = int(sqrt(np.shape(A)[0]))
		N = n**2
		a = np.zeros((n,n))
		b = np.zeros((n,n))
		c = np.zeros((n,n))
		k = np.zeros(n)
		i=0
		b[:,:] = A[i*n:n*(i+1),i*n:n*(i+1)]
		c[:,:] = A[i*n:n*(i+1),n*(i+1):n*(i+2)]
		H.append(-1*np.dot(np.linalg.inv(b),c))
		k[:] = Up[0:n]
		g.append(np.dot(np.linalg.inv(b),k))
		for i in xrange(1,n-1):
			a[:,:] = A[i*n:n*(i+1),n*(i-1):i*n]
			b[:,:] = A[i*n:n*(i+1),i*n:n*(i+1)]
			c[:,:] = A[i*n:n*(i+1),n*(i+1):n*(i+2)]
			tmp.append(np.linalg.inv(b+np.dot(a,H[i-1])))
			H.append(-1*np.dot(tmp[i-1],c))
			k[:] = Up[i*n:(i+1)*n]
			g.append(np.dot(tmp[i-1],(k-np.dot(a,g[i-1]))))
		i = n-1
		a[:,:] = A[i*n:n*(i+1),n*(i-1):i*n]
		b[:,:] = A[i*n:n*(i+1),i*n:n*(i+1)]
		k[:] = Up[i*n:(i+1)*n]
		tmp.append(np.linalg.inv(b+np.dot(a,H[i-1])))
		g.append(np.dot(tmp[-1],(k-np.dot(a,g[i-1]))))
		
		x = np.zeros(np.shape(Up))
		x[N-1-n:N-1] = g[-1][:]
		for i in xrange(n-1,0,-1):
			# print x[i*n:(i+1)*n]
			x[(i-1)*n:i*n] = g[i] + np.dot(H[i-1],x[i*n:(i+1)*n])
		return x


	def Tridiag(self, up, a, b,c):
		#Specialized gaussian elimination for tridiagonal matrices (BE1D)
		m = len(up)
		temp = np.zeros(m)
		x = np.zeros(m)
		btemp = b[0]
		x[0] = up[0]/btemp;
		for i in xrange(1,m):
			#Forward substitution
			temp[i] = c[i-1]/btemp
			btemp = b[i]-a[i-1]*temp[i];
			x[i] = (up[i] -a[i-1]*x[i-1])/btemp;

		for i in xrange(m-2,-1,-1):
			#Backward substitution
			x[i] -= temp[i+1]*x[i+1];

		return x

	def Source(self,x,t):
		return 0

	def Solve(self,Up):
		Up += self.dt*self.Source(self.x,self.t)
		if self.preconditioned:
			U = self.modified_Block_Tridiag(self.H,self.D_,self.Aa,Up)
		else:
			U = self.Tridiag(Up,np.diagonal(self.A,-1),np.diagonal(self.A,0),np.diagonal(self.A,1))
		self.t += self.dt
		return U