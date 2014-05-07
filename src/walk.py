# -*- coding: utf-8 
# simple random walk program

import numpy as np
import random
from walker import Walker


class Walk:
	"""
	Random walk class at the moment in 1D and 2D only
	"""
	def __init__(self,steplength,dx,d=1):
		self.dx = dx

		self.d = d 	#dimension of the walk
		self.steplength = steplength
		self.walkers = []
		self.nwalkers = 0
		self.left = self.right = 0


		self._x1 = 1 + self.dx/2.0
		self._x0 = 0 - self.dx/2.0

	def AddWalker(self,pos):
		tmp = Walker(pos)
		self.walkers.append(tmp)


	def Solve(self):

		self.left = self.right = 0
		l_limit = 3*self.dx
		r_limit = 4*self.dx
		indices = []
		walkers_leaving_area = []
		counter = 0

		for walker in self.walkers:
			direction = (-1)**random.randint(0,1)
			if walker.r > l_limit and walker.r<r_limit:
				if direction < 0:
					self.left += 1
				else:
					self.right +=1
			walker.r += direction*self.steplength
			self.checkpos(walker)

		self.nwalkers = len(self.walkers)


	
	
	def checkpos(self,walker):
		"""Implements reflecting boundaries"""
		_x0 = self._x0
		_x1 = self._x1
		if walker.r < _x0:
			walker.r += 2*(_x0-walker.r)
		elif walker.r > _x1:
			walker.r += 2*(_x1-walker.r)



