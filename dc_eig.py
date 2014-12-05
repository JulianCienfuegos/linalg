from numpy import array, zeros, ones, outer,  linspace, diag, dot, copy
from numpy import concatenate, finfo, multiply
from scipy.optimize import fsolve
from numpy.linalg import norm, eig
from scipy.linalg import block_diag 

def MakeTridiagonalMatrix(main, offset_one):
	"""This function will make a tridiagonal 2D array (NOT a matrix)
	which has the main array on its main diagonal and the offset_one 
	array on its super and sub diagonals.
	"""
	return diag(main) + diag(offset_one, k = -1) + diag(offset_one, k = 1)
	
def dc_eig(T):
	""" This function is an implementation of Dr. G. Hiary's divide and conquer 
	algorithm.
	"""
	eps = finfo(float).eps
	(n,m) = T.shape
	print n, m
	print n
	if(n==1):
		Q = array([1])
		L = T
		return Q, L
	else:
		m = n/2
		v = zeros(n)
		v[m-1:m+1] = 1
		rho = T[m-1][m]
		R = rho*outer(v,v)
		TT = T - R
		T1 = TT[0:m, 0:m]
		T2 = TT[m:, m:]
		print T1
		print T2
		(Q1, L1) = dc_eig(T1)
		(Q2, L2) = dc_eig(T2)
		
		D = block_diag(L1, L2)
		QQ = block_diag(Q1.T, Q2.T)
		u=dot(QQ,v)

		f = lambda x: 1 + rho*dot((multiply(u, 1/(diag(D)-x))),u)
		Ds = copy(diag(D)) # make a deep copy of the diag
		Ds.sort() # sort now, otherwise the diagonal would be sorted in place
		Ds[:] = Ds[::-1] # reverse the list
		L = zeros((n,n))
		
		if(rho > 0):
			for i in range(0,n-1): # Because Python is crazy, this will only go up to n-2.
				print i
				I = array([Ds[i+1]+eps, Ds[i] - eps])
				L[i+1][i+1] = fsolve(f, I)
			x0 = Ds[1] + eps
			if( f(x0) > 0):
				print 'Error 1 : cannot resolve zero of f'
				return
			x1 = x0
			while (f(x1) < 0):
				x1 = x1 + 0.1			
			L[0][0] = fsolve(f, array([x0, x1]))
		else:
			for i in range(0,n-1):
				print i
				I = array([Ds[i+1]+eps, Ds[i] - eps])
				L[i][i] = fsolve(f, I)
			x0 = Ds[n] - eps
			if(f(x0) > 0):
				print('Error 2: cannot resolve zero of f')
				return
			x1 = x0
			while (f(x1) < 0):
				x1 = x1 - 0.1
			L[n][n] = fsolve(f, array([x1, x0]))
		Qprime = full((n, n), nan)
		for i in range (0, n):
			alpha = L[i][i]
			qi = u/(diag(D) - alpha) # I think python defaults to itemwise division
			qi = qi/norm(qi)
			Qprime[:,i] = qi
			Q = block_diag(Q1, Q2)*Qprime
		return Q,L


main = linspace(-10,10, 5)
off = linspace(-1,6,4)
T = MakeTridiagonalMatrix(main, off)
print T
w,v = eig(T)
print w
(Q,L) = dc_eig(T)
print Q
print T
