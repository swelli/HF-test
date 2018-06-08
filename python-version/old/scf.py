#!/user/bin/env python3
import numpy as np
from scipy.integrate import nquad
from scipy.optimize import minimize

distCutoff=5
maxIntegralError=1e-2

def integrate(function):
    func = lambda x,y,z: function(np.sqrt(x**2+y**2+z**2))
    return nquad(func, [(-distCutoff,distCutoff),(-distCutoff,distCutoff),(-distCutoff,distCutoff)], opts={"epsabs":maxIntegralError})[0]

# test:
# --------------------------
# 6-31G basis set for helium
theta_1s = lambda r: 0.0237660*np.exp(-38.4216340*r) + 0.1546790*np.exp(-5.7780300*r) + 0.4696300*np.exp(-1.2417740*r)
theta_1s_prime = lambda r: 1.0000000*np.exp(-0.2979640*r)
# inital guess
phi_1s_alpha = lambda r: 1.0*theta_1s(r)+0.0*theta_1s_prime(r)
phi_1s_beta = lambda r: 0.0*theta_1s(r)+1.0*theta_1s_prime(r)

b = [theta_1s, theta_1s_prime]
C = np.array([[1,0],[0,1]])
def H(k,l):
	return
def S(k,l):
	return
def S(k,l,m,n):
	return integral(lambda r: )

def lagrange(C, l_norm):
	C = C.reshape((2,2))
	function = 0
	constraint = 0
	for i in range(2):
		for k in range(2):
			for l in range(2):
				function += C[i,k]*C[i,l]*H(k,l)
		for j in range(2):
			overlap = 0
			for k in range(2):
				for l in range(2):
					overlap += C[i,k]*C[j,l]*S(k,l)
					for m in range(2):
						for n in range(2):
							function += C[i,k]*C[j,l]*C[i,m]*C[j,n]*K(k,l,m,n)
			if i==j:
				constraint += 1 - overlap
			else:
				constraint -= overlap
	return function - l_norm * constraint

# Caluclate S
S_11 = integrate(lambda r: phi_1s_alpha(r)*phi_1s_alpha(r))
S_12 = integrate(lambda r: phi_1s_alpha(r)*phi_1s_beta(r))
S_22 = integrate(lambda r: phi_1s_beta(r)*phi_1s_beta(r))

S = np.array([[S_11,S_12],[S_12,S_22]])
print(S)
