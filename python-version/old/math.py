#!/usr/bin/env python3
import numpy as np
import sympy as sp
import scipy.integrate
from sympy.parsing import sympy_parser
import pdb
distCutoff=5
maxIntegralError=1e-2

class Quadrature:
	# class attributes (= static/shared)
	method = None
	bounds = None
	epsabs = None
	def __init__(self, method=scipy.integrate.nquad, bounds=[(-np.inf,np.inf),(-np.inf,np.inf),(-np.inf,np.inf)], epsabs=1.49e-8):
		self.method = method
		self.bounds = bounds
		self.epsabs = epsabs
	def __call__(self, function):
		return self.method(function, bounds, epsabs=epsabs)

class MathFunction:
	'''
	Class of functions that support symbolic math
	'''
	def __init__(self, expression):
		if isinstance(expression, str):
			self.expression = sympy_parser.parse_expr(expression)
		else:
		# elif isinstance(expression, sympy.core):
			self.expression = expression
		# else:
			# raise TypeError("Input must either be a sympy expression or string")
		# cache
		self.parsed = None
	def symbols(self):
		'''
		Returns all sympy symbols in alphabetically order
		'''
		return tuple(sp.ordered(self.expression.free_symbols))
	def replace(self, *args):
		'''
		Inplace version of sympy's <expr>.replace method
		'''
		self.expression = self.expression.replace(*args)
	def subs(self, *args):
		'''
		Inplace version of sympy's <expr>.subs method
		'''
		self.expression = self.expression.subs(*args)
	def __call__(self, *args):
		# load cache
		if self.parsed is None:
			self.parsed =  sp.lambdify(self.symbols(), self.expression, modules=np)
		return self.parsed(*args)
	def __str__(self):
		return str(self.expression)
	def __pow__(self, other):
		if isinstance(other, MathFunction):
			return MathFunction(self.expression ** other.expression)
		else:
			return MathFunction(self(*args) ** other)
	def __ipow__(self, other):
		# delete cache
		self.parsed = None
		if isinstance(other, MathFunction):
			self.expression **= other.expression
		else:
			self.expression **= other
	def __mul__(self, other):
		if isinstance(other, MathFunction):
			return MathFunction(self.expression * other.expression)
		else:
			return MathFunction(self.expression * other)	
	def __rmul__(self, other):
		return self.__mul__(other)	
	def __imul__(self, other):
		self.parsed = None
		if isinstance(other, MathFunction):
			self.expression *= other.expression
		else:
			self.expression *= other
		return self
	def __add__(self, other):
		if isinstance(other, MathFunction):
			return MathFunction(self.expression + other.expression)
		else:
			return MathFunction(self.expression + other)	
	def __radd__(self, other):
		return self._add__(other)	
	def __iadd__(self, other):
		# delete cache
		self.parsed = None
		if isinstance(other, MathFunction):
			self.expression += other.expression
		else:
			self.expression += other
		return self
	def __sub__(self, other):
		if isinstance(other, MathFunction):
			return MathFunction(self.expression - other.expression)
		else:
			return MathFunction(self.expression - other)	
	def __rsub__(self, other):
		return self.__sub__(other)	
	def __isub__(self, other):
		# delete cache
		self.parsed = None
		if isinstance(other, MathFunction):
			self.expression -= other.expression
		else:
			self.expression -= other
		return self


class Wavefunction(MathFunction):
	'''
	Parent class of all wavefunction
	'''
	quadrature = Quadrature()
	def __init__(self, function):
		super().__init__(function)
	def normalize(self):
		self.expression /= self.quadrature(self**2)
class GaussianPrimitive(MathFunction):
	def __init__(self, d, a):
		super().__init__("d*exp(-a*r**2)")
		self.subs([("d", d),("a", a)])

class GaussianBasis(self, ):

quadrature = Quadrature()
a = Wavefunction("d_a*exp(-a_a*r)")
b = Wavefunction("d_b*exp(-a_b*r)")
# print(a.expression, a.symbols())
# print(b.expression, b.symbols())
for i in range(20*10):
	a+=b
	b+=a
a*=b
a.replace("r", "(x1**2+y1**2+z1**2)**0.5")
print(a.expression, a.symbols())
# b = Wavefunction(lambda x,a=2:x*a)
# print((2*a*a**a(1,a=3))(1))
# print(a(1))
# pdb.set_trace()
# print(a.expression.subs([("d_a",3),("a_a",2),("a_b",2)]))

a.subs([("a_a",1),("a_b",1),("d_a",1),("d_b",1)])
print(a.symbols(), a(2,2,2))

'''
class QuantumObject:
	def __init__(self, wavefunction=None):
		self.wavefunction = function
		

class SingleElectron(quantumObject):
	def __init__(self):
		super().__init__()

def Integrate(function):
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

# Caluclate S
S_11 = integrate(lambda r: phi_1s_alpha(r)*phi_1s_alpha(r))
S_12 = integrate(lambda r: phi_1s_alpha(r)*phi_1s_beta(r))
S_22 = integrate(lambda r: phi_1s_beta(r)*phi_1s_beta(r))

S = np.array([[S_11,S_12],[S_12,S_22]])
print(S)
'''
