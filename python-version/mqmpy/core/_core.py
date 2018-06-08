'''
General math objects
'''

import sympy as sp

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

