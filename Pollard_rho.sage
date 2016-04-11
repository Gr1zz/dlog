def Pollard_rho(target, p, q, g):
	Fq = GF(q)
	x, a, b = GF(p)(1), Fq(0), Fq(0)
	x2, a2, b2 = x, a, b
	# iterative function used for Floyd's cycle detection
	def walk(x, a, b):
		if (int(x) % 3 == 0):
			return x*x, a*2, b*2
		elif (int(x) % 3 == 1):
			return x*target, a, b+1
		else:
			return x*g, a+1, b
	while (true):
		#print "xab=", x, a, b
		x, a, b = walk(x, a, b)
		x2, a2, b2 = walk(x2, a2, b2)
		x2, a2, b2 = walk(x2, a2, b2)
		#print "xab=", x, a, b
		#print x2, a2, b2
		if (x == x2):
			#print "yes=", x, x2
			r = b2-b
			if (gcd(r, q) != 1):
				return "Failure."
			result = (a-a2) / r
			return result
def Pollard(n):
	x = 2
	y = 2
	d = 1
	def g(x, n):
		return x**2+1%n
	while (d == 1):
		x = g(x, n)
		y = g(g(y, n), n)
		d = gcd(abs(x-y),n)
	if (d == n):
		return 0
	else:
		return d




