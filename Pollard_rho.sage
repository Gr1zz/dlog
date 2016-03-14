def Pollard_rho2(target, p, q, g):
	Fq = GF(q)
	[x, a, b] = [GF(p)(1), Fq(0), Fq(0)]
	[x2, a2, b2] = [x, a, b] #f(xab, target, p, g)

	# iterative function used for Floyd's cycle detection
	def walk(x, a, b):
		if (Mod(x,3) == 0):
			x *= x
			a *= 2
			b *= 2
		if (Mod(x,3) == 1):
			x *= target
			b = b+1
		if (Mod(x,3) == 2):
			x *= g
			a = a+1
		return [x, a, b]
	def walk2(x, a, b):
		if (Mod(x,3) == 0):
			return [x*x, a*2, b*2]
		if (Mod(x,3) == 1):
			return [x*target, a, b+1]
		if (Mod(x,3) == 2):
			return [x*g, a+1, b]
	while (true):
		[x,a,b] = walk(x, a, b)
		[x2, a2, b2] = walk(x2, a2, b2)
		[x2, a2, b2] = walk(x2, a2, b2)
		if (x == x2):
			r = b2-b
			if (gcd(r, q) != 1):
				return "failure"
			result = (a-a2) / r % q
			return result

def Pollard_rho(public_key, order, generator, modulus, a=0, b=0):
    """ This is the implementation of the algorithm introduced in
    the Handbook of Applied Cryptography chapter 3.6.3
    It is definitely not the latest improvement, 
    nor the parallelizable version.
    You need to know the order of the generator to use it
    """
    # initialization
    alpha = generator # to keep the same variables as in the book
    beta = public_key

    if a != 0 or b != 0: # <- wrong
        x = Mod(power_mod(alpha, a, modulus) * power_mod(beta, b, modulus), modulus)
    else:
        x = Mod(1, modulus)

    x = [x, x]
    a = Mod(a, order)
    a = [a, a]
    b = Mod(b, order)
    b = [b, b]
    
    # iteration function
    def iteration(x, a, b):
        if Mod(x, 3) == 1: # x in S_1 (chosen from example)
            x = beta * x
            b = b + 1

        elif Mod(x, 3) == 0: # x in S_2
            x = x * x
            a = 2 * a
            b = 2 * b

        else: # x in S_3
            x = alpha * x
            a = a + 1

        return x, a, b

    # loop
    while True:
        # iteration function
        x[0], a[0], b[0] = iteration(x[0], a[0], b[0])

        x[1], a[1], b[1] = iteration(x[1], a[1], b[1])
        x[1], a[1], b[1] = iteration(x[1], a[1], b[1])

        # detect collision
        if x[0] == x[1]:
            r = b[0] - b[1]
            if r != 0:
                return r^-1 * (a[1] - a[0])
            else:
                break

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
