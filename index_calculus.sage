from collections import Counter
import random

def is_smooth_ecm(n, B):
	factors = Counter()
	e = ECM()
	nextFactor = 1
	algoCounter = 0
	bound = 2^20
	trials = is_smooth_trial(n, bound)
	firstFactors = trials[1]
	n = trials[2]
	print trials
	print "bound=", bound
	bound = e.recommended_B1(len(str(n)))
	print "bound=", bound
	while not (n == 1):
		if (algoCounter == 0):
			algo = "P-1"
		elif (algoCounter == 1):
			algo = "P+1"
		else:
			algo = "ECM"
		print "#", algoCounter, algo
		nextFactor = e.one_curve(n, B1=bound, algorithm=algo)[0]
		print "next = ", nextFactor, "n=", n, "B1=", bound
		isPrime = is_prime(nextFactor)
		if (isPrime):
			if (nextFactor > B):
				return [false]
			else:
				n = n.divide_knowing_divisible_by(nextFactor)
				factors[nextFactor] += 1
		else:
			lastFactors = is_smooth_trial(n, B)
			if (lastFactors[0]):
				factors.update(lastFactors[1])
			else:
				return [false]
		algoCounter += 1
			
	return [true,factors]

def is_smooth_trial(n, B):
	factors = Counter()
	while not (n == 1):
		nextFactor = trial_division(n, B)
		if (nextFactor > B):
			return [false, factors, n]
		nextFactorPow = valuation(n, nextFactor)
		n = n.divide_knowing_divisible_by(nextFactor^nextFactorPow)
		factors[nextFactor] += nextFactorPow
	return [true,factors, n]

def is_smooth(n, B):
	""" Naive version of smoothness detection
	If you want to use this function, you have to adapt Kraitchik's code
	e.g. :
		i = indexes[p_i[0]]
		relations[k, i] = Fq(p_i[1])
	"""
	factors = factor(n)
	if (factors[len(factors)-1][0] > B):
		return [false]    
	return [true, factors]

def Kraitchik(target, p, q, g):
	""" Returns the discrete log of target in base g.
	INPUT :

	* "target" -- the target, e.g a Diffie-Hellman public key
	* "p" -- the modulus of the group
	* "q" -- the order of g, belongs to GF(p)
	* "g" -- a generator of GF(p)
	"""
	Fp = GF(p)
	Fq = GF(q)
	B = ceil(exp(0.5*sqrt(2*log(p)*log(log(p)))))
	base = list(primes(B+1))
	# Precompute indexes
	indexes = {prime:base.index(prime) for prime in base}
	S = len(base)
	relations = matrix(Fq, S+1, S, sparse=True)
	min = ceil(log(q))
	max = ceil(sqrt(q))
	row = []
	k = 0
	while (k < S+1):
		while (true):
			a = Fq.random_element()#(random.randrange(min, max))
			b = Fq.random_element()#(random.randrange(min, max))
			if not (a,b) in row:
				break
		# Fast modular exponentiation
		z = Fp(g)^a*Fp(target)^b
		#z = g^a * target^b % p
		isSmooth = is_smooth_trial(ZZ(z), B)
		if (isSmooth[0]):
			row.append((a,b))
			for p_i in isSmooth[1]:
				#i = indexes[p_i[0]]
				#relations[k, i] = Fq(p_i[1])
				i = indexes[p_i]
				relations[k, i] = Fq(isSmooth[1][p_i])
			k = k+1
	ker = relations.left_kernel().basis()[0]
	Z = 1 ; A = 0 ; B = 0
	for ker_i, row_i in zip(ker, row):
		A += ker_i*row_i[0]
		B += ker_i*row_i[1]
	return -A*Fq(B**-1)

def power_m(g_base, a, p_mod):
  x=1
  bits = "{0:b}".format(a)
  for i, bit in enumerate(bits):
    if bit=='1': x = (((x**2)*g_base)%p_mod)
    elif bit=='0': x = ((x**2)%p_mod)
  return x%p_mod
