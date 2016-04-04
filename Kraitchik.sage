from collections import Counter
import random

testSmooth = [435682255480385335523828724249141091727624667875593489604343,
	3663840832853790873152444875116402874459315374909942896789450,
	654652562720341179717613494638996718646726003058221250978957,
	587072201011339627968436784943805434544568356060797599113084,
	237727123467893239512117415871911714087369865229864581775858,
	71452117609987925031778992261150994093558813445367719069074,
	102586112126343194451136748240085342925284766227291358170664,
	612284158337267326194396004206752864417216479008382232680127,
	168997433684626505415773457671544968917142822182742003555113,
	546097631335380720506264048780419976944832913083173983886747,
	654515694029474037890397753393110284433628162118803758847121,
	28686269421783111864294608644064015953031760727908881675314110
]

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
		if (nextFactor > B):
			return [false]
		if (is_prime(nextFactor)):
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

def foo():
	targets =[173114611229025298050625,
 620015254827899774343889,
 3571738524033324718284324,
 19983674595088883924731161,
 50842279960751853178444161,
 54431372065519868346301284,
 107817957581486271163822969,
 235606539908871452207227441
]
	Fq = GF(1036840546802783)
	t0 = time.time()
	max = len(targets)
	e = ECM()
	for i in range(0, max):
		e.one_curve(targets[i], B1=2^30)
		print i
	print time.time()-t0

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
		isSmooth = is_smooth_ecm(ZZ(z), B)
		if (isSmooth[0]):
			row.append((a,b))
			for p_i in isSmooth[1]:
				i = indexes[p_i]
				relations[k, i] = Fq(isSmooth[1][p_i])
			k = k+1
	ker = relations.left_kernel().basis()[0]
	Z = 1 ; A = 0 ; B = 0
	for ker_i, row_i in zip(ker, row):
		A += ker_i*row_i[0]
		B += ker_i*row_i[1]
	return -A*Fq(B**-1)
