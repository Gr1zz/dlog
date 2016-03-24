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
	algoCounter = 2^20
	nextFactor = 1
	while not (n == 1):
		if (algoCounter > 2):
			nextFactor = trial_division(n, B)
		else:
			if (algoCounter == 2):
				algo = "P-1"
			elif (algoCounter == 1):
				algo = "P+1"
			else:
				algo = "ECM"
			nextFactor = e.one_curve(n, B1=B, algorithm=algo)[0]	
		if (nextFactor > B):
			return [false]
		if (is_prime(nextFactor)):
			n = n.divide_knowing_divisible_by(nextFactor)
			factors[nextFactor] += 1
		algoCounter -= 1
			
	return [true,factors]

def is_smooth_trial(n, B):
	factors = Counter()
	while not (n == 1):
		nextFactor = trial_division(n, B)
		if (nextFactor > B):
			return [false]
		n = n.divide_knowing_divisible_by(nextFactor)
		factors[nextFactor] += 1
	return [true,factors]

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
