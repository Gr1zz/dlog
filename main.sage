import time

load("dlog.sage")
load("Pollard_rho.sage")
load("index_calculus.sage")
DEBUG = False
VERBOSE = True

listPrimes = [
	[1319, 659, 5], #0, ~11 bits
	[6659, 3329, 5], #1, ~13 bits
	[32603, 16301, 3], #2, ~15 bits
	[99839, 49919, 4], #3, ~17 bits
	[312839, 156419, 4], #4, ~19 bits
	[1785803, 892901, 4], #5, ~21 bits
	[6771419, 3385709, 3], #6, ~23 bits
	[27612779, 13806389, 3], #7, ~25 bits
	[111619679, 55809839, 2], #8, ~27 bits
	[526361819, 263180909, 3], #9, ~29 bits
	[2111297939, 1055648969, 4], #10, ~31 bits
	[5848937243, 2924468621, 3], #11, ~33 bits
	[24570203447, 12285101723, 2], #12, ~35 bits
	[102595058507, 51297529253, 3], #13, ~37 bits
	[291849846923, 145924923461, 3], #14, ~39 bits
	[1232673178823, 616336589411, 4], #15, ~41 bits
	[5707721278283, 2853860639141, 3], #16, ~43 bits
	[22464712289039, 11232356144519, 2], #17, ~45 bits
	[125223411188483, 62611705594241, 3], #18, ~47 bits
	[499957513308299, 249978756654149, 3], #19, ~49 bits
	[2073681093605567, 1036840546802783, 4], #20, ~51 bits
	[8949958374119339, 4474979187059669, 3], #21, ~53 bits
	[25724961786521087, 12862480893260543, 2] #22, ~55 bits
]


def tests(algo="IC", fieldRange=[], repetitions=1):
	""" Tests a set of targets.
	INPUT :

	* "algo" -- can be "Kraitchik" (default), "BSGS", "Pollard_rho" or "bruteforce"
	* "fieldRange" -- the range of the fields we want to test. Full field list by default.

	EXAMPLES :

	sage: DEBUG=True; tests("Kraitchik")

	sage: tests("BSGS", [2, 11, 4])

	sage: tests("Pollard_rho", range(8,11))

	"""
	targets =[173114611229025298050625,
 620015254827899774343889,
 3571738524033324718284324,
 19983674595088883924731161,
 50842279960751853178444161,
 54431372065519868346301284,
 107817957581486271163822969,
 235606539908871452207227441
]
	#targets = findTargets(5)
	for i in range(0, repetitions):
		test(algo, fieldRange)

def test(algo="IC", fieldRange=[]):
	times = {}
	if fieldRange == []:
		fieldRange = range (0, len(listPrimes))
	for i in fieldRange:
		F = GF(listPrimes[i][0])
		target = F.random_element()^2
		t0 = time.clock()
		if (algo == "bruteforce"):
			result = bruteforce(F(target), listPrimes[i][0], listPrimes[i][1], listPrimes[i][2])
		elif (algo == "BSGS"):
			result = BSGS(F(target), listPrimes[i][0], listPrimes[i][1], listPrimes[i][2])
		elif (algo == "Pollard_rho"):
			result = Pollard_rho(F(target), listPrimes[i][0], listPrimes[i][1], listPrimes[i][2])
		elif (algo == "IC"):
			result = index_calculus(F(target), listPrimes[i][0], listPrimes[i][1], listPrimes[i][2])
		else:
			return "Unknown algo"

		timer = time.clock() - t0
		times[listPrimes[i][0]] = [result, timer]
		if (VERBOSE):
			print "log(", target, ") = ", result
			if (DEBUG):
				print  "verif :", discrete_log_rho(F(target), F(listPrimes[i][2]))
			print "[p,q,r] =", listPrimes[i]
			print "time =", timer
			print "-------------------------------\n"
	formatResults(times)

def formatResults(list):
	sortedList = sorted(list.items(), key=operator.itemgetter(0))
	#print "Bits\t\tTime\t\tResult" 
	print "Times"	
	for i in range(0, len(sortedList)):
		item = sortedList[i]
		bits = ceil(log_b(item[0],2))
		time = round(item[1][1],5)
		#print bits, ",\t\t", time, ",\t\t", item[1][0] 
		print time
	print "\n"

def findTargets(n):
	"""
	"""
	l = set()
	testedTargets = set()
	while (len(l) < n):
		while (true):
			testTarget = ceil(random()*listPrimes[len(listPrimes)-1][0])
			if not (testTarget in testedTargets):
				break
		j = 0
		while (j < len(listPrimes)):
			K = GF(listPrimes[j][0], modulus="primitive")
			if (K(testTarget) != 0):
				order = K(testTarget^2).multiplicative_order()
			else:
				break
			if (order != listPrimes[j][1]):
				break
			j=j+1
		if (j == len(listPrimes)):
			l.add(testTarget^2)
	return l

# Generates prime numbers of size k in base 'base'
def generate(k, base):
	while (true):
		r = base^(k-1) + randint(0, base^(k-1))
		q = next_prime(r)
		p = 2*q+1
		if (is_prime(p)):
			break
	g = 2
	while (true):
		if (pow(g, q, p) == 1):
			return [p,q,g]
		else:
			g = g+1

