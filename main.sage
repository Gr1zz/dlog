import time
load("dlog.sage")
load("Pollard_rho.sage")
load("Kraitchik.sage")
DEBUG = True
VERBOSE = True

listPrimes = [
	[1319, 659, 5], # ~10 bits
	[6659, 3329, 5], # ~12 bits
	[32603, 16301, 3], # ~14 bits
	[99839, 49919, 4], # ~16 bits
	[312839, 156419, 4], # ~18 bits
	[1785803, 892901, 4], # ~20 bits
	[66999767, 33499883, 4], # ~25 bits
	[2111297939, 1055648969, 4], # ~30 bits
	[47260079003, 23630039501, 4], # ~35 bits
	[1232673178823, 616336589411, 4] # ~40 bits
	#[39476820129707, 19738410064853, 4], # ~45 bits
	#[2073681093605567, 1036840546802783, 4] # ~50 bits

]

# Tests a set of computed targets.
def tests():
	targets =[173114611229025298050625,
 620015254827899774343889,
 3571738524033324718284324,
 19983674595088883924731161,
 50842279960751853178444161,
 54431372065519868346301284,
 107817957581486271163822969,
 235606539908871452207227441
]
	for i in range(0, len(targets)):
		test(targets[i])

def test(target):
	times = {}
	for i in range (0, len(listPrimes)):
		F = GF(listPrimes[i][0])
		t0 = time.clock()
		result = Pollard_rho(F(target), listPrimes[i][0], listPrimes[i][1], listPrimes[i][2])
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
	print "Bits\t\tTime\t\tResult" 
	for i in range(0, len(sortedList)):
		item = sortedList[i]
		bits = ceil(log_b(item[0],2))
		time = round(item[1][1],5)
		print bits, ",\t\t", time, ",\t\t", item[1][0] 
	print "\n"

# trouve n cibles (si possible), entre min et max.
# findTargets(7,2,30) = {3, 4, 9, 12, 16, 25, 27}
def findTargets(n):
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

# génère des nombres premiers aléatoire de taille k en base "base"
def gen(k, base):
	while (true):
		r = base^(k-1) + randint(0, base^(k-1))
		q = next_prime(r)
		p = 2*q+1
		if (is_prime(p)):
			break
	while (true):
		g = GF(p).random_element()
		if (g^q == 1):
			break
	return([p,q,g])

