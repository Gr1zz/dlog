import time
DEBUG = True

listPrimes = [
	[1319, 659, 5], # ~10 bits
	[6659, 3329, 5], # ~12 bits
	[32603, 16301, 3], # ~14 bits
	[99839, 49919, 4], # ~16 bits
	[312839, 156419, 4], # ~18 bits
	[1785803, 892901, 4], # ~20 bits
	[66999767, 33499883, 4], # ~25 bits
	[2111297939, 1055648969, 4], # ~30 bits
	[47260079003, 23630039501, 4] # ~35 bits
	#[1232673178823, 616336589411, 4], # ~40 bits
	#[39476820129707, 19738410064853, 4], # ~45 bits
	#[2073681093605567, 1036840546802783, 4] # ~50 bits

]

""" Tests a set of computed targets.
    You can remove some, but remember that Pollard_rho gives very different
    results depending on the target. """
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
		result = Pollard_rho(target, listPrimes[i][0], listPrimes[i][1], listPrimes[i][2])
		timer = time.clock() - t0
		times[listPrimes[i][0]] = [result, timer]
		if (DEBUG):
			print "log(", target, ") = ", result
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

def Pollard_rho(target, p, q, g):
	Fq = GF(q)
	[x, a, b] = [GF(p)(1), Fq(0), Fq(0)]
	[x2, a2, b2] = [x, a, b] #f(xab, target, p, g)

	# iterative function used for Floyd's cycle detection
	def walk(x, a, b):
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
