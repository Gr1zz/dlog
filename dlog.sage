# Literally a bruteforce
def bruteforce(target, p, q, g):
	h = 1
	i = 0
	while (h != target):
		i=i+1
		h = h*g % p
	return i

# Literally a bruteforce with a space-time trade-off.
def BSGS(target, p, q, g):
	m = ceil(sqrt(q))
	l = {}
	fast_pow = 1
	for j in xrange(0, m):                                 
		l[fast_pow] = j
		fast_pow = fast_pow * g % p
	c = g^-m % p
	y = target
	for i in xrange(0, m):
		if l.has_key(y):		
			return l[y] + i*m % q
		else:
			y = y*c % p
	return "Failed"

