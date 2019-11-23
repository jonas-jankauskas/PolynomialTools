#This SAGE script contains mathematical functions used to
#calculate number of roots of integer polynomial inside/on the unit disk
#maximas/minimas on the unit circle, and other things

def rec(f):
    return R([c for c in reversed(f.list())])

def diff_abs_f2_g2(f, g):
    return f*rec(f)*x^max(g.degree()-f.degree(), 0) - g*rec(g)*x^max(f.degree()-g.degree(), 0)
	
def selfrec2chebyshev(f):
	k = f.degree() // 2
	cfs = f.list()
	g = 0*x+cfs[k]
	for j in range(k):
		g += 2*cfs[j]*chebyshev_T(k-j, x)
	return g
	
def f2_2chebyshev(f): 
    return selfrec2chebyshev(f*rec(f))
	
def min11(f):
    if f.degree() > 0:
        critical_points = [t for t in (f.derivative()).roots(ring=RR, multiplicities=False) if (abs(t) <= 1)]
    else:
        critical_points = []
    return min([abs(f(1)), abs(f(-1))] + [abs(f(t)) for t in critical_points])

def max11(f):
    if f.degree() > 0:
        critical_points = [t for t in (f.derivative()).roots(ring=RR, multiplicities=False) if (abs(t) <= 1)]
    else:
        critical_points = []
    return max([abs(f(1)), abs(f(-1))] + [abs(f(t)) for t in critical_points])

def circle_min(f):
    return sqrt(min11(f2_2chebyshev(f)))

def circle_max(f):
    return sqrt(max11(f2_2chebyshev(f)))

def rotate(sequence, shift=1):
    return sequence[-shift:] + sequence[:-shift]

def all_rotations(sequence):
    rotations = []
    for j in range(len(sequence)):
        rotations.append(rotate(sequence, j))
    return rotations

def skew(sequence):
    n = len(sequence)
    new_seq = [(-1)^j*sequence[n-j-1] for j in range(n)]
    return new_seq
	
def num_real_zeros(f):
	return sum([w[1] for w in f.roots(ring=RR)])

#warning: if f has a zero at x=-1 or x=1, the function might report the count incorrectly
def num_real_zeros11(f):
	return sum([w[1] for w in f.roots(ring=RR) if abs(w[0]) < 1])

def filter_nonreciprocal(poly_list):
    new_list = []
    for f in poly_list:
        g = R([c for c in reversed(f.list())])
        if not((f == g) or (g in new_list)):
            new_list.insert(len(new_list), f)
    return new_list


def clear11(f):
	if f == 0:
		m = +Infinity
	else:
		m = 0
		while (f.quo_rem(x-1)[1] == 0):
			f = f.quo_rem(x-1)[0]
			m +=1
		while (f.quo_rem(x+1)[1] == 0):
			f = f.quo_rem(x+1)[0]
			m +=1
	return (f, m)
	
def numzeros_selfrec(f):
		g, u = clear11(f);
		u += 2*num_real_zeros11(selfrec2chebyshev(g))
		return((f.degree()-u)//2, u)
		
def numzeros(f):
	g = f.gcd(rec(f))
	n, u = numzeros_selfrec(g)
	h = f.quo_rem(g)[0]
	n += sum([w[1] for w in h.roots(ring=CC) if abs(w[0]) < 1])
	return(n, u)

#Deprecated: zeros on the boundary of the disk might no be accounted proprely using real arithmetics	

# def num_zeros(f):
    # return  sum([w[1] for w in f.roots(ring=CC, multiplicities=True) if abs(w[0]) < 1])




