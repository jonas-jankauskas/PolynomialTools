load('Lib/CirclePolyTools.sage')

#Generates list of all possible {0,1} lists of a given length
def gen01(length):
    new_seqs = [[]]
    for j in range(length):
        new_seqs = ([seq+[0] for seq in new_seqs] + [seq+[1] for seq in new_seqs])
    return new_seqs

#Generates list of  all possible {0,1} lists of a given length that start with 1
def gen10(length):
    seqs = gen01(length-1)
    return [[1]+seq for seq in seqs]

#Generates list of  all possible {0,1} lists of a given length that start and end with 1
def gen101(length):
    seqs = gen01(length-2)
    return [[1]+seq+[1] for seq in seqs]


#Generates list of  all possible {0,1} lists that start and end with 1 of a given length excluding reversed and palindrome sequences
def gen101asym(length):
    seqs = gen10(length // 2)
    new_seqs = []
    if ((length % 2) == 0):
        for j in range(len(seqs)):
            for k in range(j+1, len(seqs)):
                new_seqs.append(seqs[j]+[c for c in reversed(seqs[k])])
    else:
        for j in range(len(seqs)):
            for k in range(j+1, len(seqs)):
                new_seqs.append(seqs[j]+[0]+[c for c in reversed(seqs[k])])
                new_seqs.append(seqs[j]+[1]+[c for c in reversed(seqs[k])])
    return new_seqs

#Mod 2 sublists

#Even-terms
def even_part(sequence):
    return sequence[::2]

#Odd-terms
def odd_part(sequence):
    return sequence[1::2]
	
#Mod 3 sublists
def zero_mod3_part(sequence):
	return sequence[::3]

def one_mod3_part(sequence):
	return sequence[1::3]
	
def two_mod3_part(sequence):
	return sequence[2::3]

#List of asymetric Newman polynomials of a given degree (reciprocated polynomials are excluded)
def gen_newman_asym(degree=0):
    return [R(seq) for seq in gen101asym(degree+1)]

#Generates a list of polynomials from their coefficient lists
def gen_poly_list(sequences):
    return [R(seq) for seq in sequences]

#Tests the values of sequence at points at z=1,-1, i, -i, exp(2*pi*i/3), -exp(2*pi*i/3)
def test_simple(sequence, min_val):
    s0 = even_part(sequence); a = sum(s0)
    s1 = odd_part(sequence); b = sum(s1)
    s00 = even_part(s0); aa = sum(s00)
    s01 = odd_part(s0); ab = sum(s01)
    s10 = even_part(s1); ba = sum(s10)
    s11 = odd_part(s1); bb = sum(s11)
	
    t0 = zero_mod3_part(sequence); c = sum(t0)
    t1 = one_mod3_part(sequence); d = sum(t1)
    t2 = two_mod3_part(sequence); e = sum(t2)
    
    t00 = even_part(t0); ca = sum(t00)
    t01 = odd_part(t0); cb = sum(t01)
    t10 = even_part(t1); da = sum(t10)
    t11 = odd_part(t1); db = sum(t11)
    t20 = even_part(t2); ea = sum(t20)
    t21 = odd_part(t2); eb = sum(t21)
    
    s = a+b
    ss = abs(a - b)
    sss = (aa-ab)^2 + (ba-bb)^2
	
    t = c^2 + d^2 + e^2 - c*d - d*e - e*c
    tt = (ca-cb)^2 + (da-db)^2 + (ea-eb)^2 - (ca-cb)*(da-db) - (da-db)*(ea-eb) - (ea-eb)*(ca-cb)
	
    return ((s >= min_val) and (ss >= min_val) and (sss >= min_val^2) and (t >= min_val^2) and (tt >= min_val^2))
#insert
#and (s <= len(sequence)-3)

def get_min_val(f, points):
    return min([abs(f(t)) for t in points])

def filter_min_val(poly, points, min_val):
    return [p for p in poly if (get_min_val(p[1], points) > min_val)]

def iterate_fairey(sequence=[0, 1]):
    new_seq = []
    for j in range(len(sequence)-1):
        new_seq.append(sequence[j])
        new_seq.append((sequence[j].numerator()+sequence[j+1].numerator())/(sequence[j].denominator()+sequence[j+1].denominator()))
    new_seq.append(sequence[len(sequence)-1])
    return new_seq

def iterate_fairey_sym(sequence=[-1, 0, 1]):
    new_seq = [c for c in sequence if c >=0]
    new_seq = iterate_fairey(new_seq)
    new_seq = [-c for c in reversed(new_seq)] + new_seq[1:]
    return new_seq

def filter_simple(sequences, min_val):
    return [seq for seq in sequences if test_simple(seq, min_val)]
	
def filter_by_square_amplitude(N, sequences, min_sq):
    new_seqs = []
    ft = FastFourierTransform(N)
    for seq in sequences:
        for j in range(len(seq)):
		    ft[j] = seq[j]
        for j in range(len(seq), N):
            ft[j] = 0
        ft.forward_transform()
        if min([abs(CDF(w)) for w in ft]) >= min_sq:
            new_seqs.append(seq)
    return new_seqs

def gen_Littlewood(degree):
    littlewood_list = [1+0*x]
    for j in range(degree):
        littlewood_list = littlewood_list + [x^(j+1)+ f for f in littlewood_list] + [x^(j+1)- f for f in littlewood_list]
    return litlewood_list








