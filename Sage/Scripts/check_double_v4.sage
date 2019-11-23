#mathematical tools for calculations of roots of polynomials
load("../Lib/CirclePolyTools.sage")

#save file directory and filename prefix
savedir = 'Sage/'
savename = 'SPN'

#pattern data - I enter this manually for each pattern
prefix = [1]
rep1 = [1]
padd = []
rep2 = [-1]
suffix = []


#maximal number of repetitions of first and second periods
#the table size (no. of rows and columns), also controls the speed :) 

max1 = 15
max2 = 15

#matrix that stores pairs (N(p), U(p)) computed for polynomial p for each cell in the table
rootmatrix = []
#matrix that stores increments in N(p) as you go 1 step down or 1 step right 
diffmatrix = []
#error matrix that stores the difference (dN(p)/dm, dN(p)/dn) between the guessed formula and the observed value of N(p)
errormatrix = []

#most frequenquent and second most frequent observed errors 
maxdiff = (0, 0)
oldmaxdiff = maxdiff
diffreqs = {maxdiff: 0}

maxerr = 0
oldmaxerr = maxerr
errfreqs = {maxerr: 0}

#pattern to text string conversions
def pm2str(pattern):

	s = ''
		
	for k in pattern:
		if k==-1:
			s = s + '-'
		elif k==1:
			s = s + '+'
	
	return s
	
def zo2str(pattern):

	s = ''
		
	for k in pattern:
		if k==0:
			s = s + '0'
		elif k==1:
			s = s + '1'
	
	return s
	
def pattern2str(prefix, rep1, padd, rep2, suffix):

	s = pm2str(prefix);
	if len(rep1) > 0:
		s += '(' + pm2str(rep1)+ ')'
		
	s += pm2str(padd)
	
	if len(rep2) > 0:
		s += '(' + pm2str(rep2)+ ')'
	
	s += pm2str(suffix)
	
	return s

#data file to save the results
savefile = open(savedir + savename + '_' + '_' + pattern2str(prefix, rep1, padd, rep2, suffix)+'.txt', 'w');


#count roots inside/on the unit disk
for m in range(0, max1+1):
	rootmatrix.append([])
	for n in range(0, max2+1):
		pattern = prefix + rep1*m + padd + rep2*n+suffix
		poly = R(pattern)
		rootmatrix[m].append(numzeros(poly))

#calculate step-down and step-right increments in N(p)
for m in range(0, max1):
	diffmatrix.append([])
	for n in range(0, max2):
		diffmatrix[m].append((rootmatrix[m+1][n][0] - rootmatrix[m][n][0], rootmatrix[m][n+1][0] - rootmatrix[m][n][0]))
		if diffmatrix[m][n] in diffreqs:
			diffreqs[diffmatrix[m][n]] += 1
		else:
			diffreqs[diffmatrix[m][n]] = 1
		if diffreqs[diffmatrix[m][n]] > diffreqs[maxdiff]:
			oldmaxdiff = maxdiff
			maxdiff = diffmatrix[m][n]

#find the most/second most frequent ocurring increments, guess the linear equation and calculate errors			
for m in range(0, max1+1):
	errormatrix.append([])
	for n in range(0, max2+1):
		errormatrix[m].append(rootmatrix[m][n][0] - maxdiff[0]*m - maxdiff[1]*n)
		if errormatrix[m][n] in errfreqs:
			errfreqs[errormatrix[m][n]] += 1
		else:
			errfreqs[errormatrix[m][n]] = 1
		if errfreqs[errormatrix[m][n]] > errfreqs[maxerr]:
			oldmaxerr = maxerr
			maxerr = errormatrix[m][n]

#write the results
savefile.write('k-matrix= \n')
for row in rootmatrix:
	savefile.write('[')
	for element in row:
		savefile.write(' '+str(element[0]))
	savefile.write(']\n')
savefile.write('\n')

savefile.write('u-matrix= \n')
for row in rootmatrix:
	savefile.write('[')
	for element in row:
		savefile.write(' '+str(element[1]))
	savefile.write(']\n')
savefile.write('\n')

savefile.write('dk/dm-matrix= \n')
for row in diffmatrix:
	savefile.write('[')
	for element in row:
		savefile.write(' '+str(element[0]))
	savefile.write(']\n')
savefile.write('\n')

savefile.write('dk/dl-matrix= \n')
for row in diffmatrix:
	savefile.write('[')
	for element in row:
		savefile.write(' '+str(element[1]))
	savefile.write(']\n')
savefile.write('\n')

savefile.write('Diffreqs   = '+str(diffreqs)+'\n')
savefile.write('Maxdiff    = '+str(maxdiff)+': '+str(diffreqs[maxdiff])+'\n')
savefile.write('Oldmaxdiff = '+str(oldmaxdiff)+': '+str(diffreqs[oldmaxdiff])+'\n')

savefile.write('\n error-matrix= \n')
for row in errormatrix:
	savefile.write('[')
	for element in row:
		savefile.write(' '+str(element-maxerr))
	savefile.write(']\n')
savefile.write('\n')

savefile.write('Maxerr    ='+str(maxerr)+': '+str(errfreqs[maxerr])+'\n')
savefile.write('Oldmaxerr ='+str(oldmaxerr)+': '+str(errfreqs[oldmaxerr])+'\n')

savefile.write('Zeros  eq.:  k = ' + str(maxdiff[0])+ '*m' + ' + ' + str(maxdiff[1])+ '*l' + ' + ' + str(maxerr)+'\n')
savefile.write('Degree eq.:  n = ' + str(len(rep1))+ '*m' + ' + ' + str(len(rep2))+ '*l' + ' + ' + str(len(prefix)+len(padd)+len(suffix)-1)+'\n')
savefile.write('Determinant: det  = ' + str(maxdiff[0]*len(rep2)-maxdiff[1]*len(rep1))+'\n')
#savefile.write('Pattern no.       = ' + str(l)+'\n')		
savefile.close()
