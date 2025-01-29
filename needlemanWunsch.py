import numpy as np
from blosum import BLOSUM as bl #Import Substitution Matrix

def linearize(l): #Flatten Nested List
	l2 = []
	for i in l:
		for j in i:
			if type(j) == list:
				l2.append(j)
			else:
				l2.append(i)
				break
	return l2

def align(seq1, seq2, mat, d):
	#mat - Substitution Matrix
	d = abs(d) # Linear gap penalty
	#Create Dynamic programming Table
	m = len(seq1)
	n = len(seq2)
	table = np.zeros((m,n))
	back = [[[] for i in range(n)] for i in range(m)]
	
	#Base Case
	for i in range(1,m):
		table[i][0] = -i*d
		back[i][0] = [(i-1,0)]
	for i in range(1,n):
		table[0][i] = -i*d
		back[0][i] = [(0,i-1)]
	table[0][0] = 0
	back[0][0] = [(0,0)]
	
	#Fill Table with bactracking arrows
	
	for i in range(1,m):
		for j in range(1,n):
			s1 = table[i-1][j-1] + mat[seq1[i]][seq2[j]] #Diagonal
			s2 = table[i][j-1] - d #Vertical
			s3 = table[i-1][j] - d #Horizontal
			mx = max(s1, s2, s3)
			table[i][j] = mx
			back[i][j] = []
			if s1 == mx:
				back[i][j].append((i-1,j-1))
			if s2 == mx:
				back[i][j].append((i,j-1))
			if s3 == mx:
				back[i][j].append((i-1,j))
	del i, j, mx, s1, s2, s3
	
	traceback = [[(m-1,n-1)]]
	flag = True
	while flag:
		flag = False
		for i in range(len(traceback)):
			last = traceback[i][-1]
			if last == (0,0):
				continue
			b = back[traceback[i][-1][0]][traceback[i][-1][1]]
			if len(b) == 1:
				traceback[i] = traceback[i] + b
			else:
				t = []
				for j in b:
					t.append(traceback[i] + [j])
				traceback[i] = t
				traceback = linearize(traceback)
		for i in traceback:
			if i[-1] != (0,0):
				flag = True
	
	
	count = len(traceback)
	
	for t in traceback:	
		aseq1c = ""
		aseq2c = ""
		seq1c = seq1[1:]
		seq2c = seq2[1:]
		score = 0
		match = 0
		mismatch = 0
		gaps = 0
		mem = t[-1]
		
		for i in t[::-1][1:]:
			if (mem[0]+1, mem[1]) == i:
				aseq1c+=seq1c[0]
				seq1c = seq1c[1:]
				aseq2c+='-'
				score -= d
				gaps+=1
			elif (mem[0], mem[1]+1) == i:
				aseq1c+='-'
				aseq2c+=seq2c[0]
				seq2c = seq2c[1:]
				score -= d
				gaps+=1
			elif (mem[0]+1, mem[1]+1) == i:
				if seq1c[0] == seq2c[0]:
					match+=1
				else:
					mismatch+=1
				score += mat[seq1c[0]][seq2c[0]]
				aseq1c+=seq1c[0]
				seq1c = seq1c[1:]
				aseq2c+=seq2c[0]
				seq2c = seq2c[1:]
			else:
				print(i)
			mem = i
		
		
		print(aseq1c)
		print(aseq2c)
		print("Score = ", score, end = ',\t')
		print("Identical Matches = ", match, end = ',\t')
		print("Mismatches = ", mismatch, end = ',\t')
		print("Gaps = ", gaps, end = '\n\n')
	return score


seq1 = "*"+"".join([i.strip() for i in open("Sequence1.txt","r").readlines()])
seq2 = "*"+"".join([i.strip() for i in open("Sequence2.txt","r").readlines()])
align(seq1,seq2,bl(62),4) #BLOSUM62 with -4 linear gap  penalty