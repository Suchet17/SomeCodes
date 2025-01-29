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


#Read sequences
seq1 = "*"+"".join([i.strip() for i in open("Sequence1.txt","r").readlines()])
seq2 = "*"+"".join([i.strip() for i in open("Sequence2.txt","r").readlines()])

#Readi Substitution Matrix, Define Gap Penalty
mat = bl(62)
d = 4 #Keep this positive, it subtracts later

#Create Dynamic programming Table
m = len(seq1)
n = len(seq2)
table = np.zeros((m,n))
back = [[[] for i in range(n)] for i in range(m)]

#Base Case
for i in range(1,m):
	back[i][0] = [(0,0)]
for i in range(1,n):
	back[0][i] = [(0,0)]
back[0][0] = [(0,0)]
start = (0,0)


#Fill Table with bactracking arrows
for i in range(1,m):
	for j in range(1,n):
		s1 = table[i-1][j-1] + mat[seq1[i]][seq2[j]] #Diagonal
		s2 = table[i][j-1] - d #Vertical
		s3 = table[i-1][j] - d #Horizontal
		mx = max(s1, s2, s3)
		back[i][j] = []
		if mx >= table[start[0]][start[1]]:
			start = (i,j)
		if mx <= 0:
			back[i][j].append((0,0))
			table[i][j] = 0
			continue
		else:
			table[i][j] = mx
		if s1 == mx:
			back[i][j].append((i-1,j-1))
		if s2 == mx:
			back[i][j].append((i,j-1))
		if s3 == mx:
			back[i][j].append((i-1,j))
del i, j, mx, s1, s2, s3


traceback = [[(start[0],start[1])]]
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
	score = 0
	match = 0
	mismatch = 0
	gaps = 0
	mem = t[-1]
	
	for i in t[::-1][1:]:
		if (mem[0]+1, mem[1]) == i:
			aseq1c+=seq1[i[0]]
			aseq2c+='-'
			score -= d
			gaps+=1
		elif (mem[0], mem[1]+1) == i:
			aseq1c+='-'
			aseq2c+=seq2[i[1]]
			score -= d
			gaps+=1
		elif (mem[0]+1, mem[1]+1) == i:
			if seq1[i[0]] == seq2[i[1]]:
				match+=1
			else:
				mismatch+=1
			score += mat[seq1[i[0]]][seq2[i[1]]]
			aseq1c+=seq1[i[0]]
			aseq2c+=seq2[i[1]]
		elif i == (0,0):
			print(i,mem)
		mem = i
	
	print(aseq1c)
	print(aseq2c)
	print("Score = ", score, end = ',\t')
	print("Identical Matches = ", match, end = ',\t')
	print("Mismatches = ", mismatch, end = ',\t')
	print("Gaps = ", gaps, end = '\n\n')