from sets import Set
import random
import os

ymax = 491

def getXM():
	f = open('mhad435.10m', 'r')
	x = []
	for line in f.readlines():
		if line.startswith('%'):
			continue
		y = [float(value) for value in line.split()]
		x.append(y)
		
	f.close()
	return x

def constructXYM(srcx):
	li = getXM()
	returnlist = []
	for x,m in li:
		y = float(random.randrange(0, ymax))/10
		returnlist.append([srcx+x,y,m])
	return returnlist

def getRandomPairs(srcx, k):
	li = constructXYM(srcx)
	tmpSet = Set()
	res = []
	for i in range(2*k):
		while True:
			idx = random.randrange(len(li))			
			if idx not in tmpSet:
				res.append(li[idx])
				break
	writeToFile(res)

def writeToFile(li, out='out.txt'):
	print len(li)
	f = open(out, 'w')
	for i in li:
		f.write(' '.join([str(j) for j in i]))
		f.write('\n')
	f.close()
