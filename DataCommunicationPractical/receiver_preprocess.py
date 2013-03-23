from sets import Set
import random
import os

ymax = 491

def getXM():
	f = open('mhad435.10m', 'r')
	map = {}
	for line in f.readlines():
		if line.startswith('%'):
			continue
		y = [float(value) for value in line.split()]
		map[y[0]] = y[1]
		
	f.close()
	return map

def constructXYM(srcx=0):
	map = getXM()

	xylist = []
	f = open('hadsund.dhm')
	for line in f.readlines():
		p = [float(value) for value in line.split()]
		p[1] += 2.4
		if abs(p[0]-srcx) in map:
			p.append(map[abs(p[0]-srcx)])
			xylist.append(p)

	return xylist

def getRandomPairs(k, srcx=0):
	li = constructXYM(srcx)
	tmpSet = Set()
	res = []
	for i in range(2*k):
		while True:
			idx = random.randrange(len(li))			
			if idx not in tmpSet:
				res.append(li[idx])
				tmpSet.add(idx)
				break
	writeToFile(res)

def writeToFile(li, out='out.txt'):
	print len(li)
	f = open(out, 'w')
	for i in li:
		f.write(' '.join([str(j) for j in i]))
		f.write('\n')
	f.close()
