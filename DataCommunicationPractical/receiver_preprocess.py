f = open('mhad435.10m', 'r')
x = []
for line in f.readlines():
	if line.startswith('%'):
		continue
	y = [float(value) for value in line.split()]
	x.append(y)
	
f.close()
print x