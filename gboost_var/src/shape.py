import re
import sys

argvs = sys.argv
argc = len(argvs)

if (argc != 2):
	print 'Usage: # python %s filename' % argvs[0]
	quit()

output = argvs[1] + '_reg'
f = open(output,'w')

for line in open(argvs[1], 'r'):
	tmp = line.split(' ')
	backward = []
	dellist = []
	#print tmp
	for i in range(len(tmp)):
		if ('b' in tmp[i]):
			backward.append(tmp[i - 1].strip())
			backward.append(tmp[i].strip())
			dellist.append(i - 1)
			dellist.append(i)

	dif = 0
	for i in dellist:
		del tmp[i - dif]
		dif += 1

	bparts = []
	bedge = []
	for i in range(len(backward)):
		if (i % 2 == 1):
			bparts.append(re.findall(r'[0-9]+', backward[i]))
		else:
			bedge.append(backward[i])
	bparts.append(['100000' , '100000'])
	bedge.append('100000')

	#print backward
	#print dellist
	#print bparts
	#print bedge
	
	num = 0
	insert = int(bparts[num][0])
	ans = []
	for i in range(len(tmp)):
		ans.append(tmp[i].strip()) 
		if (i == (insert * 2)):
			#print i
			ans.append(bedge[num])
			ans.append('(b' + bparts[num][1] + ')')
			num += 1
	ans.append('\n')


	printline = ""
	
	for i in ans:
		printline += i + " "
		#print i,

	f.write(printline.strip())
	f.write('\n')
