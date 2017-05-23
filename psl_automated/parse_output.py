import sys

try:
	gimap = open(sys.argv[1])
	outputfile = open(sys.argv[2])
except:
	print('script.py geneisoform.map outputfile [bool:ifyouwantcounts] > parsedoutput')
	sys.exit(1)

gi_map = {}
for line in gimap:
	line = line.rstrip().split('\t')
	gi_map[line[1]] = line[0]

countsbool = len(sys.argv)>3

counts={}
for line in outputfile :
	if '=' not in line or not line.startswith('TT'):
		continue
	line = line.rstrip().split('=')
	num = line[1]
	if num not in counts:
		counts[num] = 0
	counts[num] += 1
	if not countsbool and float(num) > 0.8:
		line0 = line[0].split(',')
		trans0 = line0[0][line0[0].find('ENST'):]
		trans1 = line0[1][1:-2]
		gene0 = gi_map[trans0]
		gene1 = gi_map[trans1]
		print('\t'.join([trans0, trans1, gene0, gene1, num]))
if countsbool:
	print(counts)

