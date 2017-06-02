import sys, csv

try:
	pslout = open(sys.argv[1])
	verification = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	print('usage: script.py pslparsedoutputfile verificationset outfilename')
	sys.exit(1)

vset = {}  # verification set
vsetneg = {}  # negative interaction
assessed = {} # assessed transcripts
for line in verification:
	line = line.rstrip().split('\t')
	if line[8] == 'positive':
		vset[(line[3], line[6])] = line[1]
	else:
		vsetneg[(line[3], line[6])] = line[1]
	if line[3] not in assessed:
		assessed[line[3]] = set()
	assessed[line[3]].add(line[1])

stats = [0, 0, 0]
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in pslout:
		line = line.rstrip().split('\t')
		if line[0] == vset[(line[2],line[3])]:
			writer.writerow(line + ['1'])  # good
			stats[0] += 1
		elif line[0] not in assessed[line[2]]:
		# else:
			writer.writerow(line + ['2'])  # ambiguous
			stats[1] += 1
		elif (line[2],line[3]) in vsetneg and line[0] == vsetneg[(line[2],line[3])]:
			writer.writerow(line + ['0'])  # wrong
			stats[2] += 1
		else:
			writer.writerow(line + ['2'])
			stats[1] += 1

print(stats, 'right', 'ambiguous', 'wrong')
