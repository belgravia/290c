import sys, csv

try:
	toil = open(sys.argv[1])  # toil data
	mapfile = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	print('usage: script.py toil_matrix mapfile outfilename')
	sys.exit(1)

transcripts = set()
transcripts_leftout = set()
for line in mapfile:
	line = line.rstrip().split('\t')
	transcripts.add(line[1])
	transcripts_leftout.add(line[1])

num=0
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	header = toil.readline().rstrip().split('\t')
	header = [h.replace('"', '') for h in header]
	writer.writerow(header)
	#writer.writerow(toil.readline().rstrip().split('\t'))  # header
	for line in toil:
		line = line.rstrip().split()
		line[0] = line[0].replace('"', '')
		if line[0] in transcripts:
			writer.writerow(line)
			num += 1
			transcripts_leftout.remove(line[0])
print(num)
sys.stderr.write(','.join(transcripts_leftout))

