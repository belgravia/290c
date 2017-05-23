import sys, csv

try:
	pathway = open(sys.argv[1])
	gtf = open(sys.argv[2])
	outfilename = sys.argv[3]
except:
	print('usage: script.py pathwayfile gtffile outfilename')
	sys.exit(1)

genes = set()
for line in pathway:
	line = line.rstrip().split('\t')
	genes.update(line[:2])
print(len(genes))
already_written = set()
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t')
	for line in gtf:
		if line.startswith('#'):
			continue
		line = line.rstrip().split('\t')
		tstart = line[8].find('ENST')
		if tstart < 0:
			continue
		transcript = line[8][tstart:]
		transcript = transcript[:transcript.find('"')]
		gene_name = line[8][line[8].find('gene_name')+11:]
		gene_name = gene_name[:gene_name.find('"')]
		if transcript in already_written:
			continue
		if gene_name not in genes:
			continue
		already_written.add(transcript)
		writer.writerow([gene_name, transcript])


