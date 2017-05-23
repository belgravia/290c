import sys, csv

all_genes = []
all_interactions = []
for line in open(sys.argv[1]):
    line = line.rstrip().split("\t")
    #for gene in line:
     #   if gene not in all_genes:
      #      all_genes.append(gene)
    all_genes += [line[0]]

#with open("transcript_obs.txt", 'w') as f:
 #   for gene in all_genes:
  #      print(gene,file=f)


with open(sys.argv[2], 'wt') as f:
	writer = csv.writer(f, delimiter='\t')
	for gene1 in all_genes:
		for gene2 in all_genes:
			if (gene1, gene2) not in all_interactions:
				#print(gene1+"\t"+gene2)
				writer.writerow([gene1, gene2])
