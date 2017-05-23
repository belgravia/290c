from subprocess import call
import sys, csv

try:
	prefix = sys.argv[1]
	gg = sys.argv[2]

except:
	print('script.py dir gginteractions [gtf] [] []')
	print('generates 7 files needed to run psl')
	print('make sure scripts are in this directory')
	print('also make sure that a template data and psl file exist in this directory')
	print('and that a ref/ and toil_data/ directory exist in the parent directory of this one')
	sys.exit(1)

call(['mkdir', prefix])
# call(['cd', prefix])
# call(['cp', gg, gg[gg.rfind('/')+1:]])  # gene interaction file
call(['cp', 'data_template', prefix+'/'+'data_'+prefix])  # data
call(['cp', 'psl_template', prefix+'/'+'psl_'+prefix])  # psl
call(['cp', 'command.sh', prefix+'/'+'command.sh'])

all_genes = set()
for line in open(gg):
	line = line.rstrip().split('\t')
	all_genes.update(line)

print('Writing observed genes')
with open(prefix+'_genes_obs.txt', 'wt') as outfile:  # gene only file
	writer = csv.writer(outfile, delimiter='\t')
	for gene in all_genes:
		writer.writerow([gene])

print('Mapping genes to transcript names')
print('Number of genes:')
call(['python', 'pathway_genes_to_transcripts.py', gg, \
	'../ref/gencode.v23.annotation.gtf', prefix+'_gene_transcripts.txt']) # gene transcript mapping

print('Writing observed transcripts')
with open(prefix+'_transcripts_obs.txt', 'wt') as outfile:  # transcript only file
	writer = csv.writer(outfile, delimiter='\t')
	for line in open(prefix+'_gene_transcripts.txt'):
		line = line.rstrip().split('\t')
		writer.writerow([line[1]])

print('Getting expression data')
call(['python', 'pull_toil_given_map.py', '../toil_data/gtex_RSEM_isoform_fpkm_breast', 
	prefix+'_gene_transcripts.txt', prefix+'_transcript_expr.txt'])

print('Normalizing expression data')
transcripts = []
exprvals = []
for line in open(prefix+'_transcript_expr.txt'):
	line = line.rstrip().split('\t')
	expr = [float(n) for n in line[1:]]
	if expr == []:  # the header
		continue
	exprvals += [sum(expr)/len(expr)]
	transcripts += [line[0]]

exprvals = [e + abs(min(exprvals)) for e in exprvals]
exprvals = [e / max(exprvals) for e in exprvals]
exprvals2 = [round(e, 2) for e in exprvals]
exprvals1 = [round(e, 1) for e in exprvals]

with open(prefix+'_transcript_expr_mean2.txt', 'wt') as outfile:  # normalized expression file
	writer = csv.writer(outfile, delimiter='\t')
	for i in range(len(transcripts)):
		writer.writerow([transcripts[i], exprvals2[i]])

with open(prefix+'_transcript_expr_mean1.txt', 'wt') as outfile:  # normalized expression file
	writer = csv.writer(outfile, delimiter='\t')
	for i in range(len(transcripts)):
		writer.writerow([transcripts[i], exprvals1[i]])

print('Writing expression boolean data')
with open(prefix+'_expression_bool.txt', 'wt') as outfile:  # expression bool
	writer = csv.writer(outfile, delimiter='\t')
	for i in range(len(transcripts)):
		if exprvals2[i] != 0:
			writer.writerow([transcripts[i]])

print('Generating pairwise transcript interaction targets')
call(['python', 'dataProcessor.py', prefix+'_transcripts_obs.txt', prefix+'_tt_interactions.txt'])  # target file

print('Moving files to directory ' + prefix)
call(['mv', prefix+'_genes_obs.txt', prefix+'/'+prefix+'_genes_obs.txt'])
call(['mv', prefix+'_gene_transcripts.txt', prefix+'/'+prefix+'_gene_transcripts.txt'])
call(['mv', prefix+'_transcripts_obs.txt', prefix+'/'+prefix+'_transcripts_obs.txt'])
call(['mv', prefix+'_tt_interactions.txt', prefix+'/'+prefix+'_tt_interactions.txt'])
call(['mv', prefix+'_transcript_expr_mean2.txt', prefix+'/'+prefix+'_transcript_expr_mean2.txt'])
call(['mv', prefix+'_transcript_expr_mean1.txt', prefix+'/'+prefix+'_transcript_expr_mean1.txt'])
call(['mv', prefix+'_transcript_expr.txt', prefix+'/'+prefix+'_transcript_expr.txt'])
call(['mv', prefix+'_expression_bool.txt', prefix+'/'+prefix+'_expression_bool.txt'])




