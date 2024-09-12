# For each inversion, we find a gene scored in the PLoS Genetics paper (cis.geno files) and output the calls
# pull genotypes from 12987 [gene_id].cis.txt files
import sys

out1 = open("inv.gg.txt","w") # report of inversion associated genes

genotyped_genes={}
in11=open("genotyped.genes.txt","r")
for line_idx, line in enumerate(in11):
	cols = line.replace('\n', '').split('\t') 
# MiIM7v11010609m.g
	genotyped_genes[cols[0]]=1
in11.close()


inv_info={}
inxs=open("Genes.within_flanking.inversions","r")
for line_idz, liner in enumerate(inxs):
	cols = liner.replace('\n', '').split('\t')

# Revsied_ms_id	type	line	chrom	brk1	brk2	Chrom	stpos	endpos	old_name	new_name	62	155	444	502	541	664	909	1034	1192	scored_pops
# 56	within	664,444	Chr_13	30242943	30247935	Chr_13	30243281	30244686	MiIM7v11036762m.g	MgIM767.13G191900	yes	yes	no	no	yes	no	yesyes	yes	6

	if line_idz>0:
		try:
			uk=inv_info[cols[0]]
		except KeyError:
			inv_info[cols[0]]={"within":[],"left_flank":[],"right_flank":[]}
		try:
			uz=genotyped_genes[cols[9]]
			inv_info[cols[0]][cols[1]].append(cols[9])
		except KeyError:
			pass

inxs.close()

for inv in inv_info:
	out1.write(inv+'\twithin\t'+str(len(inv_info[inv]["within"]))+'\t'+str( inv_info[inv]["within"] )+'\n')
	out1.write(inv+'\tleft_flank\t'+str(len(inv_info[inv]["left_flank"]))+'\t'+str( inv_info[inv]["left_flank"] )+'\n')
	out1.write(inv+'\tright_flank\t'+str(len(inv_info[inv]["right_flank"]))+'\t'+str( inv_info[inv]["right_flank"] )+'\n')
