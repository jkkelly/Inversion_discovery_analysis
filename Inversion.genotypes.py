# For each inversion, we find a gene scored in the PLoS Genetics paper (cis.geno files) and output the calls
# pull genotypes from 12987 [gene_id].cis.txt files
import sys

data={}
in11=open("plants.ordered.txt","r")
for line_idx, line in enumerate(in11):
	cols = line.replace('\n', '').split('\t') 
# 1	s5_767-P36	767	three	parent
	data[int(cols[0])-1]=line.replace('\n', '')
in11.close()


inv_info={}
inxs=open("Marker.gene.per.inversion.txt","r")
for line_idz, liner in enumerate(inxs):
	cols = liner.replace('\n', '').split('\t') 
# 24	MiIM7v11014267m.g
	inv_info[cols[0]]=cols[1]
inxs.close()


for inv_id in inv_info:
	print inv_id,inv_info[inv_id]
	out1=open(inv_id+".genotype.txt","w")
	gene_id=inv_info[inv_id]
	in1=open(gene_id+".cis.txt","r") # genotype files produced for eGWAS paper, Veltsos and Kelly (2024; PLoS Genetics)
	for line_idx, line in enumerate(in1):
		cols = line.replace('\n', '').split('\t') 
# 1	0	1	0	0	0	0	0	0	0	0	0	0
		gc=0
		for j in range(4,len(cols)):
			gc+=int(cols[j])
		out1.write(data[line_idx]+'\t'+str(gc)+'\n')

	in1.close()
	out1.close()

