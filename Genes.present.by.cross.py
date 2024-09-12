# determine the IM767 genes that were present in the alt lines (captured by Liftoff)
import sys

out1=open("Genes_scored_by_cross.txt","w")

data={}
in1=open("IM767.gene.names.positions.txt","r") # 
for line_id, line in enumerate(in1):
	cols = line.replace('\n', '').split('\t') 
# Chrom	stpos	endpos	old_name	new_name
# Chr_01	13982	16715	MiIM7v11000002m.g	MgIM767.01G000100
	data[cols[3]]={}
	data[cols[3]]['info']=line.replace('\n', '')
	for lins in ['62','155','444','502','541','664','909','1034','1192']:
		data[cols[3]][lins]="no"

in1.close()


for lins in ['62','155','444','502','541','664','909','1034','1192']:

	in2=open("transcript."+lins+".filelist.txt","r") # this is the read counts mapped to each allele of each gene (expression data)
# gene_id	767_s1_767-P1	909_s1_767-P1	767_s1_767-P2	909_s1_767-P2	767_s1_767-P6	909_s1_767-P6	
# MiIM7v11006787m.g	135.0	0.0	165.113	11.886
	for line_id, liner in enumerate(in2):
		colx = liner.replace('\n', '').split('\t') 
		if line_id>0:
			data[colx[0]][lins]="yes"
	in2.close()


out1.write("Chrom\tstpos\tendpos\told_name\tnew_name")
for lins in ['62','155','444','502','541','664','909','1034','1192']:
	out1.write('\t'+lins)
out1.write('\tscored_pops\n')

for gene in data:
	out1.write(data[gene]['info'])
	cc=0
	for lins in ['62','155','444','502','541','664','909','1034','1192']:
		out1.write('\t'+data[gene][lins])
		if data[gene][lins]=="yes":
			cc+=1
	out1.write('\t'+str(cc)+'\n')
