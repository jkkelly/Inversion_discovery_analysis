# identify genes in inversions
import sys

out1=open("Genes.within_flanking.inversions","w")

byinv={}
in1=open("revised_inv_list.txt","r")
for line_idx, line in enumerate(in1):
	cols = line.replace('\n', '').split('\t')
# Revsied_ms_id	line	chrom	brk1	brk2	
# 1	444	Chr_01	2448600	2501325
	if line_idx==0:
		sx_head=line.replace('\n', '')
	else:
		try:
			byinv[cols[0]][0]+=","+cols[1]

		except KeyError:
			byinv[cols[0]]=[cols[1],cols[2],int(cols[3]),int(cols[4])]
in1.close()


left={}
right={}
for inv in byinv:
	chrom=byinv[inv][1]
	stpos=byinv[inv][2]
	endpos=byinv[inv][3]
	left[inv]=[-10000000,'']
	right[inv]=[10000000,'']
	headtext=inv+'\twithin\t'+byinv[inv][0]+'\t'+byinv[inv][1]+'\t'+str(byinv[inv][2])+'\t'+str(byinv[inv][3])

	ct=0
	in2=open("Genes_scored_by_cross.txt","r")
	for line_id, liner in enumerate(in2):
		colx = liner.replace('\n', '').split('\t') 
# Chrom	stpos	endpos	old_name	new_name	62	155	444	502	541	664	909	1034	1192	scored_pops
# Chr_11	3119244	3123078	MiIM7v11029256m.g	MgIM767.11G057200	yes	yes	yes	yes	yes	yes	yes	no	yes	8
		if line_id==0:
			sy_header=liner

		elif line_id>0:
			cno= colx[0]
			if cno==chrom:
				if int(colx[1])>=stpos and int(colx[1])<=endpos:
					out1.write(headtext+'\t'+liner)
					ct+=1
				elif int(colx[2])>=stpos and int(colx[2])<=endpos:
					out1.write(headtext+'\t'+liner)
					ct+=1
				elif int(colx[2])-stpos < 0 and (int(colx[2])-stpos)>left[inv][0]:
					left[inv]=[int(colx[2])-stpos,liner]

				elif int(colx[2])-endpos > 0 and (int(colx[2])-endpos)<right[inv][0]:
					right[inv]=[int(colx[2])-endpos,liner]

	in2.close()
	out1.write(inv+'\tleft_flank\t'+byinv[inv][0]+'\t'+byinv[inv][1]+'\t'+str(byinv[inv][2])+'\t'+str(byinv[inv][3])+'\t'+left[inv][1])
	out1.write(inv+'\tright_flank\t'+byinv[inv][0]+'\t'+byinv[inv][1]+'\t'+str(byinv[inv][2])+'\t'+str(byinv[inv][3])+'\t'+right[inv][1])
	print (headtext,ct)

in1.close()


out1.write(sx_head+'\t'+sy_header)


