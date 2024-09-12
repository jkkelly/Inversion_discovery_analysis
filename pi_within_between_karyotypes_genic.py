# output pi(RR), pi(RI), pi(II) for each gene and each inversion # I = inverted relative to 767

out2=open("pi.byinv.genic.detailed.txt","w")
out2.write('INV_ID\tchrom\tpos1\tpos2\talt_lines\tgene\tRRn\tRRpi\tRIn\tRIpi\tIIn\tIIpi\n')

inv_info={}
inxs=open("Genes.within_flanking.inversions","r")
for line_idz, liner in enumerate(inxs):
	colx = liner.replace('\n', '').split('\t') 
# Revsied_ms_id	type	line	chrom	brk1	brk2	Chrom	stpos	endpos	old_name	new_name	62	155	444	502	541	664	909	1034	1192	scored_pops
# 56	within	664,444	Chr_13	30242943	30247935	Chr_13	30243281	30244686	MiIM7v11036762m.g	MgIM767.13G191900	yes	yes	no	no	yes	no	yes	yes	yes	6

	if line_idz>0:
		inv_id = colx[0]
		try:
			uk1=inv_info[inv_id]
			inv_info[inv_id][colx[1]].append(colx[9])
		except KeyError:
			chrom = colx[3]
			inv_info[inv_id]={"ch":chrom,"brk1":int(colx[4]),"brk2":int(colx[5]),"alt_lines":colx[2],"left_flank":[],"right_flank":[],"within":[]}
			inv_info[inv_id][colx[1]].append(colx[9])

inxs.close()

for inv_id in inv_info:

	genes=inv_info[inv_id]["within"] #   colx[5].split(",")
	chrom=inv_info[inv_id]["ch"]
	ALTs = inv_info[inv_id]["alt_lines"].split(",")

	bygene_stats_inv={}
	gtype={}
	for g in genes:
		bygene_stats_inv[g]={"RR":[0,0.0],"RI":[0,0.0],"II":[0,0.0]}
		gtype[g]="within"

	if len(inv_info[inv_id]["left_flank"])==1:
		g=inv_info[inv_id]["left_flank"][0]
		bygene_stats_inv[g]={"RR":[0,0.0],"RI":[0,0.0],"II":[0,0.0]}
		gtype[g]="left_flank"

	if len(inv_info[inv_id]["right_flank"])==1:
		g=inv_info[inv_id]["right_flank"][0]
		bygene_stats_inv[g]={"RR":[0,0.0],"RI":[0,0.0],"II":[0,0.0]}
		gtype[g]="right_flank"



	gst={}
	gend={}

	in11=open("pairwise_pi/"+chrom+".gene.pairwise.pi.txt","r")
	for line_idx, line in enumerate(in11):
		cols = line.replace('\n', '').split('\t') 
		# 	Chr_01	2089564	2093370	MiIM7v11000666m.g	62	767	3807	0.0115576569477
		try:
			g=cols[3]
			uk=gtype[g]
			gst[g]=cols[1]
			gend[g]=cols[2]
			if int(cols[6])>0:

				if cols[4] in ALTs and cols[5] in ALTs:
					bygene_stats_inv[g]["II"][0]+=int(cols[6])
					bygene_stats_inv[g]["II"][1]+= float(cols[6])*float(cols[7])
				elif cols[4] in ALTs or cols[5] in ALTs:
					bygene_stats_inv[g]["RI"][0]+=int(cols[6])
					bygene_stats_inv[g]["RI"][1]+= float(cols[6])*float(cols[7])
				else:
					bygene_stats_inv[g]["RR"][0]+=int(cols[6])
					bygene_stats_inv[g]["RR"][1]+= float(cols[6])*float(cols[7])
		except KeyError:
			pass

	in11.close()


	for g in bygene_stats_inv:
		out2.write(inv_id+'\t'+chrom+'\t'+gst[g]+'\t'+gend[g]+'\ta_'+inv_info[inv_id]["alt_lines"]+'\t'+g+'\t'+gtype[g] )
		if bygene_stats_inv[g]["RR"][0]>0:
			rx = bygene_stats_inv[g]["RR"][1]/float(bygene_stats_inv[g]["RR"][0])
			out2.write( '\t'+str(bygene_stats_inv[g]["RR"][0])+'\t'+str(rx) )
		else:
			out2.write( '\t0\tNA' )
		if bygene_stats_inv[g]["RI"][0]>0:
			rx = bygene_stats_inv[g]["RI"][1]/float(bygene_stats_inv[g]["RI"][0])
			out2.write( '\t'+str(bygene_stats_inv[g]["RI"][0])+'\t'+str(rx) )
		else:
			out2.write( '\t0\tNA' )
		if bygene_stats_inv[g]["II"][0]>0:
			rx = bygene_stats_inv[g]["II"][1]/float(bygene_stats_inv[g]["II"][0])
			out2.write( '\t'+str(bygene_stats_inv[g]["II"][0])+'\t'+str(rx) )
		else:
			out2.write( '\t0\tNA' )

		out2.write( '\n' )



