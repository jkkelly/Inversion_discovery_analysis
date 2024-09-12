# output pi(RR), pi(RI), pi(II) for each gene and each inversion # I = inverted relative to 767

# line_list = ["62","155","444","502","541","664","909","1034","1192"]

# out1=open("annotated_inv_list.txt","w")
out2=open("pi.byinv.genic.AD.txt","w")
out3=open("pi.byinv.genic.unpolarized.txt","w")

inv_info={}
AD={}
inv_pi={}
invflank_pi={}
interior_pi={}
inside_pi={}
outside_pi={}
altcc={}
inxs=open("revised_inv_list.v3.txt","r")
for line_idz, liner in enumerate(inxs):
	colx = liner.replace('\n', '').split('\t') 

# INV_ID	767orientation	Outgroup seg	alts	chrom	brk1	brk2	genes inside	b1left	b1right	b2left	b2right
# 1	Ancestral	No	444	Chr_01	2448600	2501325	11	2448600	2448600	2501325	2501325

	if line_idz>0:
		inv_id = colx[0]
		inv_info[inv_id]=[ colx[1],colx[2],colx[4],"a_"+colx[3], str(int(colx[6])-int(colx[5])) ] # ,int(colx[4]),int(colx[5]),[int(colx[7]),int(colx[8])],[int(colx[9]),int(colx[10])] ]
		alts = colx[3].split(",")
		altcc[inv_id]=len(alts)
		if colx[1]=="Ancestral" or colx[1]=="Derived":

			inv_pi[inv_id]={"AA":[0,0],"AD":[0,0],"DD":[0,0]}
			invflank_pi[inv_id]={"AA":[0,0],"AD":[0,0],"DD":[0,0]}
			inside_pi[inv_id]={"AA":[0,0],"AD":[0,0],"DD":[0,0]}
			interior_pi[inv_id]={"AA":[0,0],"AD":[0,0],"DD":[0,0]}
			outside_pi[inv_id]={"AA":[0,0],"AD":[0,0],"DD":[0,0]}
		else:
			inv_pi[inv_id]={"RR":[0,0],"RI":[0,0],"II":[0,0]}
			invflank_pi[inv_id]={"RR":[0,0],"RI":[0,0],"II":[0,0]}
			inside_pi[inv_id]={"RR":[0,0],"RI":[0,0],"II":[0,0]}
			interior_pi[inv_id]={"RR":[0,0],"RI":[0,0],"II":[0,0]}
			outside_pi[inv_id]={"RR":[0,0],"RI":[0,0],"II":[0,0]}
inxs.close()

in11=open("pi.byinv.genic.detailed.txt","r")
for line_idx, line in enumerate(in11):
	cols = line.replace('\n', '').split('\t') 
# INV_ID	chrom	pos1	pos2	alt_lines	gene	invpos	invpos_detailed	RRn	RRpi	RIn	RIpi	IIn	IIpi
# 1	Chr_01	2437546	2439679	a_444	MiIM7v11000759m.g	left_flank	left_flank	76044	0.0286939140498	19086	0.0325369380698	0	NA
# 1	Chr_01	2450545	2459992	a_444	MiIM7v11000763m.g	within	within_1	337112	0.015377678635	82697	0.0261920021282	0	NA
	if line_idx>0:
		inv_id=cols[0]
		ginv_loc = cols[6]
		loc_detail = cols[7]
		gstart = int(cols[2])
		gend = int(cols[3])
		AD767 = inv_info[inv_id][0]

		RRn=float(cols[8])
		RRa=0
		if RRn>0:
			RRa=float(cols[9])*RRn
		RIn=float(cols[10])
		RIa=0
		if RIn>0:
			RIa=float(cols[11])*RIn
		IIn=float(cols[12])
		IIa=0
		if IIn>0:
			IIa=float(cols[13])*IIn

		if ginv_loc=="within":

			if AD767=="Ancestral":
				inv_pi[inv_id]["AA"][0]+=RRn
				inv_pi[inv_id]["AA"][1]+=RRa
				inv_pi[inv_id]["AD"][0]+=RIn
				inv_pi[inv_id]["AD"][1]+=RIa
				inv_pi[inv_id]["DD"][0]+=IIn
				inv_pi[inv_id]["DD"][1]+=IIa

			elif AD767=="Derived":
				inv_pi[inv_id]["DD"][0]+=RRn
				inv_pi[inv_id]["DD"][1]+=RRa
				inv_pi[inv_id]["AD"][0]+=RIn
				inv_pi[inv_id]["AD"][1]+=RIa
				inv_pi[inv_id]["AA"][0]+=IIn
				inv_pi[inv_id]["AA"][1]+=IIa

			elif AD767=="Unclear":
				inv_pi[inv_id]["RR"][0]+=RRn
				inv_pi[inv_id]["RR"][1]+=RRa
				inv_pi[inv_id]["RI"][0]+=RIn
				inv_pi[inv_id]["RI"][1]+=RIa
				inv_pi[inv_id]["II"][0]+=IIn
				inv_pi[inv_id]["II"][1]+=IIa
			else:
				print "fail 43",AD767

			if loc_detail=="within":
				if AD767=="Ancestral":
					interior_pi[inv_id]["AA"][0]+=RRn
					interior_pi[inv_id]["AA"][1]+=RRa
					interior_pi[inv_id]["AD"][0]+=RIn
					interior_pi[inv_id]["AD"][1]+=RIa
					interior_pi[inv_id]["DD"][0]+=IIn
					interior_pi[inv_id]["DD"][1]+=IIa

				elif AD767=="Derived":
					interior_pi[inv_id]["DD"][0]+=RRn
					interior_pi[inv_id]["DD"][1]+=RRa
					interior_pi[inv_id]["AD"][0]+=RIn
					interior_pi[inv_id]["AD"][1]+=RIa
					interior_pi[inv_id]["AA"][0]+=IIn
					interior_pi[inv_id]["AA"][1]+=IIa

				elif AD767=="Unclear":
					interior_pi[inv_id]["RR"][0]+=RRn
					interior_pi[inv_id]["RR"][1]+=RRa
					interior_pi[inv_id]["RI"][0]+=RIn
					interior_pi[inv_id]["RI"][1]+=RIa
					interior_pi[inv_id]["II"][0]+=IIn
					interior_pi[inv_id]["II"][1]+=IIa
				else:
					print "fail interior",AD767



		if loc_detail=="left_flank" or loc_detail=="right_flank" or loc_detail=="within_1" or loc_detail=="within_2":  

			if AD767=="Ancestral":
				invflank_pi[inv_id]["AA"][0]+=RRn
				invflank_pi[inv_id]["AA"][1]+=RRa
				invflank_pi[inv_id]["AD"][0]+=RIn
				invflank_pi[inv_id]["AD"][1]+=RIa
				invflank_pi[inv_id]["DD"][0]+=IIn
				invflank_pi[inv_id]["DD"][1]+=IIa

			elif AD767=="Derived":
				invflank_pi[inv_id]["DD"][0]+=RRn
				invflank_pi[inv_id]["DD"][1]+=RRa
				invflank_pi[inv_id]["AD"][0]+=RIn
				invflank_pi[inv_id]["AD"][1]+=RIa
				invflank_pi[inv_id]["AA"][0]+=IIn
				invflank_pi[inv_id]["AA"][1]+=IIa

			elif AD767=="Unclear":
				invflank_pi[inv_id]["RR"][0]+=RRn
				invflank_pi[inv_id]["RR"][1]+=RRa
				invflank_pi[inv_id]["RI"][0]+=RIn
				invflank_pi[inv_id]["RI"][1]+=RIa
				invflank_pi[inv_id]["II"][0]+=IIn
				invflank_pi[inv_id]["II"][1]+=IIa

			else:
				print "fail 43",AD767

			if loc_detail=="within_1" or loc_detail=="within_2": 

				if AD767=="Ancestral":
					inside_pi[inv_id]["AA"][0]+=RRn
					inside_pi[inv_id]["AA"][1]+=RRa
					inside_pi[inv_id]["AD"][0]+=RIn
					inside_pi[inv_id]["AD"][1]+=RIa
					inside_pi[inv_id]["DD"][0]+=IIn
					inside_pi[inv_id]["DD"][1]+=IIa

				elif AD767=="Derived":
					inside_pi[inv_id]["DD"][0]+=RRn
					inside_pi[inv_id]["DD"][1]+=RRa
					inside_pi[inv_id]["AD"][0]+=RIn
					inside_pi[inv_id]["AD"][1]+=RIa
					inside_pi[inv_id]["AA"][0]+=IIn
					inside_pi[inv_id]["AA"][1]+=IIa

				elif AD767=="Unclear":
					inside_pi[inv_id]["RR"][0]+=RRn
					inside_pi[inv_id]["RR"][1]+=RRa
					inside_pi[inv_id]["RI"][0]+=RIn
					inside_pi[inv_id]["RI"][1]+=RIa
					inside_pi[inv_id]["II"][0]+=IIn
					inside_pi[inv_id]["II"][1]+=IIa

			else:

				if AD767=="Ancestral":
					outside_pi[inv_id]["AA"][0]+=RRn
					outside_pi[inv_id]["AA"][1]+=RRa
					outside_pi[inv_id]["AD"][0]+=RIn
					outside_pi[inv_id]["AD"][1]+=RIa
					outside_pi[inv_id]["DD"][0]+=IIn
					outside_pi[inv_id]["DD"][1]+=IIa

				elif AD767=="Derived":
					outside_pi[inv_id]["DD"][0]+=RRn
					outside_pi[inv_id]["DD"][1]+=RRa
					outside_pi[inv_id]["AD"][0]+=RIn
					outside_pi[inv_id]["AD"][1]+=RIa
					outside_pi[inv_id]["AA"][0]+=IIn
					outside_pi[inv_id]["AA"][1]+=IIa

				elif AD767=="Unclear":
					outside_pi[inv_id]["RR"][0]+=RRn
					outside_pi[inv_id]["RR"][1]+=RRa
					outside_pi[inv_id]["RI"][0]+=RIn
					outside_pi[inv_id]["RI"][1]+=RIa
					outside_pi[inv_id]["II"][0]+=IIn
					outside_pi[inv_id]["II"][1]+=IIa


in11.close()


for inv_id in inv_info:
	q_derived=-99
	if inv_info[inv_id][0]=="Ancestral":
		q_derived = float(altcc[inv_id])/10.0
	elif inv_info[inv_id][0]=="Derived":
		q_derived = 1.0 - float(altcc[inv_id])/10.0
	print (inv_id+'\t'+inv_info[inv_id][0]+'\t'+inv_info[inv_id][1],q_derived)

	if inv_info[inv_id][0]=="Ancestral" or inv_info[inv_id][0]=="Derived":

		out2.write( inv_id+'\t'+inv_info[inv_id][0]+'\t'+inv_info[inv_id][1]+'\t'+inv_info[inv_id][2]+'\t'+inv_info[inv_id][3]+'\t'+inv_info[inv_id][4]+'\t'+str(q_derived) )

		out2.write('\t'+str(inv_pi[inv_id]["AA"][0])+'\t'+str(inv_pi[inv_id]["AA"][1])+'\t'+str(inv_pi[inv_id]["AD"][0])+'\t'+str(inv_pi[inv_id]["AD"][1])+'\t'+str(inv_pi[inv_id]["DD"][0])+'\t'+str(inv_pi[inv_id]["DD"][1]))
		out2.write('\t'+str(invflank_pi[inv_id]["AA"][0])+'\t'+str(invflank_pi[inv_id]["AA"][1])+'\t'+str(invflank_pi[inv_id]["AD"][0])+'\t'+str(invflank_pi[inv_id]["AD"][1])+'\t'+str(invflank_pi[inv_id]["DD"][0])+'\t'+str(invflank_pi[inv_id]["DD"][1]))
		out2.write('\t'+str(interior_pi[inv_id]["AA"][0])+'\t'+str(interior_pi[inv_id]["AA"][1])+'\t'+str(interior_pi[inv_id]["AD"][0])+'\t'+str(interior_pi[inv_id]["AD"][1])+'\t'+str(interior_pi[inv_id]["DD"][0])+'\t'+str(interior_pi[inv_id]["DD"][1]))
		out2.write('\t'+str(inside_pi[inv_id]["AA"][0])+'\t'+str(inside_pi[inv_id]["AA"][1])+'\t'+str(inside_pi[inv_id]["AD"][0])+'\t'+str(inside_pi[inv_id]["AD"][1])+'\t'+str(inside_pi[inv_id]["DD"][0])+'\t'+str(inside_pi[inv_id]["DD"][1]))
		out2.write('\t'+str(outside_pi[inv_id]["AA"][0])+'\t'+str(outside_pi[inv_id]["AA"][1])+'\t'+str(outside_pi[inv_id]["AD"][0])+'\t'+str(outside_pi[inv_id]["AD"][1])+'\t'+str(outside_pi[inv_id]["DD"][0])+'\t'+str(outside_pi[inv_id]["DD"][1])+'\n')

	else:
		q_alt = float(altcc[inv_id])/10.0
		out3.write( inv_id+'\t'+inv_info[inv_id][0]+'\t'+inv_info[inv_id][1]+'\t'+inv_info[inv_id][2]+'\t'+inv_info[inv_id][3]+'\t'+inv_info[inv_id][4]+'\t'+str(q_alt) )

		out3.write('\t'+str(inv_pi[inv_id]["RR"][0])+'\t'+str(inv_pi[inv_id]["RR"][1])+'\t'+str(inv_pi[inv_id]["RI"][0])+'\t'+str(inv_pi[inv_id]["RI"][1])+'\t'+str(inv_pi[inv_id]["II"][0])+'\t'+str(inv_pi[inv_id]["II"][1]))
		out3.write('\t'+str(invflank_pi[inv_id]["RR"][0])+'\t'+str(invflank_pi[inv_id]["RR"][1])+'\t'+str(invflank_pi[inv_id]["RI"][0])+'\t'+str(invflank_pi[inv_id]["RI"][1])+'\t'+str(invflank_pi[inv_id]["II"][0])+'\t'+str(invflank_pi[inv_id]["II"][1]))
		out3.write('\t'+str(interior_pi[inv_id]["RR"][0])+'\t'+str(interior_pi[inv_id]["RR"][1])+'\t'+str(interior_pi[inv_id]["RI"][0])+'\t'+str(interior_pi[inv_id]["RI"][1])+'\t'+str(interior_pi[inv_id]["II"][0])+'\t'+str(interior_pi[inv_id]["II"][1]))
		out3.write('\t'+str(inside_pi[inv_id]["RR"][0])+'\t'+str(inside_pi[inv_id]["RR"][1])+'\t'+str(inside_pi[inv_id]["RI"][0])+'\t'+str(inside_pi[inv_id]["RI"][1])+'\t'+str(inside_pi[inv_id]["II"][0])+'\t'+str(inside_pi[inv_id]["II"][1]))
		out3.write('\t'+str(outside_pi[inv_id]["RR"][0])+'\t'+str(outside_pi[inv_id]["RR"][1])+'\t'+str(outside_pi[inv_id]["RI"][0])+'\t'+str(outside_pi[inv_id]["RI"][1])+'\t'+str(outside_pi[inv_id]["II"][0])+'\t'+str(outside_pi[inv_id]["II"][1])+'\n')



