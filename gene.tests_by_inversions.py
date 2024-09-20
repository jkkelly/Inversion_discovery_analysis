import numpy as np
from math import sqrt

Min_genes = 5
out1=open("exp.stats.by.inv.txt","w")
out1a=open("Alpha.by.gene.txt","w")

anc_der={}
in1=open("revised_inv_list.v3.txt","r")
for line_idx, line in enumerate(in1):
	cols = line.replace('\n', '').split('\t') 

# INV_ID	767orientation	Outgroup seg	alts	chrom	brk1	brk2	genes inside	b1left	b1right	b2left	b2right
# 1	Ancestral	No	444	Chr_01	2448600	2501325	11	2448600	2448600	2501325	2501325
	if line_idx>0:
		inv_id = cols[0]
		anc_der[inv_id]=cols[1]
in1.close()

within_inv={}
alpha= {}
alphaNI= {} # NI = not inversion
not767= {}
inx  =open("pi.byinv.genic.detailed.txt", "r")
for line_idx, line in enumerate(inx):
	cols = line.replace('\n', '').split('\t')

# INV_ID	chrom	pos1	pos2	alt_lines	gene	invpos	invpos_detailed	RRn	RRpi	RIn	RIpi	IIn	IIpi
# 1	Chr_01	2437546	2439679	a_444	MiIM7v11000759m.g	left_flank	left_flank	76044	0.0286939140498	19086	0.0325369380698	0	NA
# 1	Chr_01	2450545	2459992	a_444	MiIM7v11000763m.g	within	within_1	337112	0.015377678635	82697	0.0261920021282	0	NA

	if line_idx>0:
		not767[cols[0]]=(cols[4].split("_")[1]).split(",")
		try:
			within_inv[cols[5]].append(cols[0])
		except KeyError:
			within_inv[cols[5]]=[cols[0]]

		alpha[cols[5]+"_"+cols[0]]={"info":cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[6]+'\t'+cols[7]+'\t'+anc_der[cols[0]], "RR":[], "RI":[],"y":[]}
			
inx.close()


inx  =open("Gene_by_cross.tests.v2.txt", "r") # use to determine the IM line samples
for line_id, liner in enumerate(inx):
	cols = liner.replace('\n', '').split('\t')
# MiIM7v11000002m.g	62	37,29.722291373	41,26.7994250422	14,33.8948579473	0.03275253811449811	0.07180544868900296	-0.1724043125362287
	try:
		relinvs = within_inv[ cols[0] ]
		for inv_id in relinvs:
			if cols[1] in not767[ inv_id ]: # RI
				alpha[cols[0]+"_"+inv_id]["RI"].append(float(cols[5])) # raw alpha
			else: # RR
				alpha[cols[0]+"_"+inv_id]["RR"].append(float(cols[5])) # raw alpha
			alpha[cols[0]+"_"+inv_id]["y"].append(float(cols[5])) # raw alpha
	except KeyError:
		try:
			alphaNI[cols[0]]['info']+= ","+cols[5]
			alphaNI[cols[0]]['y'].append( float(cols[5]) )
		except KeyError:
			alphaNI[cols[0]]={}
			alphaNI[cols[0]]['info']="Not_within_inversion\t"+cols[5]
			alphaNI[cols[0]]['y']=[float(cols[5])]

inx.close()

big1=[ [],[],[],[],[],[] ]
lx1=[]
lx2=[]
for g in alpha:
	n0 = len(alpha[g]["RR"])
	n1 = len(alpha[g]["RI"])

	if len(alpha[g]["y"])>=Min_genes and n0>0 and n1>0:

		factoids=alpha[g]["info"].split("\t")
		inv_id=factoids[0]
		posr = factoids[5]
		ad_by_gene = factoids[7]

		out1a.write(g+'\t'+alpha[g]["info"])
		m=np.average(alpha[g]["y"])
		m0=np.average(alpha[g]["RR"])
		m1=np.average(alpha[g]["RI"])
		var_tot = sum((xi - m) ** 2 for xi in alpha[g]["y"]) / float(len(alpha[g]["y"])-1.0)
		lx2.append(var_tot)
		absy=[]
		for xi in alpha[g]["y"]:
			absy.append(abs(xi))
		ma = np.average(absy)
		lx1.append( ma )

		pwd=[-99,-99,-99,-99,-99] # within RR, within RI, between RR/RI, within Anc, within Derived
		pwd_within=[0,0.0]
		if n0>1:
			
			dist=[0,0.0]
			for i in range(n0-1):
				for j in range(i+1,n0): 
					dist[0]+=1
					dist[1]+=abs(alpha[g]["RR"][i]-alpha[g]["RR"][j])
			pwd[0]=dist[1]/float(dist[0])
			pwd_within[0]+=dist[0]
			pwd_within[1]+=dist[1]
		if n1>1:

			dist=[0,0.0]
			for i in range(n1-1):
				for j in range(i+1,n1): 
					dist[0]+=1
					dist[1]+=abs(alpha[g]["RI"][i]-alpha[g]["RI"][j])
			pwd[1]=dist[1]/float(dist[0])
			pwd_within[0]+=dist[0]
			pwd_within[1]+=dist[1]

		if ad_by_gene=="Ancestral":
			pwd[3]=pwd[0]
			pwd[4]=pwd[1]
		elif ad_by_gene=="Derived":
			pwd[3]=pwd[1]
			pwd[4]=pwd[2]

		dist=[0,0.0]
		for i in range(n0):
			for j in range(n1): 
				dist[0]+=1
				dist[1]+=abs(alpha[g]["RR"][i]-alpha[g]["RI"][j])
		pwd[2]=dist[1]/float(dist[0]) # between cross types
		pwd[1]=pwd_within[1]/float(pwd_within[0]) # within cross types

		out1a.write('\t'+str(ma)+'\t'+str(m0)+'\t'+str(m1)+'\t'+str(var_tot)+'\t'+str(pwd[2]))
		out1a.write('\t'+str(n0)+'\t'+str(n1)+'\t'+str(pwd[1])+'\t'+str(pwd[3])+'\t'+str(pwd[4])+"\n")
		
		if posr=="within": # all genes within inversions
			big1[0].append(abs(m1)-abs(m0))
			big1[1].append(pwd[2])
			big1[2].append(pwd[1])
			if pwd[3]>=0.0:
				big1[3].append(pwd[3])
			if pwd[4]>=0.0:
				big1[4].append(pwd[4])

print( "INV genes: mean absolute",np.average(lx1),"Vg cis (raw estimates) ",np.average(lx2) )
print( "means across genes, six stats" )
for j in range(5):
	print( j,len(big1[j]),np.average(big1[j]) )


lx1=[]
lx2=[]
lx3=[]
for g in alphaNI: # results for colinear genome
	if len(alphaNI[g]["y"])>=Min_genes:
		absy=[]
		for xi in alphaNI[g]["y"]:
			absy.append(abs(xi))
		ma = np.average(absy)

		m = sum(alphaNI[g]["y"]) / float(len(alphaNI[g]["y"]))
		var_res = sum((xi - m) ** 2 for xi in alphaNI[g]["y"]) / float(len(alphaNI[g]["y"])-1.0)
		lx1.append(ma)
		lx2.append(var_res)
		lx3.append(sqrt(var_res))
print( len(lx1),"Background genes: mean absolute",np.average(lx1),"Vg cis (raw estimates) ",np.average(lx2)," sqrt vg ",np.average(lx3) )



