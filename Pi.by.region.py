import sys
inv_id = sys.argv[1]

boundary=20000
winsize =2000

path1 = "../INV_v_expression/"

out1 = open(inv_id+".pairwise.pi.genomic_windows.txt","w")
out2 = open(inv_id+".AA_AD_DD.pi.genomic_windows.txt","w")

crosslist=["62","155","444","502","541","664","909","1034","1192"]

polarity={}
src  =open("revised_inv_list.v3.txt","r")
for line_idz, liner in enumerate(src):
	cols = liner.replace('\n', '').split('\t')
# INV_ID	767orientation	Outgroup seg	alts	chrom	brk1	brk2	genes inside	b1left	b1right	b2left	b2right
# 1	Ancestral	No	444	Chr_01	2448600	2501325	11	2448600	2448600	2501325	2501325

	if cols[0]==inv_id:
		out2.write(cols[0]+'\t'+cols[1]+'\t'+cols[2]+'\t'+cols[3]+'\t'+cols[4]+'\t'+cols[5]+'\t'+cols[6]+'\n')
		out2.write('chrom\twindow_center\tnAA\tpi_AA\tnAD\tpi_AD\tnDD\tpi_DD\n')
		alts = cols[3].split(",")
		chrom= cols[4]
		stp= int( cols[5] )-boundary
		endp=int( cols[6] )+boundary
		AD767=cols[1]
		polarity["767"]=AD767
		for ln in crosslist:
			if AD767=="Ancestral":
				polarity[ln]="Ancestral"
				if ln in alts:
					polarity[ln]="Derived"
			elif AD767=="Derived":
				polarity[ln]="Derived"
				if ln in alts:
					polarity[ln]="Ancestral"
			else:
				polarity[ln]="Unclear"

src.close()

# print "polarity",polarity
seqlen=endp-stp

pi_pw={}
pi_767={}
pi_AA={}
pi_AD={}
pi_DD={}
for w in range(int(seqlen/winsize)+1):
	pi_pw[w]={}
	pi_AA[w]=[0,0]
	pi_AD[w]=[0,0]
	pi_DD[w]=[0,0]
	for j in range(8):
		pi_pw[w][j]={}
		for k in range(j+1,9):
			pi_pw[w][j][k]=[0,0]
	pi_767[w]={}
	for j in range(9):
		pi_767[w][j]=[0,0]

# print stp,endp,"LAST WIN",w

src  =open(path1+chrom+".calls.txt", "r")
for line_idx, line in enumerate(src):
	cols = line.replace('\n', '').split('\t')
	if line_idx==0:
# chrom	pos	62	155	444	502	541	664	909	1034	1192
		pass
	elif int(cols[1])>= stp and int(cols[1])<= endp:
# Chr_01	1	1	0	0	0	0	0	0	0	0
		relpos=int(cols[1])-stp
		w = int(relpos/winsize)

		for j in range(8):
			if int(cols[j+2])>0:
				pi_767[w][j][0]+=1
				if int(cols[j+2])>1:
					pi_767[w][j][1]+=1
				for k in range(j+1,9):
					if int(cols[k+2])>0:
						pi_pw[w][j][k][0]+=1
						if int(cols[j+2]) != int(cols[k+2]):
							pi_pw[w][j][k][1]+=1

		if int(cols[8+2])>0: # v 767 for last line
			pi_767[w][8][0]+=1
			if int(cols[8+2])>1:
				pi_767[w][8][1]+=1
src.close()

for w in pi_pw:
	pos = int( stp + (w+0.5)*winsize )
	for j in range(9):
		if pi_767[w][j][0]>0:
			rx=float(pi_767[w][j][1])/float(pi_767[w][j][0])
			out1.write(chrom+'\t'+str(pos)+'\t'+crosslist[j]+'\t767\t'+str(pi_767[w][j][0])+'\t'+str(rx)+'\n')
			if AD767=="Ancestral":
				if polarity[crosslist[j]]=="Ancestral":
					pi_AA[w][0]+=pi_767[w][j][0]
					pi_AA[w][1]+=pi_767[w][j][1]
				else:
					pi_AD[w][0]+=pi_767[w][j][0]
					pi_AD[w][1]+=pi_767[w][j][1]
			elif AD767=="Derived":
				if polarity[crosslist[j]]=="Ancestral":
					pi_AD[w][0]+=pi_767[w][j][0]
					pi_AD[w][1]+=pi_767[w][j][1]
				else:
					pi_DD[w][0]+=pi_767[w][j][0]
					pi_DD[w][1]+=pi_767[w][j][1]
			elif AD767=="Unclear":
				pi_AA[w][0]+=pi_767[w][j][0] # AA is all contrasts
				pi_AA[w][1]+=pi_767[w][j][1]



		else:
			out1.write(chrom+'\t'+str(pos)+'\t'+crosslist[j]+'\t767\t'+str(pi_767[w][j][0])+'\tNA\n')

	for j in range(8):
		for k in range(j+1,9):
			if pi_pw[w][j][k][0]>0:
				rx=float(pi_pw[w][j][k][1])/float(pi_pw[w][j][k][0])
				out1.write(chrom+'\t'+str(pos)+'\t'+crosslist[j]+'\t'+crosslist[k]+'\t'+str(pi_pw[w][j][k][0])+'\t'+str(rx)+'\n')
				if polarity[crosslist[j]]=="Ancestral" and polarity[crosslist[k]]=="Ancestral":
					pi_AA[w][0]+=pi_pw[w][j][k][0]
					pi_AA[w][1]+=pi_pw[w][j][k][1]
				elif polarity[crosslist[j]]=="Derived" and polarity[crosslist[k]]=="Derived":
					pi_DD[w][0]+=pi_pw[w][j][k][0]
					pi_DD[w][1]+=pi_pw[w][j][k][1]
				elif polarity[crosslist[j]]=="Ancestral" or polarity[crosslist[j]]=="Derived": # not unclear must be AD contrast
					pi_AD[w][0]+=pi_pw[w][j][k][0]
					pi_AD[w][1]+=pi_pw[w][j][k][1]
				elif polarity[crosslist[j]]=="Unclear" and polarity[crosslist[k]]=="Unclear":
					pi_AA[w][0]+=pi_pw[w][j][k][0]
					pi_AA[w][1]+=pi_pw[w][j][k][1]
				else:
					print "error 76",j,k,polarity

			else:
				out1.write(chrom+'\t'+str(pos)+'\t'+crosslist[j]+'\t'+crosslist[k]+'\t'+str(pi_pw[w][j][k][0])+'\tNA\n')



for w in pi_pw:
	pos = int( stp + (w+0.5)*winsize )
	out2.write(chrom+'\t'+str(pos))
	if pi_AA[w][0]>0:
		rx=float(pi_AA[w][1])/float(pi_AA[w][0])
		out2.write('\t'+str(pi_AA[w][0])+'\t'+str(rx))
	else:
		out2.write('\t'+str(pi_AA[w][0])+'\tNA')
	if pi_AD[w][0]>0:
		rx=float(pi_AD[w][1])/float(pi_AD[w][0])
		out2.write('\t'+str(pi_AD[w][0])+'\t'+str(rx))
	else:
		out2.write('\t'+str(pi_AD[w][0])+'\tNA')
	if pi_DD[w][0]>0:
		rx=float(pi_DD[w][1])/float(pi_DD[w][0])
		out2.write('\t'+str(pi_DD[w][0])+'\t'+str(rx))
	else:
		out2.write('\t'+str(pi_DD[w][0])+'\tNA')
	out2.write('\n')

