# Read Chr_*.calls.txt from read.mummer1.py
import sys
chrom = sys.argv[1]

out1 = open(chrom+".gene.pairwise.pi.txt","w")

pix={}
genelist={}
src  =open("IM767v1.genes.txt", "r")
for line_idx, line in enumerate(src):
	cols = line.replace('\n', '').split('\t')
# Chr_01	13982	16715	+	ID=MiIM7v11000002m.g
	if cols[0]==chrom:
		genelist[cols[4].split("=")[1]]=cols[1]+'\t'+cols[2]
		for j in range(int(cols[1]),int(cols[2])+1):
			pix[j]=cols[4].split("=")[1]
src.close()

crosslist=["62","155","444","502","541","664","909","1034","1192"]
pi_pw={}
pi_767={}
for g in genelist:
	pi_pw[g]={}
	for j in range(8):
		pi_pw[g][j]={}
		for k in range(j+1,9):
			pi_pw[g][j][k]=[0,0]

	pi_767[g]={}
	for j in range(9):
		pi_767[g][j]=[0,0]

src  =open(chrom+".calls.txt", "r")
for line_idx, line in enumerate(src):
	cols = line.replace('\n', '').split('\t')
	if line_idx==0:
# chrom	pos	62	155	444	502	541	664	909	1034	1192
		pass
	else:
# Chr_01	1	1	0	0	0	0	0	0	0	0
		try:
			g=pix[int(cols[1])]
			for j in range(8):
				if int(cols[j+2])>0:
					pi_767[g][j][0]+=1
					if int(cols[j+2])>1:
						pi_767[g][j][1]+=1
					for k in range(j+1,9):
						if int(cols[k+2])>0:
							pi_pw[g][j][k][0]+=1
							if int(cols[j+2]) != int(cols[k+2]):
								pi_pw[g][j][k][1]+=1
			if int(cols[8+2])>0: # v 767 for last line
				pi_767[g][8][0]+=1
				if int(cols[8+2])>1:
					pi_767[g][8][1]+=1

		except KeyError:
			pass

src.close()

for g in genelist:
	for j in range(9):
		if pi_767[g][j][0]>0:
			rx=float(pi_767[g][j][1])/float(pi_767[g][j][0])
			out1.write(chrom+'\t'+genelist[g]+'\t'+g+'\t'+crosslist[j]+'\t767\t'+str(pi_767[g][j][0])+'\t'+str(rx)+'\n')
		else:
			out1.write(chrom+'\t'+genelist[g]+'\t'+g+'\t'+crosslist[j]+'\t767\t'+str(pi_767[g][j][0])+'\tNA\n')

	for j in range(8):
		for k in range(j+1,9):
			if pi_pw[g][j][k][0]>0:
				rx=float(pi_pw[g][j][k][1])/float(pi_pw[g][j][k][0])
				out1.write(chrom+'\t'+genelist[g]+'\t'+g+'\t'+crosslist[j]+'\t'+crosslist[k]+'\t'+str(pi_pw[g][j][k][0])+'\t'+str(rx)+'\n')
			else:
				out1.write(chrom+'\t'+genelist[g]+'\t'+g+'\t'+crosslist[j]+'\t'+crosslist[k]+'\t'+str(pi_pw[g][j][k][0])+'\tNA\n')
