# test each gene in each cross
import sys
from sklearn.linear_model import LinearRegression
import numpy as np

MinMeanCPM=0.5
filesub = sys.argv[1]

out2 = open("Gene_by_cross.tests."+filesub+".txt","w")

reg = LinearRegression() # Create an instance of the LinearRegression class

cross_for_plant={}
in11=open("plants.ordered.txt","r")
for line_idx, line in enumerate(in11):
	cols = line.replace('\n', '').split('\t') 
# 1	s5_767-P36	767	three	parent
	cross_for_plant[cols[1]]=cols[2]
in11.close()

rpp={}
src  =open("raw.reads.per.plant.txt", "r") # use to determine the IM line samples
for line_idx, line in enumerate(src):
	cols = line.replace('\n', '').split('\t')
# 62	s3_62-168	1776928.997
	rpp[cols[1]]=float(cols[2])

src.close()

Lines_in_seq=["62","155","444","502","541","664","909","1034","1192"]
inx=open(filesub,"r") # filesub is a subset of "Genes_scored_by_cross.txt"
for line_id, liner in enumerate(inx):
	colx = liner.replace('\n', '').split('\t')
# Chrom	stpos	endpos	old_name	new_name	62	155	444	502	541	664	909	1034	1192	scored_pops
# Chr_11	3119244	3123078	MiIM7v11029256m.g	MgIM767.11G057200	yes	yes	yes	yes	yes	yes	yes	no	yes	8
	if line_id>0:
		rlines=[]
		old_name=colx[3]
		for j in range(9):
			if colx[j+5]=="yes":
				rlines.append(Lines_in_seq[j])
		if len(rlines)>0: # data from at least one cross
			chrom=colx[0]
			pos=int(colx[1])

			ally=[]
			for cross in rlines: # test if average expression of plants is >= MinMeanCPM
				relevant_numbers={}
				src  =open("transcript."+cross+".filelist.txt", "r")
				for line_idx, line in enumerate(src):
					cols = line.replace('\n', '').split('\t')
					if line_idx==0:
				# gene_id	767_s1_767-P1	62_s1_767-P1	767_s1_767-P2	62_s1_767-P2	...
						for j in range(1,len(cols),2):
							if cols[j][0:2]=="62": # not really needed since only first allele is read
								plt = cols[j][3:]
							else:
								plt = cols[j][4:]

							relevant_numbers[plt]=[j,j+1]

					else:
				# MiIM7v11006787m.g	135.0	0.0	165.113	11.886	16.0	0.0	27.0	0.0	59.0	0.0	3.0	...
						if cols[0]==old_name:
							for plt in relevant_numbers:
								r767=1000000.0*float(cols[relevant_numbers[plt][0]])/rpp[plt]
								rALT=1000000.0*float(cols[relevant_numbers[plt][1]])/rpp[plt]
								ally.append(float(r767+rALT))

				src.close()

			print( old_name,chrom,pos,len(rlines),np.average(ally) )
			if np.average(ally)>=MinMeanCPM: # gene has enough expression

				marker={}
				for lins in Lines_in_seq:
					marker[lins]=[-99,10000000]

				in11=open("marker.for.gene.txt","r")
				for line_idx, line in enumerate(in11):
					cols = line.replace('\n', '').split('\t') 
				# MiIM7v11000007m.g	47990	664	Chr_01	start	47990
					if cols[3]==chrom:
						if cols[4]=="start":
							if abs(int(cols[5])-pos) < marker[cols[2]][1]:
								marker[cols[2]]=[cols[5],abs(int(cols[5])-pos)]

						elif cols[5]=="end":
							if abs(int(cols[4])-pos) < marker[cols[2]][1]:
								marker[cols[2]]=[cols[4],abs(int(cols[4])-pos)]

						else:
							if abs(int(cols[4])-pos) < marker[cols[2]][1]:
								marker[cols[2]]=[cols[4],abs(int(cols[4])-pos)]
							if abs(int(cols[5])-pos) < marker[cols[2]][1]:
								marker[cols[2]]=[cols[5],abs(int(cols[5])-pos)]
				in11.close()


				vals_a={}
				in11=open("all.F2_geno_PP.txt","r")
				for line_idx, line in enumerate(in11):
					cols = line.replace('\n', '').split('\t') 
					# s5_1192-145	Chr_01	47990	BB	BB:1.0
					if cols[1]==chrom:
						if cols[2] == marker[ cross_for_plant[cols[0]] ][0]:
							call=cols[4].split(":")[0]
							vals_a[cols[0]] = 0
							if call=="AB":
								vals_a[cols[0]] = 1
							elif call=="BB":
								vals_a[cols[0]] = 2	
				in11.close()

				for cross in rlines:
					relevant_numbers={}
					src  =open("transcript."+cross+".filelist.txt", "r")
					for line_idx, line in enumerate(src):
						cols = line.replace('\n', '').split('\t')
						if line_idx==0:
					# gene_id	767_s1_767-P1	62_s1_767-P1	767_s1_767-P2	62_s1_767-P2	...
							for j in range(1,len(cols),2):
								if cols[j][0:2]=="62": # not really needed since only first allele is read
									plt = cols[j][3:]
								else:
									plt = cols[j][4:]

								try:
									amigood=vals_a[plt] # must a genotyped f2
									relevant_numbers[plt]=[j,j+1]

								except KeyError:
									pass
						else:
					# MiIM7v11006787m.g	135.0	0.0	165.113	11.886	16.0	0.0	27.0	0.0	59.0	0.0	3.0	...
							if cols[0]==old_name:
								x1=[]
								x2=[]
								y=[]
								y_bc=[]
								datums=[[],[],[]]
								for plt in relevant_numbers:
									r767=1000000.0*float(cols[relevant_numbers[plt][0]])/rpp[plt]
									rALT=1000000.0*float(cols[relevant_numbers[plt][1]])/rpp[plt]
									# out1.write(gene+"\t"+plantx[plt]+"\t"+str(r767)+"\t"+str(rALT)+"\t"+str(r767+rALT)+'\n')
									x1.append([ vals_a[plt] ])
									if vals_a[plt]==0:
										x2.append([0,0])
									elif vals_a[plt]==1:
										x2.append([1,1])
									elif vals_a[plt]==2:
										x2.append([2,0])
									else:
										print "failure 45",genotype[plt]
									y.append(float(r767+rALT))
									# y_bc.append(0.5+float(r767+rALT)) # for boxcox
									datums[ vals_a[plt] ].append(float(r767+rALT))

								out2.write(old_name+"\t"+cross)
								for i in range(3):
									if len(datums[i])>0:
										m1=sum(datums[i])/float(len(datums[i]))
										out2.write("\t"+str(len(datums[i]))+","+str(m1))
									else:
										out2.write("\t0,NA")

								mean_exp=sum(y)/float(len(y))
								if mean_exp>0.0:
									for i in range(len(y)):
										y[i]=y[i]/mean_exp
									Y = np.array(y)
									X = np.array(x1)
									reg.fit(X, Y) # Fit the additive model to the data
									out2.write("\t"+str(reg.coef_[0]))

									X = np.array(x2)
									reg.fit(X, Y) # Fit additive/dominance model to the data
									out2.write("\t"+str(reg.coef_[0])+"\t"+str(reg.coef_[1])+"\n") # Print the coefficients of the model
					src.close()
