# v1:: read transcript* files from Salmon mapping (transcript*filelist), determine total RNA reads for each plant

out1 = open("raw.reads.per.plant.txt","w")
crosslist=["62","1034","1192","502","541","664","909","155","444"]

for cross in crosslist:
	data={}
	print( "working ",cross )
	relevant_numbers={}
	src  =open("transcript."+cross+".filelist.txt", "r")
	for line_idx, line in enumerate(src):
		cols = line.replace('\n', '').split('\t')
		if line_idx==0:
	# gene_id	767_s1_767-P1	62_s1_767-P1	767_s1_767-P2	62_s1_767-P2	...
			for j in range(1,len(cols),2):
				if cols[j][0:2]=="62":
					plt = cols[j][3:]
				else:
					plt = cols[j][4:]

				relevant_numbers[plt]=[j,j+1]
				data[ plt ]=0.0 # sum of all reads

		else:
	# MiIM7v11006787m.g	135.0	0.0	165.113	11.886	16.0	0.0	27.0	0.0	59.0	0.0	3.0	...

			for plt in relevant_numbers:
				data[ plt ]+= ( float(cols[relevant_numbers[plt][0]])+float(cols[relevant_numbers[plt][1]]) )
	src.close()
	for plt in data:
		out1.write(cross+'\t'+plt+'\t'+str(data[plt])+'\n')

