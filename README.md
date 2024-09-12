# Inversion_discovery_analysis
python programs related to the identification and description of inversion polymorphisms in Mimulus guttatus 

Four programs align genes to inversions
1. Determine the IM767 genes that were present in each alternative line assembly
	python Genes.present.by.cross.py ---> Genes_scored_by_cross.txt
2. Identify genes within or flanking each inversion
	python Genes.within.inversions.py --->  Genes.within_flanking.inversions
3. Relate inversion coordinates to the specific genes within the genetic maps
	python Inversions.to.genotypedgenes.py ---> inv.gg.txt
4. Determine the inversion genotype for all 64 inversion in each of the 1588 plants measured for gene expression in the eGWAS experiment
	python Inversion.genotypes.py ---> [inv_id].genotype.txt (64 files)


