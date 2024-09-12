# Inversion_discovery_analysis
python programs related to the identification and description of inversion polymorphisms in Mimulus guttatus 

Four programs align genes to inversions
1. Determine the IM767 genes that were present in each alternative line assembly
	python Genes.present.by.cross.py ---> Genes_scored_by_cross.txt
2. Identify genes within or flanking each inversion
	python Genes.within.inversions.py --->  Genes.within_flanking.inversions
3. Relate inversion coordinates to the specific genes within the genetic maps
	python Inversions.to.genotypedgenes.py ---> inv.gg.txt
4. Determine the inversion genotype for all 64 inversion in each of the 1588 plants measured for gene expression in Veltsos and Kelly (2024)
	python Inversion.genotypes.py ---> [inv_id].genotype.txt (64 files)

Four programs conduct the nucleotide sequence analyses
1.  Determine the number of aligned positions and differing positions of across the ten assemblies.  Summarizes are output for each gene on each chromosome.
	python Pi.by.gene.py ---> chrom+".gene.pairwise.pi.txt"
2.  Determine nucleotide diversity for each gene classifying lines according to whether standard or inverted within each inversion.
	python pi_within_between_karyotypes_genic.py --->   pi.byinv.genic.detailed.txt
3.  Apply polarization (Ancestral/Derived/Unclear) to each inversion and extract mean values (per inversion) for each statistic.		
	python pi_within_between_ancestral_derived_genic.py ---> pi.byinv.genic.AD.txt
4.  Calculate nucleotide diversity statistics in genomic windows from upstream to downstream of each inversion.

   	python Pi.by.region.py [inv_ID]
    
    		---> [inv_ID].pairwise.pi.genomic_windows.txt
    
    		---> [inv_ID].AA_AD_DD.pi.genomic_windows.txt.gz

The expression level input files were published with Veltsos and Kelly (2024) The quantitative genetics of gene expression in Mimulus guttatus. PLOS Genetics 20:e1011072. DOI: 10.1371/journal.pgen.1011072.  These analyses were conducted on a draft assembly and annotation of the genome.  The v1 assembly is archived by JGI, the annotation is included here as a gff3 file. 



