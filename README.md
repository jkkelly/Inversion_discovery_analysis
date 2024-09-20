# Inversion_discovery_analysis
python programs related to the identification and description of inversion polymorphisms in Mimulus guttatus 

Two programs to run Anchorwave aligning each chromosome of our nine alternative genomes to the IM767 reference genome: 00_anchorwave_preprocessing.py, 02_Anchor_wave_parallel.py

Four programs align genes to the inversions.  These locations are essential to the analysis of sequence variation and expression level variation.
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
   
   	python Pi.by.region.py [inv_ID] ---> [inv_ID].pairwise.pi.genomic_windows.txt, [inv_ID].AA_AD_DD.pi.genomic_windows.txt.gz (64 of each file)

Four programs analyze expression variation of genes



1.  Retrieve read counts from the original mapping of RNAseq reads to the genomes

	python reads.per.plant.py ---> raw.reads.per.plant.txt

2.  Perform a linear regression of standardized expression level of each F2 within each cross for each gene using the count of the alternative allele as the independent variable.

      python Test.all.genes.minCPM.py [datafile] ---> Gene_by_cross.tests.v2.txt

3.  Reading the "inversion relevance" of each gene, test whether genes within inversions exhibit greater cis genetic variation than in the remainder of the genome.  For all genes internal to inversions, partition the pairwise differences in expression effects within and between orientatins.

      python gene.tests_by_inversions.py  ---> Alpha.by.gene.txt

4.  The partitioning calculations from the previous program (gene.tests_by_inversions.py) is nested within a loop that permutes estimates across orientations within genes.
  
  
	python gene.tests_by_inversions.permute.py ---> perm.stats


The expression level input files were published with Veltsos and Kelly (2024) The quantitative genetics of gene expression in Mimulus guttatus. PLOS Genetics 20:e1011072. DOI: 10.1371/journal.pgen.1011072.  The programs listed above were applied using a draft assembly and annotation of the genome of the IM767 genome.  This "v1 assembly" is archived by JGI, the annotation of the v1 is included here as a gff3 file. 



