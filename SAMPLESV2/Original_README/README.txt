README

SAMPLES from 1000 Genomes: HG00119, HG00188, HG00258, HG00513, HG01383


1. downloaded sliced BAM/BAI files for each sample for tp53 gene with coordinates 17:7565097-7590868
	data slicer found at: http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer

2. downloaded sliced VCF file with all 5 samples with coordinates 17:7565097-7590868
	data slicer found at: http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer

3. using MuSiC_REMIX, bam_to_bed, converted each sample’s BAM file to BED file
	parameters: delta = 10 outFile = SAMPLE_NAME.bed

4. using MuSiC_REMIX, beds_to_covered_regions, converted the BED files of all samples to single adequately covered regions BED file
	parameters: delta = 0.75 outFile = adequately_covered_regions.bed

5. quality-checked all BED files using IGV. Visually made sure that edges were aligned.

6. created Region-Of-Interest (ROI) file from tp53 exon locations
	name = roi.bed
	exon locations found at: http://www.rcsb.org/pdb/gene/TP53

7. using MuSiC_REMIX, intersect_roi_and_popCovBed, created merged bed file from roi.bed and adequately_covered_regions.bed. As it so happened, all ROI regions were within the adequately_covered_regions BED.
	name = roi_popCov_intersect.bed

8. remix MuSiC_REMIXs calculate_mutation_rate.py
	mutation rate = ( SIGMA( #samples_with_mutation_within_ROI ) ) / ( ROI_basepair_size * #samples_total )