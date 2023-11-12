# MTBC<sub>0</sub>
Mycobacterium tuberculosis complex - Estimated Most recent common ancestor genomic sequence (MTBC<sub>0</sub>)

Repository to accompany the following article:
Harrison et al. An imputed ancestral reference genome for the Mycobacterium tuberculosis complex better captures structural genomic diversity for reference-based alignment workflows. Submitted to Microbial Genomics.

Contents:  
- MTBC<sub>0</sub>_v1.1.fasta - the inputed ancestral sequence of the MTBC  
- MTBC<sub>0</sub>.v1.1.cactus.hal - Hierachechal alignment format file containg the alignmnet of 30 closed genomes and estimated ancestral seqeunces
- MTBC<sub>0</sub>.cactus - Cactus control file
- MTBC<sub>0</sub>v1.1_PGAP_annot.gff - PGAP sequence annotation for the MTBC<sub>0</sub> sequence
- MTBC0_Goigetal_regions_toDiscard.bed - H37Rv-based annotation of regions to exclude from SNP calls (e.g. PE/PPE, IS, etc), from Goig et al., 2020 (Comas et al. TB Pipeline), translated to MTBC<sub>0</sub> coordinates using halLiftover  
- liftover/ - collection of liftovers of genome annotations (M. canettii, L8 and H37Rv) onto MTBC<sub>0</sub>, including a 1:1 position liftover for H37Rv in .bed and spreadsheet formats
- other.reference.squences/ - directory containing H37Rv (NC_000962.3, downloaded from NCBI) and the ancestral sequence of Comas et al. (2010)
- scripts/SNP_pipeline.sh - The GATK-based SNP pipeline used in this article
- scripts/Figure1.script - The R script and H37RvonMTBC0 alignment in MAF format to generated windowed sequence identity tracks for Figure 1
- scripts/bed2gff.R - helper script to translate a bed file (Created with gff2bed) back to GFF
- trees/ - raw phylogenetic trees used to generate Supplementary Figures 1, 2 and 3
