# MTBC<sub>0</sub>
Mycobacterium tuberculosis complex - Estimated Most recent common ancestor genomic sequence (MTBC<sub>0</sub>)

Repository to accompany the following article:
Harrison et al. An imputed ancestral reference genome for the Mycobacterium tuberculosis complex better captures structural genomic diversity for reference-based alignment workflows. To be sumitted to Microbial Genomics.

Contents:  
- MTBC0.fasta - the inputed ancestral sequence of the MTBC  
- MTBC0.cactus.hal - Hierachechal alignment format file containg the alignmnet of 30 closed genomes and estimated ancestral seqeunces
- MTBC0.cactus - Cactus control file
- NC_000962onMTBC0.gff3 - H37Rv (NCBI NC000962.3) genome annotation translated to MTBC<sub>0</sub> coordiantes using halLiftover  
- CP048071onMTBC0.gff3 - L8 (NCBI CP048071) genome annotation translated to MTBC<sub>0</sub> coordinates using halLiftover
- scripts/SNP_pipeline.sh - The GATK-based SNP pipeline used in this article
