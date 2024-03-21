## Script to call SNPs and INDELs against a given MTBC reference genome
## Input $1 = a file name with a list of SRA accession numbers, one per line to call SNPs from

## Required tools and versions
#trimmomatic/0.39
#kraken2/2.1.2
#seqtk/1.3
#bwa/0.7.17
#picard/2.26.3
#bedtools/2.30.0
#samtools/1.16.1
#varscan/2.4.2
#gatk/4.2.4.0
#bcftools/1.11
#vcftools/0.1.16

## Chose reference sequence
REFERENCE="MTBC0_v1.1.fasta"

## Exclude file in BED format of regions of MTBC0 within which NOT to call SNPs/INDELS (e.g. PE/PPE)
## Derived from Goig et al. "ThePipeline" TB pipeline, lifted over to MTBC0 via halLIftover
EXCLUDE_FILE="MTBC0.REGIONStoDISCARD.bed"

## Directory containing FASTQ files
FASTQDIR="/pathtofastqfiles/"
FILTERED_FASTQ_DIR="/pathtofilteredfastqfiles/"
SCRATCH_DIR="/pathtoscratchfilespace/"

for i in `cat $1`
do
 echo -n "Processing "
 echo ${i}
 if test -f ${i}.gatk.raw.vcf
  then
  echo "Already processed, skipping"
 else 

 ### Determine if these are paired reads or not
 if test -f "${FILTERED_FASTQ_DIR}${i}.P1.filtered.fastq"
 then
  echo "already filtered"
 else 
  if test -f "${FASTQDIR}${i}_pass_1.fastq.gz"
  then
   ## The code below for read filtering and alignment mostly from Goig et al. SNP pipeline flowchart, with some modifications and adapted for kraken2.
   echo "Processing paired reads...."
   echo -n "trimmomatic... "
   java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 8 -phred33 -trimlog ${SCRATCH_DIR}${i}.trimlog ${FASTQDIR}${i}_pass_[12].fastq.gz ${SCRATCH_DIR}${i}.P1.clean.fastq ${SCRATCH_DIR}${i}.U1.clean.fastq ${SCRATCH_DIR}${i}.P2.clean.fastq ${SCRATCH_DIR}${i}.U2.clean.fastq TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:20
   echo "done."
   echo -n "kraken2... "
   kraken2 --threads 8 --db /home/luke/MTB_MRCA/krakenDB/KrakenDB_default_20221111 --use-names --report ${i}.kraken.report --paired ${SCRATCH_DIR}${i}.P[12].clean.fastq > ${SCRATCH_DIR}${i}.kraken.output
   echo "done"
   echo -n "Filtering out MTBC reads... "
   fgrep "Mycobacterium" ${SCRATCH_DIR}${i}.kraken.output | cut -f2 > ${SCRATCH_DIR}${i}.P1.filtered.readlist
   sed -e 's/\.1$/\.2/g' ${SCRATCH_DIR}${i}.P1.filtered.readlist > ${SCRATCH_DIR}${i}.P2.filtered.readlist
   seqtk subseq ${SCRATCH_DIR}${i}.P1.clean.fastq ${SCRATCH_DIR}${i}.P1.filtered.readlist | sed -e 's/\.1 / /g' > ${FILTERED_FASTQ_DIR}${i}.P1.filtered.fastq
   seqtk subseq ${SCRATCH_DIR}${i}.P2.clean.fastq ${SCRATCH_DIR}${i}.P2.filtered.readlist | sed -e 's/\.2 / /g' > ${FILTERED_FASTQ_DIR}${i}.P2.filtered.fastq
   echo "done."
  elif test -f "${FASTQDIR}${i}_pass.fastq.gz"
   then
    echo "Processing single end reads...."
    echo -n "trimmomatic... "
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 8 -phred33 -trimlog ${SCRATCH_DIR}${i}.trimlog ${FASTQDIR}${i}_pass.fastq.gz ${SCRATCH_DIR}${i}.P1.clean.fastq TRAILING:10 SLIDINGWINDOW:25:20 MINLEN:50
    echo "done."
    echo -n "kraken2... "
    kraken2 --threads 8 --db /home/luke/MTB_MRCA/krakenDB/KrakenDB_default_20221111 --use-names --report ${i}.kraken.report ${SCRATCH_DIR}${i}.P1.clean.fastq > ${SCRATCH_DIR}${i}.kraken.output
    echo "done"
    echo -n "Filtering out MTBC reads... "
    fgrep "Mycobacterium" ${SCRATCH_DIR}${i}.kraken.output | cut -f2 > ${SCRATCH_DIR}${i}.P1.filtered.readlist
    sed -e 's/\.1$/\.2/g' ${i}.P1.filtered.readlist > ${i}.P2.filtered.readlist
    seqtk subseq ${SCRATCH_DIR}${i}.P1.clean.fastq ${SCRATCH_DIR}${i}.P1.filtered.readlist | sed -e 's/\.1 / /g' > ${FILTERED_FASTQ_DIR}${i}.P1.filtered.fastq
    echo "done."
  else
   echo "${i} not found..."
   echo ${i} >> SRAs.not.found.txt
   continue
  fi
 fi
 echo "Mapping..."
 echo -n "bwa mem mapping... "
 bwa mem -t 8 ${REFERENCE} ${FILTERED_FASTQ_DIR}${i}.P[12].filtered.fastq | samtools view -bt ${REFERENCE} - | samtools sort -o ${i}.sort.bam
 echo "done."
 ## clean up temp files
 rm ${SCRATCH_DIR}${i}.P[12].clean.fastq ${SCRATCH_DIR}${i}.U[12].clean.fastq ${SCRATCH_DIR}${i}.kraken.output ${SCRATCH_DIR}${i}.P[12].filtered.readlist ${SCRATCH_DIR}${i}.trimlog

 java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${i}.sort.bam O=${i}.nd.sort.bam M=${i}.dup.metrix ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT 
 echo -n "done. genomeCoverage and stats..."
 genomeCoverageBed -ibam ${i}.nd.sort.bam -d -g ${REFERENCE} > ${i}.coverage
 samtools stats ${i}.nd.sort.bam > ${i}.samtools.stats
 echo "done."
 ### Variant calling
 echo -n "samtools index... "
 samtools index ${i}.nd.sort.bam 
 echo "done."
 java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${i}.nd.sort.bam O=${i}.nd.sort.renamed.bam SORT_ORDER=coordinate RGID=1 RGLB=NA RGPU=NA RGPL=illumina RGSM=${i} CREATE_INDEX=FALSE
 mv ${i}.nd.sort.renamed.bam ${i}.nd.sort.bam
 samtools index ${i}.nd.sort.bam
 gatk HaplotypeCaller --reference ${REFERENCE} --input ${i}.nd.sort.bam --native-pair-hmm-threads 8 --sample-ploidy 1 -ERC GVCF -O ${i}.gatk.raw.vcf
 ## Clean up remainder of temp files
 rm ${i}.sort.bam 
 fi 
done

## Combine GATK GVCFs into a multiple sample file gVCF file and call variants
for i in `cat $1`
do
 sed -i "s/Sample1/${i}/g" $i.gatk.raw.vcf
 GVCF_FILE="$GVCF_FILE --variant ${i}.gatk.raw.vcf"
done
gatk CombineGVCFs --call-genotypes true -R ${REFERENCE} ${GVCF_FILE} -O "Combined.GVCFs.gatk.raw.gvcf"
gatk GenotypeGVCFs --sample-ploidy 1 --call-genotypes true -R ${REFERENCE} -V "Combined.GVCFs.gatk.raw.gvcf" -O "Combined.called.GVCFs.gatk.raw.gvcf"

## Now filter called variants
## Filter SNPs.
gatk SelectVariants --reference ${REFERENCE} -V "Combined.called.GVCFs.gatk.raw.gvcf" --select-type SNP -O "Combined.called.GVCFs.gatk.raw.SNPs.gvcf"
gatk VariantFiltration --reference ${REFERENCE} -V "Combined.called.GVCFs.gatk.raw.SNPs.gvcf" --filter-name "LowQD" --filter-expression "QD<2.0" --filter-name "HighFS" --filter-expression "FS>60.0" --filter-name "MQRankSum" --filter-expression "MQRankSum<-12.5" --filter-name "Low40MQ" --filter-expression "MQ<40.0" --filter-name "LowReadPRS" --filter-expression "ReadPosRankSum<-8.0" --filter-name "LowDepth" --filter-expression "DP<10" -cluster 3 -window 50 -O "Combined.called.GVCFs.gatk.SNPs.gvcf"
 gatk SelectVariants --reference ${REFERENCE} -V "Combined.called.GVCFs.gatk.SNPs.gvcf" --exclude-non-variants --exclude-filtered -select-type SNP -O "Combined.called.GVCFs.gatk.fSNPs.gvcf"
 vcftools --exclude-bed ${EXCLUDE_FILE} --vcf "Combined.called.GVCFs.gatk.fSNPs.gvcf" --recode --stdout > "Combined.called.GVCFs.gatk.ffSNPs.gvcf"

## Filter INDELS.
 gatk SelectVariants --reference ${REFERENCE} -V "Combined.called.GVCFs.gatk.raw.gvcf" --select-type INDEL -O "Combined.called.GVCFs.gatk.raw.INDEL.gvcf"
 gatk VariantFiltration --reference ${REFERENCE} -V "Combined.called.GVCFs.gatk.raw.INDEL.gvcf" --filter-name "LowQD" --filter-expression "QD<2.0" --filter-name "HighFS" --filter-expression "FS>200.0" --filter-name "HighSOR" --filter-expression "SOR>10.0" --filter-name "LowReadPRS" --filter-expression "ReadPosRankSum<-20.0" -O "Combined.called.GVCFs.gatk.INDEL.gvcf"
 gatk SelectVariants --reference ${REFERENCE} -V "Combined.called.GVCFs.gatk.INDEL.gvcf" --exclude-non-variants --exclude-filtered -select-type INDEL -O "Combined.called.GVCFs.gatk.fINDEL.gvcf"
 vcftools --exclude-bed ${EXCLUDE_FILE} --vcf "Combined.called.GVCFs.gatk.fINDEL.gvcf" --recode --stdout > "Combined.called.GVCFs.gatk.ffINDEL.gvcf"


### Now create a SNP alignment of the called ffSNPs
bgzip -c "Combined.called.GVCFs.gatk.ffSNPs.gvcf" > "Combined.called.GVCFs.gatk.ffSNPs.gvcf.gz"
tabix -f -p vcf "Combined.called.GVCFs.gatk.ffSNPs.gvcf.gz"
for i in `cat $1`
do
 bcftools consensus -a X -f ${REFERENCE} -I  -p "${i} " -M ? -s $i -o ${i}.ffSNPs.fasta "Combined.called.GVCFs.gatk.ffSNPs.gvcf.gz"
 sed -i "s/\(>.*\) MTBC.*/\1/g" ${i}.SNPs.fasta
 sed -i 's/*/-/g' ${i}.ffSNPs.fasta
 cat ${i}.ffSNPs.fasta >> ffSNPs.alignment.fasta
done


