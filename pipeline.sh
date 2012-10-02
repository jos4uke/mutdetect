#! /bin/bash

# Equipe Dev 
# Script provisoire pour la détection de mutation 


ARGS=3
GENOME_INDEX="/projects/ARABIDOPSIS/INDEX/Col0_bwa_0.6.1/TAIR10_All_Chrom.fas"
LOGFILE=$3_$(date '+%Y_%m_%d_%H_%M')_log.txt
WORKING_DIR=$(pwd)
GENOME_FASTA=

# TEST if enough args

if [ $# -ne "$ARGS" ]
then
  echo "Usage: `$0` SEQfile1 SEQfile2 ECHname"
  exit $?
fi

if [[ -e $1 || -e $2 ]]; then
	echo "# $(date '+%Y%m%d %r') Input Files exists ! Let's check sequences quality ..." >> $LOGFILE
else 
	echo "File does not exists"
	exit $?
fi 

mkdir  $3_1_Qual_Raw_Reads $3_2_Qual_Raw_Reads
fastqc $1 -o $3_1_Qual_Raw_Reads
fastqc $2 -o $3_2_Qual_Raw_Reads

# Trimmomatic

java -classpath /usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE -threads 4 -phred33 -trimlog LogTrim $1 $2 $1_paired.fq $1_single.fq $2_paired.fq $2_single.fq  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:50

mkdir $3_1_Qual_Trim_Reads $3_2_Qual_Trim_Reads
fastqc $1_paired.fq -o $3_1_Qual_Trim_Reads
fastqc $2_paired.fq -o $3_2_Qual_Trim_Reads

if [[ -e $1_paired.fq && -e $2_paired.fq ]]; then
	echo "# $(date '+%Y%m%d %r') Input Files exists ! Run bwa aln ..." >> $LOGFILE
else 
	echo "File does not exists"
	exit $?
fi 
   

# Using BWA to map reads
# UP to two mismatches
# To alliw, by exmple, one polymorphism + the mutation in the same read

bwa aln -n 2 -R 4 -t 4 $GENOME_INDEX $1_paired.fq > $3_1.sai

bwa aln -n 2 -R 4 -t 4 $GENOME_INDEX $2_paired.fq > $3_2.sai

# TEST if sai files exists

if [[ -e $3_1.sai && -e $3_2.sai ]]; then
    echo "# .sai Files exists ! Run bwa sampe ..." >> $LOGFILE
else
    echo ".sai files do not exists"
	exit $?
fi 

echo "# bwa sampe -n 4 -N 0 -a 1000 $INDEX $3_1.sai $3_2.sai $1 $2 > $3_coller0_No_Filtre.sam" >> $LOGFILE

    
bwa sampe -n 4 -a 1000 -N 0 $GENOME_INDEX $3_1.sai $3_2.sai $1_paired.fq $2_paired.fq > $3.sam

# TEST if sam file exists

if [[ -e $3.sam ]]; then
    echo "# .sam file exists ! Start to filter alignments ..." >> $LOGFILE
else 
    echo "sam file does not exist"
    exit $?
fi


# Count the number of reads in sam file
echo "$(date '+%Y%m%d %r') $3.sam $(wc -l $3.sam | awk '{print $1}') " >> $LOGFILE

head -8 $3.sam > header.txt
sed '1,8d' $3.sam > $3_tmp.sam

# Remove non aligned reads
grep -v "*" $3_tmp.sam > $3_mapped.sam

# Count the number of aligned reads in sam file
echo "$(date '+%Y%m%d %r') $3_mapped.sam $(wc -l $3_mapped.sam | awk '{print $1}') " >> $LOGFILE

# Remove reads with more than 2 mismatches
grep "NM:i:[012]" $3_mapped.sam > $3_mapped_NM2.sam

# Count the number of aligned reads in sam file
echo "$(date '+%Y%m%d %r') $3_mapped_NM2.sam $(wc -l $3_mapped_NM2.sam | awk '{print $1}') " >> $LOGFILE

# Remove suboptimal hits => Do not delete alignment but Modify reported alignment

awk '{
if($0 ~ "X1:i:[123]"){
whereX0=match($0,"X0:i:")
X0=substr($0,whereX0+5,1);

whereX1=match($0,"X1:i:")
X1=substr($0,whereX1+5,1);

whereNM=match($0,"NM:i:")
NM=substr($0,whereNM+5,1);

if(X0 == 1){
sub("\t"$NF,"",$0)
}
else{
split($NF,XAfield,":");
split(XAfield[3],XArray,";");
for(i=1;i<=(X0+X1-1);i++){
split(XArray[i],SHarray,",")
if(SHarray[4] > NM){
sub(XArray[i]";","",$0);
}}}}
sub("X1:i:[123]","X1:i:0",$0)
print
}' $3_mapped_NM2.sam > $3_mapped_NM2_BH.sam

echo "$(date '+%Y%m%d %r') $3_mapped_NM2_BH.sam $(wc -l $3_mapped_BH.sam | awk '{print $1}') " >> $LOGFILE

# Remove hits with 0<MAPQ<20 

awk '{if($5<0 || $5>20) print}' $3_mapped_NM2_BH.sam >  $3_mapped_NM2_BH_MAPQ_tmp.sam
cat header.txt $3_mapped_NM2_BH_MAPQ_tmp.sam >  $3_mapped_NM2_BH_MAPQ.sam

echo "$(date '+%Y%m%d %r') $3_mapped_NM2_BH_MAPQ.sam $(wc -l $3_mapped_NM2_BH_MAPQ.sam | awk '{print $1-8}') " >> $LOGFILE

    # Convert SAM to BAM

samtools view -bS  $3_mapped_NM2_BH_MAPQ.sam > $3.bam

    # SORT BAM

samtools sort $3.bam $3_sorted

    # INDEX BAM

samtools index $3_sorted.bam

    # GENERATE BCF (PER SITE INFOS)
    
samtools faidx $GENOME_FASTA
    
samtools mpileup -B -Q 20 -u -f $GENOME_FASTA $3_sorted.bam > $3.bcf
    
    # VARIANT CALLING

bcftools view -cvg $3.bcf > $3.vcf

    # VARIANT FILTERING
    
#    vcfutils.pl varFilter $3.vcf > $3_filtre.vcf
    
    # COMPRESS VCF 

#    bgzip $3_filtre.vcf

#    tabix -p vcf $3_filtre.vcf.gz
    
    # SNP EFFECT PREDICTION

#    java -jar snpEff.jar athaliana130 $3_filtre.vcf > $3_filtre_snpeff.txt
#fi