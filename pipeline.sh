#! /bin/bash

# Equipe Dev 
# Script provisoire pour la détection de mutation 

# Inclusion de la librairie de fonctions

PREFIX="/projects/ARABIDOPSIS/SCRIPTS/PIPELINE"
. $PREFIX/lib/pipeline_lib.inc

# Positionnement des variables

ARGS=3
DATE=$(date '+%Y_%m_%d_%H_%M')
LOGFILE=$3_$DATE\_log.txt
WORKING_DIR=$(pwd)
PIPELINE_DEFAULT_CONFIG="/projects/ARABIDOPSIS/SCRIPTS/PIPELINE/pipeline_default.config"

LOG_DIR="log"
TRIMMING_DIR="01_Trimming"
MAPPING_DIR="02_Mapping"
FILTER_DIR="03_Filter"

TRIMMING_TMP=$TRIMMING_DIR/tmp
MAPPING_TMP=$MAPPING_DIR/tmp

# DECLARE GLOBAL VARIABLE
declare -A PARAMETERS_TABLE

# TEST if enough args

if [ $# -ne "$ARGS" ]
then
  echo "Usage: `$0` SEQfile1 SEQfile2 ECHname"
  exit $?
fi

# TEST if files exist 

if [[ -e $1 || -e $2 ]]; then
	echo "# $(date '+%Y%m%d %r') Input Files exists ! Let's check sequences quality ..." >> $LOG_DIR/$LOGFILE
	#tmp=${1%\.*};
	#FASTQ1_NAME=${tmp##*/}
	#tmp=${2%\.*};
	#FASTQ2_NAME=${tmp##*/}
else 
	echo "File does not exist"
	exit $?
fi 

# TODO TEST if the fastq files have the same number of reads
# if it is not the case exit 


# TEST if pipeline_default.config exists and put the parameters into a hash table

if [[ -e $PIPELINE_DEFAULT_CONFIG ]]; then
	echo "# $(date '+%Y%m%d %r') [get_pipeline_default_parameters] default config file exists ! Let's check parameters validity ..." >> $LOG_DIR/$LOGFILE
	get_pipeline_default_parameters $PIPELINE_DEFAULT_CONFIG
else 
	echo "File does not exist"
	exit $?
fi 

# Report parameters in LOGFILE


########################
# SECTION TRIMMING
#######################

# If they do not already exist, create directories to store QC and Filtering results

if [[ ! -e $TRIMMING_DIR ]]; then
    mkdir $TRIMMING_DIR 
fi

if [[ ! -e  $TRIMMING_DIR/$3_1_Qual_Raw_Reads ]]; then
    mkdir $TRIMMING_DIR/$3_1_Qual_Raw_Reads
fi

if [[ ! -e  $TRIMMING_DIR/$3_2_Qual_Raw_Reads ]]; then
    mkdir $TRIMMING_DIR/$3_2_Qual_Raw_Reads 
fi

# Check raw reads quality

echo "# $(date '+%Y%m%d %r') [fastqc] QC directories exist ! Let's check raw reads quality ..." >> $LOG_DIR/$LOGFILE

fastqc $1 -o $TRIMMING_DIR/$3_1_Qual_Raw_Reads 2>$TRIMMING_DIR/$TRIMMING_DIR_$DATE.log >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log
fastqc $2 -o $TRIMMING_DIR/$3_2_Qual_Raw_Reads 2>$TRIMMING_DIR/$TRIMMING_DIR_$DATE.log >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log

# TODO: Ici recuperer les infos et graphes de qualité

# Trim reads by using Trimmomatic

echo "# $(date '+%Y%m%d %r') [Trimmomatic]  Let's trim raw reads ..." >> $LOG_DIR/$LOGFILE

if [[ ${PARAMETERS_TABLE["QUAL_ENCODING"]} -eq 33 ]]; then
    java -classpath /usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads ${PARAMETERS_TABLE["trimmo_thread_number"]} \
	-phred33 \
	-trimlog $TRIMMING_DIR/LogTrim \
	$1 \
	$2 \
	$TRIMMING_DIR/$3_1_paired.fq $TRIMMING_DIR/$3_1_single.fq \
	$TRIMMING_DIR/$3_2_paired.fq $TRIMMING_DIR/$3_2_single.fq  \
	LEADING:${PARAMETERS_TABLE["trimmo_leading_qual_min"]}  \
	TRAILING:${PARAMETERS_TABLE["trimmo_trailing_qual_min"]}  \
	SLIDINGWINDOW:${PARAMETERS_TABLE["trimmo_slinding_window_size"]}:${PARAMETERS_TABLE["trimmo_sliding_window_qual"]} \
	MINLEN:${PARAMETERS_TABLE["trimmo_min_length"]} 2>$TRIMMING_TMP \
	>> $TRIMMING_TMP
    else 
    java -classpath /usr/local/src/Trimmomatic-0.22/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticPE \
	-threads ${PARAMETERS_TABLE["trimmo_thread_number"]} \
	-phred64 \
	-trimlog $TRIMMING_DIR/LogTrim \
	$1 \
	$2 \
	$TRIMMING_DIR/$3_1_paired.fq $TRIMMING_DIR/$3_1_single.fq \
	$TRIMMING_DIR/$3_2_paired.fq $TRIMMING_DIR/$3_2_single.fq  \
	LEADING:${PARAMETERS_TABLE["trimmo_leading_qual_min"]}  \
	TRAILING:${PARAMETERS_TABLE["trimmo_trailing_qual_min"]}  \
	SLIDINGWINDOW:${PARAMETERS_TABLE["trimmo_slinding_window_size"]}:${PARAMETERS_TABLE["trimmo_sliding_window_qual"]} \
	MINLEN:${PARAMETERS_TABLE["trimmo_min_length"]} 2>$TRIMMING_TMP \
	>> $TRIMMING_TMP
fi

cat $TRIMMING_TMP >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log
cat $TRIMMING_TMP >> $LOG_DIR/$LOGFILE

# Check trimmed reads Quality

if [[ ! -e  $TRIMMING_DIR/$3_1_Qual_Trim_Reads ]]; then
    mkdir $TRIMMING_DIR/$3_1_Qual_Trim_Reads
fi

if [[ ! -e  $TRIMMING_DIR/$3_2_Qual_Trim_Reads ]]; then
    mkdir $TRIMMING_DIR/$3_2_Qual_Trim_Reads 
fi

fastqc $TRIMMING_DIR/$3_1_paired.fq -o $TRIMMING_DIR/$3_1_Qual_Trim_Reads 2>$TRIMMING_DIR/$TRIMMING_DIR_$DATE.log >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log
fastqc $TRIMMING_DIR/$3_2_paired.fq -o $TRIMMING_DIR/$3_2_Qual_Trim_Reads 2>$TRIMMING_DIR/$TRIMMING_DIR_$DATE.log >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log


########################
# SECTION MAPPING
#######################

# Using BWA to map reads
# UP to two mismatches
# To allow, by exmple, one polymorphism + the mutation in the same read

# Create Mapping directory if does not exist

if [[ ! -e $MAPPING_DIR ]]; then
    mkdir $MAPPING_DIR 
fi

if [[ -e $TRIMMING_DIR/$3_1_paired.fq && -e $TRIMMING_DIR/$3_2_paired.fq ]]; then
    echo "# $(date '+%Y%m%d %r') [bwa aln] Input Files exist ! Run bwa aln ..." >> $LOG_DIR/$LOGFILE
    bwa aln \
	-n ${PARAMETERS_TABLE["bwa_aln_n"]} \
	-R ${PARAMETERS_TABLE["bwa_aln_R"]} \
	-t ${PARAMETERS_TABLE["bwa_aln_t"]} \
	${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $TRIMMING_DIR/$3_1_paired.fq > $MAPPING_DIR/$3_1.sai 2>$MAPPING_DIR/$MAPPING_DIR_$DATE.log &
    bwa aln \
	-n ${PARAMETERS_TABLE["bwa_aln_n"]} \
	-R ${PARAMETERS_TABLE["bwa_aln_R"]} \
	-t ${PARAMETERS_TABLE["bwa_aln_t"]} ${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $TRIMMING_DIR/$3_2_paired.fq > $MAPPING_DIR/$3_2.sai 2>$MAPPING_TMP\_1
else 
    echo "File does not exists"
    exit $?
fi 


# TEST if sai files exists

if [[ -e $MAPPING_DIR/$3_1.sai && -e $MAPPING_DIR/$3_2.sai ]]; then
    echo "# $(date '+%Y%m%d %r') [bwa sampe] .sai Files exist ! Run bwa sampe ..." >> $LOG_DIR/$LOGFILE
    bwa sampe \
	-n ${PARAMETERS_TABLE["bwa_sampe_n"]} \
	-N ${PARAMETERS_TABLE["bwa_sampe_N"]} \
	${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $MAPPING_DIR/$3_1.sai $MAPPING_DIR/$3_2.sai $TRIMMING_DIR/$3_1_paired.fq $TRIMMING_DIR/$3_2_paired.fq > $MAPPING_DIR/$3.sam 2>$MAPPING_TMP\_2
else
    echo ".sai files do not exists"
    exit $?
fi 

cat $MAPPING_TMP\_1 >> $MAPPING_DIR/$MAPPING_DIR_$DATE.log
cat $MAPPING_TMP\_2 >> $MAPPING_DIR/$MAPPING_DIR_$DATE.log

grep "^\[infer_isize\] inferred external" $MAPPING_DIR/$MAPPING_DIR_$DATE.log >> $LOG_DIR/$LOGFILE

########################
# SECTION FILTERING
#######################

if [[ -e $3.sam ]]; then
    echo "# .sam file exists ! Start to filter alignments ..." >> $LOG_DIR/$LOGFILE
else 
    echo "sam file does not exist"
    exit $?
fi

if [[ ! -e $FILTER_DIR ]]; then
    mkdir $FILTER_DIR 
fi

# Catch the sam file header
grep "^@" $MAPPING_DIR/$3.sam > $FILTER_DIR/header.txt
grep -v "^@" $MAPPING_DIR/$3.sam > $FILTER_DIR/$3_tmp.sam


# Count the initial number of reads in sam file
echo "$(date '+%Y%m%d %r') $3.sam $(wc -l $3.sam | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

head -8 $3.sam > header.txt
sed '1,8d' $3.sam > $3_tmp.sam

# Remove non aligned reads
grep -v "*" $3_tmp.sam > $3_mapped.sam

# Count the number of aligned reads in sam file
echo "$(date '+%Y%m%d %r') $3_mapped.sam $(wc -l $3_mapped.sam | awk '{print $1}') " >> $LOGFILE

# Remove reads with more than 2 mismatches
#grep "NM:i:[012][[:space:]]" $3_mapped.sam > $3_mapped_NM2.sam

awk '{if($5<0 || $5>20) print}' $3_mapped.sam >  $3_mapped_MAPQ_tmp.sam
cat header.txt $3_mapped_MAPQ_tmp.sam >  $3_mapped_MAPQ.sam

# Count the number of aligned reads in sam file
echo "$(date '+%Y%m%d %r') $3_mapped_MAPQ.sam $(wc -l $3_mapped_MAPQ.sam | awk '{print $1-8}') " >> $LOGFILE

# Filter CIGAR code for I or D < 5 (make it a variable)
#TODO


# Filter on XO and XM (3 cases)
#TODO: XM=0:XO<=2 ; XM=1:XO<=1; XM=2:XO=0 

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
}' $3_mapped_MAPQ.sam > $3_mapped_MAPQ_BH.sam

echo "$(date '+%Y%m%d %r') $3_mapped_MAPQ_BH.sam $(wc -l $3_mapped_MAPQ_BH.sam | awk '{print $1}') " >> $LOGFILE

# Remove hits with 0<MAPQ<20 

#awk '{if($5<0 || $5>20) print}' $3_mapped_NM2_BH.sam >  $3_mapped_NM2_BH_MAPQ_tmp.sam
#cat header.txt $3_mapped_NM2_BH_MAPQ_tmp.sam >  $3_mapped_NM2_BH_MAPQ.sam

#echo "$(date '+%Y%m%d %r') $3_mapped_NM2_BH_MAPQ.sam $(wc -l $3_mapped_NM2_BH_MAPQ.sam | awk '{print $1-8}') " >> $LOGFILE

    # Convert SAM to BAM

#samtools view -bS  $3_mapped_NM2_BH_MAPQ.sam > $3.bam

    # SORT BAM

#samtools sort $3.bam $3_sorted

    # INDEX BAM

#samtools index $3_sorted.bam

    # GENERATE BCF (PER SITE INFOS)
    
#samtools faidx $GENOME_FASTA
    
#samtools mpileup -B -Q 20 -u -f $GENOME_FASTA $3_sorted.bam > $3.bcf
    
    # VARIANT CALLING

#bcftools view -cvg $3.bcf > $3.vcf

    # VARIANT FILTERING
    
#    vcfutils.pl varFilter $3.vcf > $3_filtre.vcf
    
    # COMPRESS VCF 

#    bgzip $3_filtre.vcf

#    tabix -p vcf $3_filtre.vcf.gz
    
    # SNP EFFECT PREDICTION

#    java -jar snpEff.jar athaliana130 $3_filtre.vcf > $3_filtre_snpeff.txt
#fi