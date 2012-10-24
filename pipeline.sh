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
ANALYSIS_DIR="04_Analysis"


TRIMMING_TMP=$TRIMMING_DIR/tmp
MAPPING_TMP=$MAPPING_DIR/tmp
FILTER_TMP=$FILTER_DIR/tmp
ANALYSIS_TMP=$ANALYSIS_DIR/tmp


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
	-n "${PARAMETERS_TABLE['bwa_sampe_n']}" \
	-N "${PARAMETERS_TABLE['bwa_sampe_N']}" \
	${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $MAPPING_DIR/$3_1.sai $MAPPING_DIR/$3_2.sai $TRIMMING_DIR/$3_1_paired.fq $TRIMMING_DIR/$3_2_paired.fq >$MAPPING_DIR/$3.sam 2>$MAPPING_TMP\_2
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

if [[ -e $MAPPING_DIR/$3.sam ]]; then
    echo "# .sam file exists ! Start to filter alignments ..." >> $LOG_DIR/$LOGFILE
else 
    echo "sam file does not exist"
    exit $?
fi

if [[ ! -e $FILTER_DIR ]]; then
    mkdir $FILTER_DIR 
fi

# Catch the sam file header
echo "$(date '+%Y%m%d %r') [Filter start] save sam header file" >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter start] save sam header file" >> $LOG_DIR/$LOGFILE

grep "^@" $MAPPING_DIR/$3.sam > $FILTER_DIR/header.txt
grep -v "^@" $MAPPING_DIR/$3.sam > $FILTER_DIR/$3_tmp.sam

# Count the initial number of reads in sam file
echo "$(date '+%Y%m%d %r') [Filter all reads] $3.sam $(wc -l $FILTER_DIR/$3_tmp.sam | wc -l) " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter all reads] $3.sam $(wc -l $FILTER_DIR/$3_tmp.sam | wc -l) " >> $LOG_DIR/$LOGFILE

# Remove non aligned reads
grep -v "*" $FILTER_DIR/$3_tmp.sam > $FILTER_DIR/$3_mapped.sam

# Count the number of aligned reads in sam file
echo "$(date '+%Y%m%d %r') [Filter aligned reads] $FILTER_DIR/$3_mapped.sam $(wc -l $FILTER_DIR/$3_mapped.sam | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter aligned reads] $FILTER_DIR/$3_mapped.sam $(wc -l $FILTER_DIR/$3_mapped.sam | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

# FILTRE sur MAPQ
awk '{if($5==0 || $5>=20) print}' $FILTER_DIR/$3_mapped.sam > $FILTER_DIR/$3_mapped_MAPQ.sam

# Count the number of aligned reads in sam file filtered for given MAPQ threshold
echo "$(date '+%Y%m%d %r') [Filter MAPQ] $FILTER_DIR/$3_mapped_MAPQ.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ.sam | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter MAPQ] $FILTER_DIR/$3_mapped_MAPQ.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ.sam | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

# Remove reads with more than x independent events
# Filter on XO and XM (6 cases): Xo+Xm <= 2 (default in ${PARAMETERS_TABLE["nb_of_independent_event"]})
remove_reads_with_more_than_x_independent_events 2 $FILTER_DIR/$3_mapped_MAPQ.sam >$FILTER_DIR/$3_mapped_MAPQ_XOXM.sam 2>$FILTER_DIR/$3_XOXM.excluded

# Count the number of reads in sam file filtered for given x independent events threshold
echo "$(date '+%Y%m%d %r') [Filter XOXM independent events] $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter XOXM independent events] $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

echo "$(date '+%Y%m%d %r') [Filter XOXM independent events] $FILTER_DIR/$3_XOXM.excluded $(wc -l $FILTER_DIR/$3_XOXM.excluded | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter XOXM independent events] $FILTER_DIR/$3_XOXM.excluded $(wc -l $FILTER_DIR/$3_XOXM.excluded | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

# Filter CIGAR code for I or D < Y bases threshold (default 5)
remove_reads_with_indels_size_greater_than_y_bases ${PARAMETERS_TABLE["microindel_size"]} $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam >$FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam 2>$FILTER_DIR/$3_INDELS.excluded

# Count the number of reads in sam file filtered for indels greater than given size threshold
echo "$(date '+%Y%m%d %r') [Filter indels] $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter indels] $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

echo "$(date '+%Y%m%d %r') [Filter indels] $FILTER_DIR/$3_INDELS.excluded $(wc -l $FILTER_DIR/$3_INDELS.excluded | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter indels] $FILTER_DIR/$3_INDELS.excluded $(wc -l $FILTER_DIR/$3_INDELS.excluded | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

# Add the header at the end of the filtering
cat $FILTER_DIR/header.txt $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam > tmp.sam
mv tmp.sam $FILTER_DIR/$3_mapped_MAPQ.sam

echo "$(date '+%Y%m%d %r') [Filter end] Add the header to filtered sam file" >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter end] Add the header to filtered sam file" >> $LOG_DIR/$LOGFILE

    # Convert SAM to BAM
if [[ -e $FILTER_DIR/$3_mapped_MAPQ.sam ]]; then
    echo "# sam file exists ! Start to convert sam to bam ..." >> $LOG_DIR/$LOGFILE
    samtools view -bS  $FILTER_DIR/$3_mapped_MAPQ.sam > $FILTER_DIR/$3.bam 2>$FILTER_TMP\_1
    echo "$(date '+%Y%m%d %r') [Filter: Convert sam to bam] $FILTER_DIR/$3.bam" >> $LOG_DIR/$LOGFILE
    cat $FILTER_TMP\_1 >> $FILTER_DIR/$FILTER_DIR_$DATE.log
else
    echo "SAM file does not exist"
    exit $?
fi

    # SORT BAM
if [[ -e $FILTER_DIR/$3.bam ]]; then
    echo "# bam file exists ! Start to sort bam file ..." >> $LOG_DIR/$LOGFILE
    samtools sort $FILTER_DIR/$3.bam $FILTER_DIR/$3_sorted 2>$FILTER_TMP\_2
    echo "$(date '+%Y%m%d %r') [Filter: Sorted bam] $FILTER_DIR/$3_sorted.bam" >> $LOG_DIR/$LOGFILE
    cat $FILTER_TMP\_2 >> $FILTER_DIR/$FILTER_DIR_$DATE.log
else
    echo "BAM file does not exist"
    exit $?
fi

    # INDEX BAM
if [[ -e $FILTER_DIR/$3_sorted.bam ]]; then
    echo "# Sorted bam file exists ! Start to index the file ..." >> $LOG_DIR/$LOGFILE
    samtools index $FILTER_DIR/$3_sorted.bam 2>$FILTER_TMP\_3
    # TO DO TEST $?
    echo "$(date '+%Y%m%d %r') [Filter: Indexed and sorted bam] $FILTER_DIR/$3_sorted.bam.bai" >> $LOG_DIR/$LOGFILE
    cat $FILTER_TMP\_3 >> $FILTER_DIR/$FILTER_DIR_$DATE.log
else
    echo "Sorted bam file does not exist"
    exit $?
fi

########################
# Section Analysis
########################

if [[ ! -e $ANALYSIS_DIR ]]; then
    mkdir $ANALYSIS_DIR 
fi

    # GENERATE BCF (PER SITE DEPTH, INFOS)

if [[ -e $FILTER_DIR/$3_sorted.bam ]]; then
    echo "# Sorted bam file exists ! Start to index fasta genome and compute ber base depth ..." >> $LOG_DIR/$LOGFILE    
    #samtools faidx ${PARAMETERS_TABLE["REFERENCE_GENOME_FASTA"]} 2>$ANALYSIS_TMP\_1
    # This step will be done one time
    samtools mpileup -B -Q 20 -u -f ${PARAMETERS_TABLE["REFERENCE_GENOME_FASTA"]} $FILTER_DIR/$3_sorted.bam > $ANALYSIS_DIR/$3.bcf 2>$ANALYSIS_TMP\_2
    # TODO mpileup parameters
    #echo "$(date '+%Y%m%d %r') [Analysis: faidx]">> $LOG_DIR/$LOGFILE
    #cat $ANALYSIS_TMP\_1 >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Analysis: mpileup] $3.bcf" >> $LOG_DIR/$LOGFILE
    cat $ANALYSIS_TMP\_2 >> $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log
else
    echo "Sorted bam file does not exist"
    exit $?
fi

    # VARIANT CALLING

if [[ -e $ANALYSIS_DIR/$3.bcf ]]; then
    echo "# bcf file exists ! Start to call variant" >> $LOG_DIR/$LOGFILE   
    bcftools view -cvg $ANALYSIS_DIR/$3.bcf > $ANALYSIS_DIR/$3.vcf 2>$ANALYSIS_TMP\_3
    echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] $3.vcf">> $LOG_DIR/$LOGFILE
    cat $ANALYSIS_TMP\_3 >> $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log
else
    echo "BCF file does not exist"
    exit $?
fi

    # VARIANT FILTERING
# TO DO (R)    
# vcfutils.pl varFilter $3.vcf > $3_filtre.vcf
    
    # COMPRESS and INDEX VCF file 

if [[ -e $ANALYSIS_DIR/$3.vcf ]]; then
    echo "# vcf file exists ! Start to compress it ..." >> $LOG_DIR/$LOGFILE   
    bgzip -fc $ANALYSIS_DIR/$3.vcf > $ANALYSIS_DIR/$3.vcf.gz
    # TO DO Standard error
    echo "$(date '+%Y%m%d %r') [Analysis: Compress vcf] $3.vcf.gz">> $LOG_DIR/$LOGFILE
    tabix -p vcf -f $ANALYSIS_DIR/$3.vcf.gz
    # TO DO Standard error
    echo "$(date '+%Y%m%d %r') [Analysis: tabix (Index compressed vcf)] $3.vcf.gz">> $LOG_DIR/$LOGFILE
else
    echo "VCF file does not exist"
    exit $?
fi  
 
    # SNP EFFECT PREDICTION

 echo "# Start snpEff ..." >> $LOG_DIR/$LOGFILE   

java -jar ${PARAMETERS_TABLE["SNPEFF_PATH"]}/snpEff.jar ${PARAMETERS_TABLE["snpeff_data"]} \
    -c ${PARAMETERS_TABLE["SNPEFF_PATH"]}/snpEff.config \
    -i ${PARAMETERS_TABLE["snpeff_inFile_format"]} \
    -o VCF $ANALYSIS_DIR/$3.vcf > $ANALYSIS_DIR/$3_snpeff.txt 2>$ANALYSIS_TMP\_4

echo "$(date '+%Y%m%d %r') [Analysis: snpeff]">> $LOG_DIR/$LOGFILE
mv snpEff_* $ANALYSIS_DIR/.
cat $ANALYSIS_TMP\_4 >> $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log


#####################
# CLEANING
####################




####################
# Run Pipeline.Rnw
###################