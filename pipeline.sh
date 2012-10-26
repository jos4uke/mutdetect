#! /bin/bash

# Equipe Dev 
# Script provisoire pour la détection de mutation 

# TODO: 
# for all commands, get the exit status code and if different of 0, exits the pipeline with error message; update functions to exit with error code
# for all file tests, add test to check if file is empty, then exits the pipeline with error message

########################
# SECTION CONFIGURATION
#######################

# Inclusion de la librairie de fonctions

PROD_PREFIX="/usr/local"
DEV_PREFIX="$(pwd)/.."
PREFIX=$DEV_PREFIX # TO BE CHANGED WHEN SWITCHING TO PROD
. $PREFIX/lib/pipeline_lib.inc

# Positionnement des variables

ARGS=3
DATE=$(date '+%Y_%m_%d_%H_%M_%S')
LOGFILE=$3_$DATE\_log.txt
WORKING_DIR=$(pwd)
PIPELINE_SHARED=$PREFIX/share/$(basename ${0%.*})
PIPELINE_DEFAULT_CONFIG=$PIPELINE_SHARED/etc/pipeline_default.config
PROD_PIPELINE_USER_CONFIG=$WORKING_DIR/pipeline_user.config
DEV_PIPELINE_USER_CONFIG=$PREFIX/pipeline_user.config
PIPELINE_USER_CONFIG=$DEV_PIPELINE_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD

LOG_DIR="log"
TRIMMING_DIR="01_Trimming"
MAPPING_DIR="02_Mapping"
FILTER_DIR="03_Filter"
ANALYSIS_DIR="04_Analysis"
REPORT_DIR="05_Report"
FORMATTED_LOG_DIR="$REPORT_DIR/formatted_log"

TRIMMING_TMP=$TRIMMING_DIR/tmp
MAPPING_TMP=$MAPPING_DIR/tmp
FILTER_TMP=$FILTER_DIR/tmp
ANALYSIS_TMP=$ANALYSIS_DIR/tmp

ERROR_TMP="/tmp/tmp_pipeline_error_${USER}_$DATE.log"

# DECLARE GLOBAL VARIABLE

declare -A PARAMETERS_TABLE


#==============================================
# TEST if enough args else print usage message
#==============================================
[[ $# -ne "$ARGS" ]] && { printf %s "\
Program: $(basename $0)
Version: none
Contact: IJPB Bioinformatics Dev Team

Usage: $(basename $0) SEQfile1 SEQfile2 ECHname

Arguments: SEQfile1 Forward read sequences file (Illumina fastq file)
           SEQfile2 Reverse read sequences file (Illumina fastq file)
           ECHname  Prefix to use for the analysis, i.e. prefix can be the sample name or anything else suitable       

Notes: 1. this pipeline version is actually able to perform variant calling

";
exit 1; }

#=======
# BEGIN
#=======

# Create log directory

if [[ -d $LOG_DIR ]]; then
    echo "$(date '+%Y%m%d %r') [$(basename $0)] Start running the pipeline." > $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [$(basename $0)] Executed command: $0 $*" >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." >> $LOG_DIR/$LOGFILE
else
    mkdir $LOG_DIR 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Log directory] Failed Log directory, $LOG_DIR, was not created." | tee -a $ERROR_TMP 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP."
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [$(basename $0)] Start running the pipeline." > $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [$(basename $0)] Executed command: $0 $*" >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." >> $LOG_DIR/$LOGFILE	
    fi
fi


# TEST if files exist 

if [[ -e $1 && -s $1 ]]; then
    if [[ -e $2 && -s $2 ]]; then
	echo "$(date '+%Y%m%d %r') [fastq input files] OK Input files exists and are not empty." >> $LOG_DIR/$LOGFILE
    else 
	echo "$(date '+%Y%m%d %r') [fastq input files] Failed Input file, $1, does not exist or is empty" >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." >> $LOG_DIR/$LOGFILE
	exit 3 
    fi 
else
    echo "$(date '+%Y%m%d %r') [fastq input files] Failed Input file, $2, does not exist or is empty" >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." >> $LOG_DIR/$LOGFILE
    exit 3 
fi

# TEST if the fastq files have the same number of reads
# else exit (error code 4) 

if (( $(wc -l $1 | awk '{print $1}') == $(wc -l $2 | awk '{print $1}') )); then
    echo "$(date '+%Y%m%d %r') [fastq input files] OK Input files have the same number of reads." >> $LOG_DIR/$LOGFILE
else
    echo "$(date '+%Y%m%d %r') [fastq input files] Failed Input files do not have the same number of reads!" >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 4." >> $LOG_DIR/$LOGFILE
    exit 4
fi

# TEST if pipeline_default.config exists and put the parameters into a hash table

if [[ -e $PIPELINE_DEFAULT_CONFIG ]]; then
    echo "$(date '+%Y%m%d %r') [get_pipeline_default_parameters] OK $PIPELINE_DEFAULT_CONFIG default config file exists! Let's check parameters validity ..." >> $LOG_DIR/$LOGFILE
    # load default config parameters
    get_pipeline_default_parameters $PIPELINE_DEFAULT_CONFIG 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [get_pipeline_default_parameters] Failed Default parameters were not loaded. An error occurs: $(cat $ERROR_TMP)" >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." >> $LOG_DIR/$LOGFILE
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [get_pipeline_default_parameters] OK Default config parameters were loaded successfully." >> $LOG_DIR/$LOGFILE
    fi
    # check default config params type validity
    check_params_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [check_params_validity] Failed Default parameters type checking generates errors." >> $LOG_DIR/$LOGFILE
	cat $ERROR_TMP >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." >> $LOG_DIR/$LOGFILE
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [check_params_validity] OK Default config parameters type checking was done successfully." >> $LOG_DIR/$LOGFILE
    fi
    # check default config params interval validity
    check_params_interval_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [check_params_interval_validity] Failed Default parameters interval checking generates errors." >> $LOG_DIR/$LOGFILE
	cat $ERROR_TMP >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." >> $LOG_DIR/$LOGFILE
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [check_params_interval_validity] OK Default config parameters interval checking was done successfully." >> $LOG_DIR/$LOGFILE
    fi
else 
    echo "$(date '+%Y%m%d %r') [get_pipeline_default_parameters] Failed $PIPELINE_DEFAULT_CONFIG file does not exist." >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." >> $LOG_DIR/$LOGFILE
    exit 3
fi 

# TEST if pipeline_user.config exists and then override default parameters if user defined parameters exist

if [[ -e $PIPELINE_USER_CONFIG ]]; then
    echo "$(date '+%Y%m%d %r') [get_pipeline_user_parameters] OK $PIPELINE_USER_CONFIG user config file exists! Let's check parameters validity ..." >> $LOG_DIR/$LOGFILE
    # load user config parameters
    get_pipeline_user_parameters $PIPELINE_USER_CONFIG 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [get_pipeline_user_parameters] Failed User parameters were not loaded. An error occurs: $(cat $ERROR_TMP)" >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." >> $LOG_DIR/$LOGFILE
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [get_pipeline_user_parameters] OK User config parameters were loaded successfully." >> $LOG_DIR/$LOGFILE
	if [[ -s $ERROR_TMP ]]; then
	    cat $ERROR_TMP >> $LOG_DIR/$LOGFILE
	fi
    fi
    # check user config params type validity
    check_params_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [check_params_validity] Failed User parameters type checking generates errors." >> $LOG_DIR/$LOGFILE
	cat $ERROR_TMP >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." >> $LOG_DIR/$LOGFILE
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [check_params_validity] OK User config parameters type checking was done successfully." >> $LOG_DIR/$LOGFILE
    fi
    # check user config params interval validity
    check_params_interval_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [check_params_interval_validity] Failed User parameters interval checking generates errors." >> $LOG_DIR/$LOGFILE
	cat $ERROR_TMP >> $LOG_DIR/$LOGFILE
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." >> $LOG_DIR/$LOGFILE
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [check_params_interval_validity] OK User config parameters interval checking was done successfully." >> $LOG_DIR/$LOGFILE
    fi
else 
    echo "$(date '+%Y%m%d %r') [get_pipeline_user_parameters] Failed $PIPELINE_USER_CONFIG file does not exist." >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." >> $LOG_DIR/$LOGFILE
    exit 3
fi

# Report parameters in LOGFILE

echo "$(date '+%Y%m%d %r') [parameters listing] OK All pipeline config parameters were loaded and checked successfully." >> $LOG_DIR/$LOGFILE
for i in "${!PARAMETERS_TABLE[@]}"
do
    echo -e "$i=${PARAMETERS_TABLE[$i]}" >> $LOG_DIR/$LOGFILE
done


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
	SLIDINGWINDOW:${PARAMETERS_TABLE["trimmo_sliding_window_size"]}:${PARAMETERS_TABLE["trimmo_sliding_window_qual"]} \
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
	SLIDINGWINDOW:${PARAMETERS_TABLE["trimmo_sliding_window_size"]}:${PARAMETERS_TABLE["trimmo_sliding_window_qual"]} \
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
	-t ${PARAMETERS_TABLE["bwa_aln_t"]} \
	${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $TRIMMING_DIR/$3_2_paired.fq > $MAPPING_DIR/$3_2.sai 2>$MAPPING_TMP\_1
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
	"${PARAMETERS_TABLE['BWA_REFERENCE_GENOME_INDEX']}" $MAPPING_DIR/$3_1.sai $MAPPING_DIR/$3_2.sai $TRIMMING_DIR/$3_1_paired.fq $TRIMMING_DIR/$3_2_paired.fq >$MAPPING_DIR/$3.sam 2>$MAPPING_TMP\_2
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
echo "$(date '+%Y%m%d %r') [Filter all reads] $3.sam $(wc -l $FILTER_DIR/$3_tmp.sam | awk '{print $1}') " >> $FILTER_DIR/$FILTER_DIR_$DATE.log
echo "$(date '+%Y%m%d %r') [Filter all reads] $3.sam $(wc -l $FILTER_DIR/$3_tmp.sam | awk '{print $1}') " >> $LOG_DIR/$LOGFILE

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
    -o vcf $ANALYSIS_DIR/$3.vcf > $ANALYSIS_DIR/$3_snpeff.vcf 2>$ANALYSIS_TMP\_4
    
java -jar ${PARAMETERS_TABLE["SNPEFF_PATH"]}/snpEff.jar ${PARAMETERS_TABLE["snpeff_data"]} \
    -c ${PARAMETERS_TABLE["SNPEFF_PATH"]}/snpEff.config \
    -i ${PARAMETERS_TABLE["snpeff_inFile_format"]} \
    -o txt $ANALYSIS_DIR/$3.vcf > $ANALYSIS_DIR/$3_snpeff.txt 2>$ANALYSIS_TMP\_4

echo "$(date '+%Y%m%d %r') [Analysis: snpeff]">> $LOG_DIR/$LOGFILE
mv snpEff_* $ANALYSIS_DIR/.
cat $ANALYSIS_TMP\_4 >> $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log


#####################
# CLEANING
####################




####################
# Run Pipeline.Rnw
###################

# TODO: Pipeline.Rnw is a resource which should not be in the user working directory
# This resource file, like the other pipeline resources, should be placed in some location shared by all users: like /usr/local/share/<pipeline_name_dir>
# Consider to move all template and other resource files: .Rnw, .docx, etc. to this shared location
# Do not forget to create the corresponding variable path: PIPELINE_SHARED=$PREFIX/share/<pipeline_name_dir> with PREFIX=/usr/local
#
echo "$(date '+%Y%m%d %r') [Rnw filling]">> $LOG_DIR/$LOGFILE

if [[ ! -e $REPORT_DIR ]]; then
    mkdir $REPORT_DIR 
fi

if [[ ! -e  $FORMATTED_LOG_DIR ]]; then
    mkdir $FORMATTED_LOG_DIR
fi

rnw_trimming_details_subsection $LOG_DIR/$LOGFILE $FORMATTED_LOG_DIR

rnw_alignment_filtering_subsection $LOG_DIR/$LOGFILE $FORMATTED_LOG_DIR

#=====
# END
#=====
echo "$(date '+%Y%m%d %r') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
echo -n "$(date '+%Y%m%d %r') [$(basename $0)] Elapsed time: " | tee -a $LOG_DIR/$LOGFILE 2>&1
echo |awk -v time="$SECONDS" '{print strftime("%Hh:%Mm:%Ss", time, 1)}' | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y%m%d %r') [$(basename $0)] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y%m%d %r') [$(basename $0)] More information about the analysis can be found in $LOG_DIR/$LOGFILE" | tee -a $LOG_DIR/$LOGFILE 2>&1

exit 0