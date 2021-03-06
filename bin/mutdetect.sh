#! /bin/bash

# Equipe Dev 
# Script pour la détection de mutation 

VERSION=0.0.6-e3ef10e

########################
# SECTION CONFIGURATION
#######################

# Inclusion des librairies de fonctions

PROD_PREFIX="/usr/local"
DEV_PREFIX="$(pwd)/.."
PREFIX=$PROD_PREFIX # TO BE CHANGED WHEN SWITCHING TO PROD
. $PREFIX/share/mutdetect/lib/mutdetect_lib.inc
. $PROD_PREFIX/share/bash-common/lib/bash-common_lib.inc

# Positionnement des variables

ARGS=4
DATE=$(date '+%Y_%m_%d_%H_%M_%S')
LOGFILE=$3_$DATE\_log.txt
WORKING_DIR=$(pwd)
PIPELINE_SHARED=$PREFIX/share/$(basename ${0%.*})
PIPELINE_DEFAULT_CONFIG=$PIPELINE_SHARED/etc/mutdetect_default.config
PROD_PIPELINE_USER_CONFIG=$WORKING_DIR/mutdetect_user.config
DEV_PIPELINE_USER_CONFIG=$PIPELINE_SHARED/etc/mutdetect_user.config
PIPELINE_USER_CONFIG=$PROD_PIPELINE_USER_CONFIG # TO BE CHANGED WHEN SWITCHING TO PROD

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

ERROR_TMP="/tmp/tmp_mutdetect_error_${USER}_$DATE.log"

# DECLARE GLOBAL VARIABLE

declare -A PARAMETERS_TABLE
GENOME_ALIASES_LIST=()
GENOME_ALIAS=$4

###########################
# SECTION GENOMES ALIASES
###########################

# GET GENOMES/INDEXES PATHS
GENOMES_BASE_PATH=$(grep -e "^GENOMES_BASE_PATH=" $PIPELINE_DEFAULT_CONFIG | awk -F"=" '{print $2}')
INDEXES_BASE_PATH=$(grep -e "^INDEXES_BASE_PATH=" $PIPELINE_DEFAULT_CONFIG | awk -F"=" '{print $2}')
eval "BWA_INDEXES=$(grep -e "^BWA_INDEXES=" $PIPELINE_DEFAULT_CONFIG | awk -F"=" '{print $2}')"
eval "SAMTOOLS_INDEXES=$(grep -e "^SAMTOOLS_INDEXES=" $PIPELINE_DEFAULT_CONFIG | awk -F"=" '{print $2}')"
SNPEFF_PATH=$(grep -e "^SNPEFF_PATH=" $PIPELINE_DEFAULT_CONFIG | awk -F"=" '{print $2}')
eval "SNPEFF_DATA=$(grep -e "^SNPEFF_DATA=" $PIPELINE_DEFAULT_CONFIG | awk -F"=" '{print $2}')"

#echo ${GENOMES_BASE_PATH}
#echo ${INDEXES_BASE_PATH}
#echo ${BWA_INDEXES}
#echo ${SAMTOOLS_INDEXES}
#echo ${SNPEFF_DATA}

# GET GENOME ALIASES LIST
genomes_list_fasta=($(get_genomes_dir_list_with_one_fasta ${GENOMES_BASE_PATH} 2>/dev/null))
genomes_w_bwa_idx=($(get_genomes_list_with_bwa_index ${BWA_INDEXES}/$(get_tool_version "bwa") "${genomes_list_fasta[@]}" 2>/dev/null))
genomes_w_samtools_idx=($(get_genomes_list_with_samtools_index ${SAMTOOLS_INDEXES}/$(get_tool_version "samtools") "${genomes_w_bwa_idx[@]}" 2>/dev/null))
GENOME_ALIASES_LIST=($(get_genomes_list_with_snpeff_annot ${SNPEFF_DATA} "${genomes_w_samtools_idx[@]}" 2>/dev/null))

#==============================================
# TEST if enough args else print usage message
#==============================================
[[ $# -ne "$ARGS" ]] && { printf %s "\
Program: $(basename $0)
Version: $VERSION
Contact: IJPB Bioinformatics Dev Team

Usage: $(basename $0) SEQfile1 SEQfile2 ECHname Genome_alias

Arguments: SEQfile1 	Forward read sequences file (Illumina fastq file)
           SEQfile2 	Reverse read sequences file (Illumina fastq file)
           ECHname  	Prefix to use for the analysis, i.e. prefix can be the sample name or anything else suitable
           Genome_alias	Genome alias, see the available aliases in the following list

Genome aliases list:
$(for ga in "${GENOME_ALIASES_LIST[@]}"; do
	echo "- $ga"
done)

Notes: 1. this pipeline version is actually able to perform variant calling

";
exit 1; }

#=======
# BEGIN
#=======

# Create log directory

if [[ -d $LOG_DIR ]]; then
    echo "$(date '+%Y%m%d %r') [$(basename $0)] Start running the pipeline (version: $VERSION)." | tee $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Log directory] OK $LOG_DIR directory already exists. Will write log files in this directory." >> $LOG_DIR/$LOGFILE
else
    mkdir $LOG_DIR 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Log directory] Failed Log directory, $LOG_DIR, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [$(basename $0)] Start running the pipeline." | tee $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [$(basename $0)] Executed command: $0 $*" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Log directory] OK $LOG_DIR directory was created sucessfully. Will write log files in this directory." >> $LOG_DIR/$LOGFILE	
    fi
fi

# Check for genome alias
echo "$(date '+%Y%m%d %r') [$(basename $0)] Checking genome alias ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
if [[ -n $GENOME_ALIAS ]]; then
	if [[ $(in_array $GENOME_ALIAS "${GENOME_ALIASES_LIST[@]}") -eq 0 ]]; then
		echo "$(date '+%Y%m%d %r') [$(basename $0)] OK Genome alias: $GENOME_ALIAS is valid." | tee -a $LOG_DIR/$LOGFILE 2>&1
	else
    	echo "$(date '+%Y%m%d %r') [Check genome alias] Genome alias, $GENOME_ALIAS, is not valid." | tee -a $LOG_DIR/$LOGFILE 2>&1
    	echo "$(date '+%Y%m%d %r') [Check genome alias] Please refer to the pipeline usage message to get the full list of genome aliases." | tee -a $LOG_DIR/$LOGFILE 2>&1
    	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
    	exit
	fi
else
    echo "$(date '+%Y%m%d %r') [Check genome alias] Genome alias is null, this argument is mandatory." | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Check genome alias] Please refer to the pipeline usage message to get the full list of genome aliases." | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit	
fi

# TEST if fastq input files exist 

if [[ -s $1 ]]; then
    if [[ -s $2 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: fastq input files] OK Input files exists and are not empty." | tee -a $LOG_DIR/$LOGFILE 2>&1
    else 
	echo "$(date '+%Y%m%d %r') [Check config: fastq input files] Failed Input file, $2, does not exist or is empty" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 3 
    fi 
else
    echo "$(date '+%Y%m%d %r') [Check config: fastq input files] Failed Input file, $1, does not exist or is empty" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3 
fi

# TEST if the fastq files have the same number of reads
# else exit (error code 4) 

if (( $(wc -l $1 | awk '{print $1}') == $(wc -l $2 | awk '{print $1}') )); then
    echo "$(date '+%Y%m%d %r') [Check config: fastq input files] OK Input files have the same number of reads." | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y%m%d %r') [Check config: fastq input files] Failed Input files do not have the same number of reads!" | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 4." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 4
fi

# TEST if mutdetect_default.config exists and put the parameters into a hash table

if [[ -e $PIPELINE_DEFAULT_CONFIG ]]; then
    echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_default_parameters] OK $PIPELINE_DEFAULT_CONFIG default config file exists! Let's check parameters validity ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
    # load default config parameters
    get_mutdetect_default_parameters $PIPELINE_DEFAULT_CONFIG 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_default_parameters] Failed Default parameters were not loaded. An error occurs: $(cat $ERROR_TMP)" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_default_parameters] OK Default config parameters were loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
    # check default config params type validity
    check_params_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: check_params_validity] Failed Default parameters type checking generates errors." | tee -a $LOG_DIR/$LOGFILE 2>&1
	cat $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [Check config: check_params_validity] OK Default config parameters type checking was done successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
    # check default config params interval validity
    check_params_interval_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: check_params_interval_validity] Failed Default parameters interval checking generates errors." | tee -a $LOG_DIR/$LOGFILE 2>&1
	cat $ERROR_TMP | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [Check config: check_params_interval_validity] OK Default config parameters interval checking was done successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
else 
    echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_default_parameters] Failed $PIPELINE_DEFAULT_CONFIG file does not exist." | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi 

# TEST if mutdetect_user.config exists and then override default parameters if user defined parameters exist

if [[ -s $PIPELINE_USER_CONFIG ]]; then
    echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_user_parameters] OK $PIPELINE_USER_CONFIG user config file exists! Let's check parameters validity ..." | tee -a $LOG_DIR/$LOGFILE 2>&1
    # load user config parameters
    get_mutdetect_user_parameters $PIPELINE_USER_CONFIG 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_user_parameters] Failed User parameters were not loaded. An error occurs: $(cat $ERROR_TMP)" | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_user_parameters] OK User config parameters were loaded successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
	if [[ -s $ERROR_TMP ]]; then
	    cat $ERROR_TMP 2>&1 >> $LOG_DIR/$LOGFILE
	fi
    fi
    # check user config params type validity
    check_params_validity >> $LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: check_params_validity] Failed User parameters type checking generates errors." | tee -a $LOG_DIR/$LOGFILE 2>&1
	cat $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [Check config: check_params_validity] OK User config parameters type checking was done successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
    # check user config params interval validity
    check_params_interval_validity >>$LOG_DIR/$LOGFILE 2>$ERROR_TMP
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Check config: check_params_interval_validity] Failed User parameters interval checking generates errors." | tee -a $LOG_DIR/$LOGFILE 2>&1
	cat $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	echo "$(date '+%Y%m%d %r') [Check config: check_params_interval_validity] OK User config parameters interval checking was done successfully." | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
else 
    echo "$(date '+%Y%m%d %r') [Check config: get_mutdetect_user_parameters] No $PIPELINE_USER_CONFIG file does exist, will proceed with default config parameters." | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# eval parameters value for variable expansion
echo "$(date '+%Y%m%d %r') [Parameters eval] Will eval all pipeline config parameters for variable expansion" | tee -a $LOG_DIR/$LOGFILE 2>&1
for i in "${!PARAMETERS_TABLE[@]}"
do
	if [[ $(echo "${PARAMETERS_TABLE[$i]}" | grep -e "^\\$" | wc -l) -eq 1 ]]
	then
		k1=${PARAMETERS_TABLE[$i]%%/*}
		k2=${k1:1}
		#echo $k2
		eval "$k2=${PARAMETERS_TABLE[$k2]}"
		eval "k3=${PARAMETERS_TABLE[$i]}"
		#echo $k2		
		#echo $k3
		PARAMETERS_TABLE[$i]=$(echo "$k3")
		unset k1 k2 k3
		unset $(echo $k2)
		#echo $k2
		#echo $k3
	fi
done
echo "$(date '+%Y%m%d %r') [Parameters eval] Eval succesfully all pipeline config parameters for variable expansion." | tee -a $LOG_DIR/$LOGFILE 2>&1

# set genome seq/indexes/annot variables 
PARAMETERS_TABLE["REFERENCE_GENOME_FASTA"]="${PARAMETERS_TABLE["GENOMES_BASE_PATH"]}"/$GENOME_ALIAS/$(ls "${PARAMETERS_TABLE["GENOMES_BASE_PATH"]}"/$GENOME_ALIAS | grep -e "$GENOME_ALIAS\.m*fas*$")
PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]="${PARAMETERS_TABLE["BWA_INDEXES"]}"/$(get_tool_version "bwa")/$GENOME_ALIAS/$GENOME_ALIAS
PARAMETERS_TABLE["SAMTOOLS_REFERENCE_GENOME_INDEX"]="${PARAMETERS_TABLE["SAMTOOLS_INDEXES"]}"/$(get_tool_version "samtools")/$GENOME_ALIAS/$GENOME_ALIAS

#echo ${PARAMETERS_TABLE["REFERENCE_GENOME_FASTA"]}
#echo ${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]}
#echo ${PARAMETERS_TABLE["SAMTOOLS_REFERENCE_GENOME_INDEX"]}

# Report parameters in LOGFILE

echo -e "$(date '+%Y%m%d %r') [Parameters listing] OK All pipeline config parameters were loaded and checked successfully. \n\nList of parameters:\n" | tee -a $LOG_DIR/$LOGFILE 2>&1
for i in "${!PARAMETERS_TABLE[@]}"
do
    echo -e "$i=${PARAMETERS_TABLE[$i]}" | tee -a $LOG_DIR/$LOGFILE 2>&1
done

########################
# SECTION TRIMMING
#######################

# If they do not already exist, create directories to store QC and Trimming results
if [[ -d $TRIMMING_DIR ]]; then
    echo  "$(date '+%Y%m%d %r') [Trimming and quality control] Starting" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Trimming directory] OK $TRIMMING_DIR directory already exists. Will write trimming output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
else
    mkdir $TRIMMING_DIR 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Trimming directory] Failed Trimming directory, $TRIMMING_DIR, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Trimming and quality control] Start quality control and trimming process" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Trimming directory] OK $TRIMMING_DIR directory was created sucessfully. Will write trimming output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE	
    fi
fi

if [[ -d  $TRIMMING_DIR/$3_1_Qual_Raw_Reads ]]; then
    echo "$(date '+%Y%m%d %r') [Trimming: forward Raw Read Quality Control directory] OK $TRIMMING_DIR/$3_1_Qual_Raw_Reads directory already exists. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
else
    mkdir $TRIMMING_DIR/$3_1_Qual_Raw_Reads 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Trimming: forward Raw Read Quality Control directory] Failed Quality directory, $TRIMMING_DIR/$3_1_Qual_Raw_Reads, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1	
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Trimming: forward Raw Read Quality Control directory] OK $TRIMMING_DIR/$3_1_Qual_Raw_Reads directory was created sucessfully. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE	
    fi
fi

if [[ -d  $TRIMMING_DIR/$3_2_Qual_Raw_Reads ]]; then
    echo "$(date '+%Y%m%d %r') [Trimming: reverse Raw Read Quality Control directory] OK $TRIMMING_DIR/$3_2_Qual_Raw_Reads directory already exists. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
else
    mkdir $TRIMMING_DIR/$3_2_Qual_Raw_Reads 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Trimming: Reverse Raw Read Quality Control directory] Failed Quality directory, $TRIMMING_DIR/$3_2_Qual_Raw_Reads, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Trimming: reverse Raw Read Quality Control directory] OK $TRIMMING_DIR/$3_2_Qual_Raw_Reads directory was created sucessfully. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE	
    fi
fi

# Check raw reads quality

echo "$(date '+%Y%m%d %r') [Quality control: fastqc] QC directories exist ! Let's check raw reads quality ..." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE

echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Check raw reads quality in file $1" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
fastqc $1 -o $TRIMMING_DIR/$3_1_Qual_Raw_Reads 2>&1 | tee $ERROR_TMP 2>&1 >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Quality contol: fastqc] Failed An error occured during quality control of raw reads in file $1" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] OK Quality Control of raw reads in file $1 was done sucessfully." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
fi

echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Check raw reads quality in file $2" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
fastqc $2 -o $TRIMMING_DIR/$3_2_Qual_Raw_Reads 2>&1 | tee $ERROR_TMP 2>&1 >> $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Failed An error occured during quality control of raw reads in file $2" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] OK Quality Control of raw reads in file $2 was done sucessfully." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
fi

# Check for fastqc quality control failures report
if [[ $(toupper ${PARAMETERS_TABLE["bypass_fastqc_failure_report_checking"]}) == "FALSE" ]]; then
    echo "$(date '+%Y%m%d %r') [check_fastqc_quality_failure_report] Check for fastqc quality control failures in raw reads file $1" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    check_fastqc_quality_failure_report $TRIMMING_DIR/$3_1_Qual_Raw_Reads/$(basename $1)_fastqc/summary.txt 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Failed An error occured during fastqc quality failures quality checking for raw reads file $1" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	if [[ -s $ERROR_TMP ]]; then
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] A fastqc quality failure was reported for raw reads file $1: $(cat $ERROR_TMP)" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] If you do not want this pipeline to perform fasqtc quality failures checking, please change the following option in your configuration file to bypass_fastqc_failure_report_checking=TRUE" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn
	else
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] OK No fastqc quality failures reported for raw reads file $1." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
    fi

    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Check for fastqc quality control failures in raw reads file $2" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    check_fastqc_quality_failure_report $TRIMMING_DIR/$3_2_Qual_Raw_Reads/$(basename $2)_fastqc/summary.txt 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log | teee -a $LOG_DIR/$LOGFILE 2>&1
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Failed An error occured during fastqc quality failures quality checking for raw reads file $2" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	if [[ -s $ERROR_TMP ]]; then
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] A fastqc quality failure was reported for raw reads file $2: $(cat $ERROR_TMP)" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] If you do not want this pipeline to perform fasqtc quality failures checking, please change the following option in your configuration file to bypass_fastqc_failure_report_checking=TRUE" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn
	else
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] OK No fastqc quality failures reported for raw reads file $2." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
    fi
else
    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Skip No fastqc quality failures report checking for raw reads file $1." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Skip No fastqc quality failures report checking for raw reads file $2." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Trim reads by using Trimmomatic

echo "$(date '+%Y%m%d %r') [Trimming: Trimmomatic] Let's trim raw reads ..." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

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
	MINLEN:${PARAMETERS_TABLE["trimmo_min_length"]} \
	2>$ERROR_TMP >> $TRIMMING_TMP
        rtrn=$?
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
	MINLEN:${PARAMETERS_TABLE["trimmo_min_length"]} \
	2>$ERROR_TMP  >> $TRIMMING_TMP
        rtrn=$?
fi

if [[ $rtrn -ne 0 ]]; then
    cat $TRIMMING_TMP 2>&1 | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Trimming: Trimmomatic] Failed An error occured during trimming process." | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    cat $TRIMMING_TMP | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Trimming: Trimmomatic] OK Trimming of raw reads was done sucessfully." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi 

# Check trimmed reads Quality

echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Let's check trimmed reads quality ..." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

if [[ -d  $TRIMMING_DIR/$3_1_Qual_Trim_Reads ]]; then
    echo "$(date '+%Y%m%d %r') [Quality control: Forward Trimmed Read Quality Control directory] OK $TRIMMING_DIR/$3_1_Qual_Trim_Reads directory already exists. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $TRIMMING_DIR/$3_1_Qual_Trim_Reads 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Quality control: Forward Trimmed Read Quality Control directory] Failed Quality directory, $TRIMMING_DIR/$3_1_Qual_Trim_Reads, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1	
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Quality control: Forward Trimmed Read Quality Control directory] OK $TRIMMING_DIR/$3_1_Qual_Trim_Reads directory was created sucessfully. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
fi

if [[ -d  $TRIMMING_DIR/$3_2_Qual_Trim_Reads ]]; then
    echo "$(date '+%Y%m%d %r') [Quality control: Reverse Trimmed Read Quality Control directory] OK $TRIMMING_DIR/$3_2_Qual_Trim_Reads directory already exists. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $TRIMMING_DIR/$3_2_Qual_Trim_Reads 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Quality control: Reverse Trimmed Read Quality Control directory] Failed Quality directory, $TRIMMING_DIR/$3_2_Qual_Trim_Reads, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Quality control: Reverse Trimmed Read Quality Control directory] OK $TRIMMING_DIR/$3_2_Qual_Trim_Reads directory was created sucessfully. Will write fastqc output in this directory." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
fi

echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Check trimmed reads quality in file $TRIMMING_DIR/$3_1_paired.fq" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
fastqc $TRIMMING_DIR/$3_1_paired.fq -o $TRIMMING_DIR/$3_1_Qual_Trim_Reads 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Failed An error occured during quality control of trimmed reads in file $TRIMMING_DIR/$3_1_paired.fq" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] OK Quality Control of trimmed reads in file $TRIMMING_DIR/$3_1_paired.fq was done sucessfully." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Check trimmed reads quality in file $TRIMMING_DIR/$3_2_paired.fq" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fastqc $TRIMMING_DIR/$3_2_paired.fq -o $TRIMMING_DIR/$3_2_Qual_Trim_Reads 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] Failed An error occured during quality control of trimmed reads in file $TRIMMING_DIR/$3_2_paired.fq" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Quality control: fastqc] OK Quality Control of trimmed reads in file $TRIMMING_DIR/$3_2_paired.fq was done sucessfully." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Check for fastqc quality control failures report
if [[ $(toupper ${PARAMETERS_TABLE["bypass_fastqc_failure_report_checking"]}) == "FALSE" ]]; then
    echo "$(date '+%Y%m%d %r') [check_fastqc_quality_failure_report] Check for fastqc quality control failures in trimmed reads file $3_1_paired.fq" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    check_fastqc_quality_failure_report $TRIMMING_DIR/$3_1_Qual_Trim_Reads/$3_1_paired.fq_fastqc/summary.txt 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Failed An error occured during fastqc quality failures quality checking for trimmed reads file $3_1_paired.fq" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	if [[ -s $ERROR_TMP ]]; then
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] A fastqc quality failure was reported for trimmed reads file $3_1_paired.fq: $(cat $ERROR_TMP)" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] If you do not want this pipeline to perform fasqtc quality failures checking, please change the following option in your configuration file to bypass_fastqc_failure_report_checking=TRUE" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn
	else
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] OK No fastqc quality failures reported for trimmed reads file $3_1_paired.fq." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
    fi

    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Check for fastqc quality control failures in trimmed reads file $3_2_paired.fq" | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    check_fastqc_quality_failure_report $TRIMMING_DIR/$3_2_Qual_Trim_Reads/$3_2_paired.fq_fastqc/summary.txt 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Failed An error occured during fastqc quality failures quality checking for trimmed reads file $3_2_paired.fq" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	if [[ -s $ERROR_TMP ]]; then
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] A fastqc quality failure was reported for trimmed reads file $3_2_paired.fq: $(cat $ERROR_TMP)" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] If you do not want this pipeline to perform fasqtc quality failures checking, please change the following option in your configuration file to bypass_fastqc_failure_report_checking=TRUE" | tee -a $ERROR_TMP 2>&1 | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn
	else
	    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] OK No fastqc quality failures reported for trimmed reads file $3_2_paired.fq." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
    fi
else
    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Skip No fastqc quality failures report checking for trimmed reads file $3_1_paired.fq." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Quality control: check_fastqc_quality_failure_report] Skip No fastqc quality failures report checking for trimmed reads file $3_2_paired.fq." | tee -a $TRIMMING_DIR/$TRIMMING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi


########################
# SECTION MAPPING
#######################

# Using BWA to map reads
# default: 2 independent events i.e. considering mismatches and micro-indels
# To allow, for example, at most one polymorphism + one indel of length 1 in the same read (the best match)
# For the read pair mate, not the best match, this restriction is not applied and this read is "saved"
# Nevertheless, to filter those last ones, see filtering section: apply max 2 independent events and micro-indel max length of 5 by default
# default: multiple hits allowed up to 30, beyond this limit the coverage will decrease as a function of N_i (the number of reads at the position i in some repeated region) / M (the number of repetitions of the region), so take care of big gene families and highly repeated regions 

# Create Mapping directory if does not exist

if [[ -d $MAPPING_DIR ]]; then
    echo "$(date '+%Y%m%d %r') [Mapping] Starting mapping process" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Mapping: directory] OK $MAPPING_DIR directory already exists. Will write mapping output in this directory." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $MAPPING_DIR
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Mapping: directory] Failed Mapping directory, $MAPPING_DIR, was not created." | tee $ERROR_TMP 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Mapping] Starting mapping process" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Mapping: directory] OK $MAPPING_DIR directory was created sucessfully. Will write mapping output in this directory." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
fi

# Perform reads alignment
## Search suffix array for reads coordinates

echo "$(date '+%Y%m%d %r') [Mapping] Starting suffix array search process" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

if [[ -s $TRIMMING_DIR/$3_1_paired.fq ]]; then
    if [[ -s $TRIMMING_DIR/$3_2_paired.fq ]]; then
	echo "$(date '+%Y%m%d %r') [Mapping: bwa aln] Input Files exist ! Run bwa aln ..." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	bwa aln \
	    -n ${PARAMETERS_TABLE["bwa_aln_n"]} \
	    -R ${PARAMETERS_TABLE["bwa_aln_R"]} \
	    -t ${PARAMETERS_TABLE["bwa_aln_t"]} \
	    ${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $TRIMMING_DIR/$3_1_paired.fq > $MAPPING_DIR/$3_1.sai 2>$ERROR_TMP &
	bwa_aln_pid=$!
	
	bwa aln \
	    -n ${PARAMETERS_TABLE["bwa_aln_n"]} \
	    -R ${PARAMETERS_TABLE["bwa_aln_R"]} \
	    -t ${PARAMETERS_TABLE["bwa_aln_t"]} \
	    ${PARAMETERS_TABLE["BWA_REFERENCE_GENOME_INDEX"]} $TRIMMING_DIR/$3_2_paired.fq > $MAPPING_DIR/$3_2.sai 2>$MAPPING_TMP\_1
	rtrn2=$?

        # TODO: this error handling section should wait for the end of the first instance of bwa aln => cf "kill -0 bwa_aln_pid" to tell the current status of the process; and wait bwa_aln_pid will return the exit status 
	wait $bwa_aln_pid
	rtrn1=$?
	if [[ $rtrn1 -ne 0 ]]; then
	    echo "$(date '+%Y%m%d %r') [Mapping: bwa aln] Failed An error occured during suffix array search using bwa aln with file $TRIMMING_DIR/$3_1_paired.fq" | tee -a $ERROR_TMP 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn1." | tee -a $ERROR_TMP 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn1
	else
	    cat $ERROR_TMP >> $MAPPING_DIR/$MAPPING_DIR_$DATE.log
	    echo "$(date '+%Y%m%d %r') [Mapping: bwa aln] OK Suffix array search using bwa aln with file $TRIMMING_DIR/$3_1_paired.fq was done sucessfully." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
        #

	if [[ $rtrn2 -ne 0 ]]; then
	    echo "$(date '+%Y%m%d %r') [Mapping: bwa aln] Failed An error occured during suffix array search using bwa aln with file $TRIMMING_DIR/$3_2_paired.fq" | tee -a $MAPPING_TMP\_1 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn2." | tee -a $MAPPING_TMP\_1 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $MAPPING_TMP\_1." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn2
	else
	    cat $MAPPING_TMP\_1 >> $MAPPING_DIR/$MAPPING_DIR_$DATE.log
	    echo "$(date '+%Y%m%d %r') [Mapping: bwa aln] OK Suffix array search using bwa aln with file $TRIMMING_DIR/$3_2_paired.fq was done sucessfully." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi
    else 
	echo "$(date '+%Y%m%d %r') [Mapping: fastq input files] Failed Input file, $3_2_paired.fq, does not exist or is empty" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 3
    fi
else
    echo "$(date '+%Y%m%d %r') [Mapping: fastq input files] Failed Input file, $3_1_paired.fq, does not exist or is empty" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3 
fi 

## Convert SA coordinates to chromosomal coordinates for paired-end reads

echo "$(date '+%Y%m%d %r') [Mapping] Starting paired-end reads alignment process" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

if [[ -s $MAPPING_DIR/$3_1.sai ]]; then
    if [[ -s $MAPPING_DIR/$3_2.sai ]]; then
	echo "$(date '+%Y%m%d %r') [Mapping: bwa sampe] .sai Files exist ! Run bwa sampe ..." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	bwa sampe \
	    -n "${PARAMETERS_TABLE['bwa_sampe_n']}" \
	    -N "${PARAMETERS_TABLE['bwa_sampe_N']}" \
	    "${PARAMETERS_TABLE['BWA_REFERENCE_GENOME_INDEX']}" $MAPPING_DIR/$3_1.sai $MAPPING_DIR/$3_2.sai $TRIMMING_DIR/$3_1_paired.fq $TRIMMING_DIR/$3_2_paired.fq >$MAPPING_DIR/$3.sam 2>$MAPPING_TMP\_2
	rtrn=$?
	if [[ $rtrn -ne 0 ]]; then
	    echo "$(date '+%Y%m%d %r') [Mapping: bwa sampe] Failed An error occured during paired-end reads alignment using bwa sampe" | tee -a $MAPPING_TMP\_2 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $MAPPING_TMP\_2 2>&1 | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $MAPPING_TMP\_2." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	    exit $rtrn
	else
	    cat $MAPPING_TMP\_2 >> $MAPPING_DIR/$MAPPING_DIR_$DATE.log
	    grep "^\[infer_isize\] inferred external" $MAPPING_DIR/$MAPPING_DIR_$DATE.log >> $LOG_DIR/$LOGFILE
	    echo "$(date '+%Y%m%d %r') [Mapping :bwa sampe] OK Paired-end reads alignment using bwa sampe was done sucessfully." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 >> $LOG_DIR/$LOGFILE
	fi
    else
	echo "$(date '+%Y%m%d %r') [Mapping: bwa sampe] Failed Input sai file, $MAPPING_DIR/$3_2.sai, does not exist or is empty" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 3 
    fi
else
    echo "$(date '+%Y%m%d %r') [Mapping: bwa sampe] Failed Input sai file, $MAPPING_DIR/$3_1.sai, does not exist or is empty" | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." | tee -a $MAPPING_DIR/$MAPPING_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3    
fi

########################
# SECTION FILTERING
#######################

if [[ -d $FILTER_DIR ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering] Starting filtering process" 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: directory] OK $FILTER_DIR directory already exists. Will write filtering output in this directory." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $FILTER_DIR
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Filtering: directory] Failed Filtering directory, $FILTER_DIR, was not created." 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." 2>&1 | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Filtering] Starting filtering process" 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: directory] OK $FILTER_DIR directory was created sucessfully. Will write filtering output in this directory." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
fi

if [[ -s $MAPPING_DIR/$3.sam ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: sam input files] OK Sam input files exists and are not empty." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y%m%d %r') [Filtering: sam input files] Failed Input file, $MAPPING_DIR/$3.sam, does not exist or is empty" 2>&1 | tee $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2&>1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

# Catch the sam file header
echo "$(date '+%Y%m%d %r') [Filtering: sam header] Backuping sam header to file" 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
grep "^@" $MAPPING_DIR/$3.sam > $FILTER_DIR/header.txt 2>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: sam header] Failed Sam header was not backuped to file. An error occured." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: sam header] Error: $(cat $ERROR_TMP)" 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." 2>&1 | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: sam header] OK Sam header saved to file $FILTER_DIR/header.txt." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Count the initial number of reads in sam file
echo "$(date '+%Y%m%d %r') [Filtering: all reads] $3.sam $(get_total_reads_count $FILTER_DIR/$3.sam 2>/dev/null)" 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

# Remove non aligned reads
echo "$(date '+%Y%m%d %r') [Filtering: sam unmapped reads] Removing unmapped reads." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
get_mapped_reads $MAPPING_DIR/$3.sam 2>$ERROR_TMP >$FILTER_DIR/$3_mapped.sam
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: sam unmapped reads] Failed Unmapped reads were not filtered and no mapped reads were saved to sam file. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: sam unmapped reads] Error: $(cat $ERROR_TMP)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: sam unmapped reads] OK Only mapped reads were saved to sam file $FILTER_DIR/$3_mapped.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: sam unmapped reads] Info: $(cat $ERROR_TMP)." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Count the number of aligned reads in sam file
echo "$(date '+%Y%m%d %r') [Filtering: aligned reads] $FILTER_DIR/$3_mapped.sam $(wc -l $FILTER_DIR/$3_mapped.sam | awk '{print $1}') " 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

# FILTRE sur MAPQ
echo "$(date '+%Y%m%d %r') [Filtering: MAPQ] Filtering mapped reads on MAPQ." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
awk -v mapq_th=${PARAMETERS_TABLE["MAPQ_min"]} '{if($5==0 || $5>=mapq_th) print}' $FILTER_DIR/$3_mapped.sam > $FILTER_DIR/$3_mapped_MAPQ.sam 2>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: MAPQ] Failed Mapped reads were not filtered on MAPQ and no mapped and MAPQ filtered reads were saved to sam file. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: MAPQ] Error: $(cat $ERROR_TMP)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: MAPQ] OK Only mapped and MAPQ filtered (==0 && >=${PARAMETERS_TABLE['MAPQ_min']}) reads were saved to sam file $FILTER_DIR/$3_mapped.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Count the number of aligned reads in sam file filtered for given MAPQ threshold
echo "$(date '+%Y%m%d %r') [Filtering: MAPQ] $FILTER_DIR/$3_mapped_MAPQ.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ.sam | awk '{print $1}') " 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

# Remove reads with more than x independent events
# Filter on XO and XM (6 cases): Xo+Xm <= 2 (default in ${PARAMETERS_TABLE["nb_of_independent_event"]})
echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] Filtering mapped reads on number of independent events (XO+XM tags)." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
remove_reads_with_more_than_x_independent_events ${PARAMETERS_TABLE["nb_of_independent_event"]} $FILTER_DIR/$3_mapped_MAPQ.sam >$FILTER_DIR/$3_mapped_MAPQ_XOXM.sam 2>$FILTER_DIR/$3_XOXM.excluded
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] Failed Mapped reads were not filtered on independent events and no mapped and XOXM filtered reads were saved to sam file. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] Error: $(cat $ERROR_TMP)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] OK Only mapped and MAPQ filtered reads with XOXM (<=${PARAMETERS_TABLE['nb_of_independent_event']}) were saved to sam file $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] OK mapped and MAPQ filtered reads but with XOXM (>${PARAMETERS_TABLE['nb_of_independent_event']}) were saved to excluded file $FILTER_DIR/$3_XOXM.excluded." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Count the number of reads in sam file filtered for given x independent events threshold
echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam | awk '{print $1}') " 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
echo "$(date '+%Y%m%d %r') [Filtering: XOXM independent events] $FILTER_DIR/$3_XOXM.excluded $(wc -l $FILTER_DIR/$3_XOXM.excluded | awk '{print $1}') " 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

# Filter CIGAR code for I or D < Y bases threshold (default 5)
echo "$(date '+%Y%m%d %r') [Filtering: indels] Filtering mapped reads on indels size (CIGAR code)." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
remove_reads_with_indels_size_greater_than_y_bases ${PARAMETERS_TABLE["microindel_size"]} $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam >$FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam 2>$FILTER_DIR/$3_INDELS.excluded
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: indels] Failed Mapped reads were not filtered on indels size and no mapped and indels size filtered reads were saved to sam file. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: indels] Error: $(cat $ERROR_TMP)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: indels] OK Only mapped and MAPQ, XOXM filtered reads with indels size (<=${PARAMETERS_TABLE['microindel_size']}) were saved to sam file $FILTER_DIR/$3_mapped_MAPQ_XOXM.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: indels] OK mapped and MAPQ, XOXM filtered reads but with indels size (>${PARAMETERS_TABLE['microindel_size']}) were saved to excluded file $FILTER_DIR/$3_XOXM.excluded." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Count the number of reads in sam file filtered for indels greater than given size threshold
echo "$(date '+%Y%m%d %r') [Filtering: indels] $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam $(wc -l $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam | awk '{print $1}') " 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

echo "$(date '+%Y%m%d %r') [Filtering: indels] $FILTER_DIR/$3_INDELS.excluded $(wc -l $FILTER_DIR/$3_INDELS.excluded | awk '{print $1}') " 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

# Add the header at the end of the filtering
echo "$(date '+%Y%m%d %r') [Filtering: end] Appending filtered sam file to backuped sam header." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
cat $FILTER_DIR/header.txt $FILTER_DIR/$3_mapped_MAPQ_XOXM_INDELS.sam > tmp.sam 2>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: end] Failed Filtered reads sam file was not appended to sam header file. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filtering: end] Error: $(cat $ERROR_TMP)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: end] Filtered reads sam file was appended to sam header file in tmp sam file $FILTER_DIR/tmp.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi
mv tmp.sam $FILTER_DIR/$3_mapped_filtered.sam 2>>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
    echo "$(date '+%Y%m%d %r') [Filter end] Failed Tmp sam file was not renamed to $FILTER_DIR/$3_mapped_filtered.sam. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Filter end] Error: $(cat $ERROR_TMP)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit $rtrn
else
    echo "$(date '+%Y%m%d %r') [Filtering: begin conversion] OK Tmp sam file $FILTER_DIR/tmp.sam was renamed successfully to $FILTER_DIR/$3_mapped_filtered.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

echo "$(date '+%Y%m%d %r') [Filtering: begin conversion] OK Filtered reads appended to sam header file in $FILTER_DIR/$3_mapped_filtered.sam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1

# Convert SAM to BAM
if [[ -s $FILTER_DIR/$3_mapped_filtered.sam ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: samtools view] OK Starting sam conversion to bam format." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    samtools view \
	-$([[ ${PARAMETERS_TABLE["samtools_view_b"]} == "TRUE" ]] && echo -ne "b")$([[ ${PARAMETERS_TABLE["samtools_view_S"]} == "TRUE" ]] && echo -ne "S") \
	$FILTER_DIR/$3_mapped_filtered.sam > $FILTER_DIR/$3.bam 2>$FILTER_TMP\_1
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Filtering: samtools view] Failed Sam to bam conversion was not done on $FILTER_DIR/$3_mapped_filtered.sam. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: samtools view] Error: $(cat $FILTER_TMP\_1)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $FILTER_TMP\_1 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: samtools view] OK Conversion of sam file $FILTER_DIR/$3_mapped_filtered.sam to bam file $FILTER_DIR/$3.bam was done successfully." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
else
    echo "$(date '+%Y%m%d %r') [Filtering: samtools view] Failed Sam file $FILTER_DIR/$3_mapped_filtered.sam does not exist or is empty. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3" | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

# SORT BAM
if [[ -s $FILTER_DIR/$3.bam ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: samtools sort] OK Starting bam file sorting." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    samtools sort $FILTER_DIR/$3.bam $FILTER_DIR/$3_sorted 2>$FILTER_TMP\_2
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Filtering: samtools sort] Failed Bam file sorting was not done on bam file $FILTER_DIR/$3.bam. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: samtools sort] Error: $(cat $FILTER_TMP\_2)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $FILTER_TMP\_2 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: samtools sort] OK Bam file $FILTER_DIR/$3.bam was sorted successfully and saved to bam file $FILTER_DIR/$3_sorted.bam." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
else
    echo "$(date '+%Y%m%d %r') [Filtering: samtools sort] Failed Bam file $FILTER_DIR/$3.bam does not exist or is empty. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3" | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

 # INDEX BAM
if [[ -s $FILTER_DIR/$3_sorted.bam ]]; then
    echo "$(date '+%Y%m%d %r') [Filtering: samtools index] OK Starting bam file indexing." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    samtools index $FILTER_DIR/$3_sorted.bam 2>$FILTER_TMP\_3
    rtrn=$?
    if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Filtering: samtools index] Failed Bam file indexing was not done on bam file $FILTER_DIR/$3_sorted.bam. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: samtools index] Error: $(cat $FILTER_TMP\_3)" | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $FILTER_TMP\_3 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Filtering: samtools index] OK Bam file $FILTER_DIR/$3_sorted.bam was indexed successfully and saved to bam file $FILTER_DIR/$3_sorted.bam.bai." 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
else
    echo "$(date '+%Y%m%d %r') [Filtering: samtools index] Failed Bam file $FILTER_DIR/$3_sorted.bam does not exist or is empty. An error occured." | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 3" | tee -a $ERROR_TMP 2>&1 | tee -a $FILTER_DIR/$FILTER_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi

########################
# Section Analysis
########################

# Create Analysis directory if does not exist

if [[ -d $ANALYSIS_DIR ]]; then
    echo "$(date '+%Y%m%d %r') [Analysis] Start Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    echo "$(date '+%Y%m%d %r') [Analysis: directory] OK $ANALYSIS_DIR directory already exists. Will write Analysis output in this directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    mkdir $ANALYSIS_DIR 2>$ERROR_TMP
    if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: directory] Failed Analysis directory, $ANALYSIS_DIR, was not created." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code 126." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit 126
    else
	echo "$(date '+%Y%m%d %r') [Analysis] Start Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: directory] OK $ANALYSIS_DIR directory was created sucessfully. Will write Analysis output in this directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
fi

    # GENERATE BCF (PER SITE DEPTH, INFOS)

if [[ -s $FILTER_DIR/$3_sorted.bam ]]; then
    echo "$(date '+%Y%m%d %r') [Analysis] Sorted bam file exists and is not empty! Start to index fasta genome and compute per base depth ..." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1    
    #samtools faidx ${PARAMETERS_TABLE["REFERENCE_GENOME_FASTA"]} 2>$ANALYSIS_TMP\_1
    # This step will be done one time
    samtools mpileup \
	-$([[ ${PARAMETERS_TABLE["samtools_mpileup_B"]} -eq "TRUE" ]] && echo -ne "B") \
	-$([[ ${PARAMETERS_TABLE["samtools_mpileup_Q"]} ]] && echo -ne "Q" ${PARAMETERS_TABLE["samtools_mpileup_Q"]}) \
	-$([[ ${PARAMETERS_TABLE["samtools_mpileup_u"]} -eq "TRUE" ]] && echo -ne "u") \
	-$([[ ${PARAMETERS_TABLE["samtools_mpileup_f"]} -eq "TRUE" ]] && echo -ne "f") ${PARAMETERS_TABLE["SAMTOOLS_REFERENCE_GENOME_INDEX"]} $FILTER_DIR/$3_sorted.bam > $ANALYSIS_DIR/$3.bcf 2>$ERROR_TMP	
	rtrn=$?	
	if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: samtools mpileup] Failed formatage of $FILTER_DIR/$3_sorted.bam" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $ERROR_TMP | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
	echo "$(date '+%Y%m%d %r') [Analysis: samtools mpileup] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: samtools mpileup] OK Alignment was formated. Will write output in $ANALYSIS_DIR directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi

    #echo "$(date '+%Y%m%d %r') [Analysis: faidx]">> $LOG_DIR/$LOGFILE
    #cat $ANALYSIS_TMP\_1 >> $LOG_DIR/$LOGFILE
    echo "$(date '+%Y%m%d %r') [Analysis: samtools mpileup] $3.bcf" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y%m%d %r') [Analysis] Sorted bam file does not exist or is empty" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1   
    exit 3
fi

    # VARIANT CALLING

if [[ -s $ANALYSIS_DIR/$3.bcf ]]; then
    echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] bcf file exists and is not empty! Start to call variant" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1   
    bcftools view -$([[ ${PARAMETERS_TABLE["bcftools_view_c"]} -eq "TRUE" ]] && echo -ne "c")$([[ ${PARAMETERS_TABLE["bcftools_view_v"]} -eq "TRUE" ]] && echo -ne "v")$([[ ${PARAMETERS_TABLE["bcftools_view_g"]} -eq "TRUE" ]] && echo -ne "g") $ANALYSIS_DIR/$3.bcf > $ANALYSIS_DIR/$3.vcf 2>$ERROR_TMP	
	rtrn=$?
	if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] Failed variant calling" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $ERROR_TMP | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
	echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] Continue Analysis process" | tee -a $ANALYSIS_DIR/ANALYSIS_DIR_$_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1	
	echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] OK variant calling was done. Will write output in $ANALYSIS_DIR/ directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi
    echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] $3.vcf">> $LOG_DIR/$LOGFILE 2>&1
else
    echo "$(date '+%Y%m%d %r') [Analysis: bcftools view] bcf file does not exist or is empty" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1   
    exit 3
fi

    # VARIANT FILTERING
# TO DO (R)    
# vcfutils.pl varFilter $3.vcf > $3_filtre.vcf
    
    # COMPRESS and INDEX VCF file 

if [[ -s $ANALYSIS_DIR/$3.vcf ]]; then
    echo "$(date '+%Y%m%d %r') [Analysis] vcf file exists and is not empty! Start to compress it ..." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
    bgzip -$([[ ${PARAMETERS_TABLE["bgzip_f"]} -eq "TRUE" ]] && echo -ne "f")$([[ ${PARAMETERS_TABLE["bgzip_c"]} -eq "TRUE" ]] && echo -ne "c") $ANALYSIS_DIR/$3.vcf > $ANALYSIS_DIR/$3.vcf.gz 2>$ERROR_TMP	
	rtrn=$?	
	if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: bgzip] Failed compression of $3.vcf" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $ERROR_TMP | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
	echo "$(date '+%Y%m%d %r') [Analysis: bgzip] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: bgzip] OK $3.vcf was compressed. Will write output in $ANALYSIS_DIR directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    fi

    echo "$(date '+%Y%m%d %r') [Analysis: Compress vcf] $3.vcf.gz">> $LOG_DIR/$LOGFILE
    
    tabix -$([[ ${PARAMETERS_TABLE["tabix_p"]} ]] && echo -ne "p" ${PARAMETERS_TABLE["tabix_p"]}) -$([[ ${PARAMETERS_TABLE["tabix_f"]} -eq "TRUE" ]] && echo -ne "f") $ANALYSIS_DIR/$3.${PARAMETERS_TABLE["tabix_p"]}.gz 2>$ERROR_TMP	
	rtrn=$?		
	if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: tabix] Failed indexation" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $ERROR_TMP | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
	echo "$(date '+%Y%m%d %r') [Analysis: tabix] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: tabix] OK Compressed vcf was indexed. Will write output in $ANALYSIS_DIR directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi

    #echo "$(date '+%Y%m%d %r') [Analysis: tabix] (Index compressed ${PARAMETERS_TABLE['tabix_p']})] $3.${PARAMETERS_TABLE['tabix_p']}.gz" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 

else
    echo "$(date '+%Y%m%d %r') [Analysis: tabix] ${PARAMETERS_TABLE['tabix_p']} file does not exist" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
    exit 3
fi  
 
    # SNP EFFECT PREDICTION

echo "$(date '+%Y%m%d %r') [Analysis: snpEff] Start snpEff ..." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 

java -jar ${PARAMETERS_TABLE["SNPEFF_PATH"]}/snpEff.jar $GENOME_ALIAS \
    -c ${PARAMETERS_TABLE["SNPEFF_PATH"]}/snpEff.config \
	${PARAMETERS_TABLE["snpeff_options"]} \
    -i ${PARAMETERS_TABLE["snpeff_inFile_format"]} \
    -o vcf $ANALYSIS_DIR/$3.vcf > $ANALYSIS_DIR/$3_snpeff.vcf 2>$ERROR_TMP	
rtrn=$?
if [[ $? -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: snpEff] Failed snpEff" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
    else
	cat $ERROR_TMP | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1 
	echo "$(date '+%Y%m%d %r') [Analysis: snpEff] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: snpEff] OK $3_snpeff.vcf. Will write output in $ANALYSIS_DIR directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	fi

# SnpSift cmd: simplify Indels
echo "$(date '+%Y%m%d %r') [Analysis: SnpSift] Running SnpSift to simplify Indels in snpEff vcf output ..." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
java -jar ${PARAMETERS_TABLE["SNPEFF_PATH"]}/SnpSift.jar simplifyIndels $ANALYSIS_DIR/$3_snpeff.vcf > $ANALYSIS_DIR/$3_snpeff_snpsift.vcf 2>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: SnpSift] Failed SnpSift" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
else
	echo "$(date '+%Y%m%d %r') [Analysis: SnpSift] OK simplify indels in snpEff vcf output are in $3_snpeff_snpsift.vcf file" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: SnpSift] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Convert vcf to oneEffectPerLine
echo "$(date '+%Y%m%d %r') [Analysis: vcfEffOneLineN.pl] Running vcfEffOneLineN.pl script to parse one effect per line from snpsift output ..." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
perl ${PARAMETERS_TABLE["DREAMFILE_PATH"]}/vcfEffOneLineN.pl $ANALYSIS_DIR/$3_snpeff_snpsift.vcf 2>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: vcfEffOneLineN.pl] Failed vcfEffOneLineN.pl" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
else
	echo "$(date '+%Y%m%d %r') [Analysis: vcfEffOneLineN.pl] OK one effect per line in $3_snpeff_snpsift_OneLineEff.vcf file" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: vcfEffOneLineN.pl] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

# Generate the dream file
echo "$(date '+%Y%m%d %r') [Analysis: dreamFileMaker.py] Running dreamFileMaker.py script to reformat vcf input into dream file txt output ..." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
python ${PARAMETERS_TABLE["DREAMFILE_PATH"]}/dreamFileMaker.py $ANALYSIS_DIR/$3_snpeff_snpsift_OneLineEff.vcf 2>$ERROR_TMP
rtrn=$?
if [[ $rtrn -ne 0 ]]; then
	echo "$(date '+%Y%m%d %r') [Analysis: dreamFileMaker.py] Failed dreamFileMaker.py" | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] Exits the pipeline, with error code $rtrn." | tee -a $ERROR_TMP 2>&1 | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Pipeline error] More information can be found in $ERROR_TMP." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	exit $rtrn
else
	echo "$(date '+%Y%m%d %r') [Analysis: dreamFileMaker.py] OK reformated vcf output into $ANALYSIS_DIR/$3_snpeff_snpsift_OneLineEff_DF.txt file" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
	echo "$(date '+%Y%m%d %r') [Analysis: dreamFileMaker.py] Continue Analysis process" | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
fi

echo "$(date '+%Y%m%d %r') [Analysis: snpeff] Move genes and summary snpEff output files into Analysis directory." | tee -a $ANALYSIS_DIR/$ANALYSIS_DIR_$DATE.log 2>&1 | tee -a $LOG_DIR/$LOGFILE 2>&1
mv snpEff_* $ANALYSIS_DIR/.


#####################
# CLEANING
####################




####################
# Run Pipeline.Rnw
###################

# TODO: Pipeline.Rnw is a resource which should not be in the user working directory
# This resource file, like the other pipeline resources, should be placed in some location shared by all users: like /usr/local/share/<mutdetect_name_dir>
# Consider to move all template and other resource files: .Rnw, .docx, etc. to this shared location
# Do not forget to create the corresponding variable path: PIPELINE_SHARED=$PREFIX/share/<mutdetect_name_dir> with PREFIX=/usr/local
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

#exit 0
