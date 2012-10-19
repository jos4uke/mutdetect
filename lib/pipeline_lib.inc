#! /bin/bash
#
# PIPELINE FUNCTIONS LIBRARY
#

#------------------------------------------
# get_pipeline_default_parameters function
# 
get_pipeline_default_parameters()
{
    INFILE=$1
    #declare -A PARAMETERS_TABLE
    # for line in $(cat $INFILE); do
    # 	read -a p_array <<< $(echo -e $line | tr "=" "\n")
    # 	PARAMETERS_TABLE[p_array[1]]=p_array[2]
    # 	echo -e "$p_array[1] \t $p_array[2]"
    # done 
    while read line
    do
	key=$(echo "$line" | awk -F"=" '{print $1}');
	if [[ ${key} =~ "#" || -z ${key} ]]
	    then
	    continue
	    else
	    value=$(echo "$line" | awk -F"=" '{print $2;}')
	    PARAMETERS_TABLE[${key}]=${value}
	    #echo -e "$key\t$value"
	    #echo -e "$key\t${PARAMETERS_TABLE[${key}]}"

	fi    
    done < $INFILE
}

#-----------------------------------------------------------
# remove_reads_with_more_than_x_independent_events function
#
remove_reads_with_more_than_x_independent_events()
{
    INFILE=$1
    while read line
    do
	XM=$(echo -e "$line" | perl -ne 'if ($_=~m/XM:i:(\d)+/) {print $1;}')
	#echo -e "XM=$XM" 1>&2
	XO=$(echo -e "$line" | perl -ne 'if ($_=~m/XO:i:(\d)+/) {print $1;}')
	#echo -e "XO=$XO" 1>&2
	if [[ "$[$XM + $XO]" -le ${PARAMETERS_TABLE["nb_of_independent_event"]} ]]
	    then
	    echo -e "$line"
	fi
    done < $INFILE
}