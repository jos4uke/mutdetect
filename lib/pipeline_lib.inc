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
    #old=$IFS
    #IFS=$'\n'
    #for line in $(cat $INFILE)
    while read line
    do
	#XM_XO=$(echo $line | perl -ne 'if ($_=~m/XM:i:(\d+).*XO:i:(\d+)/) {print $1+$2;}')
        #echo -e "XM+XO=$XM_XO" 1>&2

	# if [[ $XM_XO -le ${PARAMETERS_TABLE["nb_of_independent_event"]} ]]
	#     then
	#     echo $line
	# fi
	echo $line | gawk -v limit=${PARAMETERS_TABLE["nb_of_independent_event"]} 'match($0, /XM:i:([0-9]+).*XO:i:([0-9]+)/, a) {if ((a[1]+a[2])<=limit) {print $0}}'
    done < $INFILE
    #IFS=$old
}