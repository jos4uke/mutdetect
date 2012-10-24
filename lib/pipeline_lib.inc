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
    
    let count=0
    rtrn=0

    while read line
    do
	let count++
	key=$(echo $line | awk -F"=" '{print $1}');
	rtrn=$?
	if [[ $rtrn -ne 0 ]]; then
	    echo -e "$(date '+%Y%m%d %r') [$0] Failed Error parsing parameter key at line $count." >&2
	    return $rtrn
	else
	    if [[ ${key} =~ "#" || -z ${key} ]]
	    then
		continue
	    else
		value=$(echo $line | awk -F"=" '{print $2;}')
		rtrn=$?
		if [[ $rtrn -ne 0 ]]; then
		    echo -e "$(date '+%Y%m%d %r') [$0] Failed Error parsing parameter value at line $count." >&2
		    return $rtrn
		else
		    PARAMETERS_TABLE[${key}]=${value}
	            #echo -e "$key\t$value"
	            #echo -e "$key\t${PARAMETERS_TABLE[${key}]}"
		fi
	    fi  
	fi
    done < $INFILE

# TODO Check parameters interval and parameters type, if problems exit
    return $?
}

#---------------------------------------
# get_pipeline_user_parameters function
# 
get_pipeline_user_parameters()
{
    INFILE=$1

    let count=0
    rtrn=0

    while read line
    do
	let count++
	key=$(echo $line | awk -F"=" '{print $1}');
	rtrn=$?
	if [[ $rtrn -ne 0 ]]; then
	    echo -e "$(date '+%Y%m%d %r') [$0] Failed Error parsing parameter key at line $count." >&2
	    return $rtrn
	else
	    if [[ ${key} =~ "#" || -z ${key} ]]
	    then
		continue
	    else
		value=$(echo $line | awk -F"=" '{print $2;}')
		rtrn=$?
		if [[ $rtrn -ne 0 ]]; then
		    echo -e "$(date '+%Y%m%d %r') [$0] Failed Error parsing parameter value at line $count." >&2
		    return $rtrn
		else
		    if [[ -n ${PARAMETERS_TABLE[${key}]} ]]; then 
			echo -e "$(date '+%Y%m%d %r') [$0] OK Override ${key} parameter default value: ${PARAMETERS_TABLE[${key}]} by user defined value: ${value}" >&2
		    fi
		    PARAMETERS_TABLE[${key}]=${value}
	            #echo -e "$key\t$value"
	            #echo -e "$key\t${PARAMETERS_TABLE[${key}]}"
		fi
	    fi  
	fi
    done < $INFILE

# TODO Check parameters interval and parameters type, if problems exit
    return $?
}

#-----------------------------------------------------------
# remove_reads_with_more_than_x_independent_events function
#
remove_reads_with_more_than_x_independent_events()
{
    LIMIT=$1
    INFILE=$2
    #old=$IFS
    #IFS=$'\n'
    #for line in $(cat $INFILE)
    #while read line
    #do
	#XM_XO=$(echo $line | perl -ne 'if ($_=~m/XM:i:(\d+).*XO:i:(\d+)/) {print $1+$2;}')
        #echo -e "XM+XO=$XM_XO" 1>&2

	# if [[ $XM_XO -le ${PARAMETERS_TABLE["nb_of_independent_event"]} ]]
	#     then
	#     echo $line
	# fi
	#echo $line | gawk -v limit=${PARAMETERS_TABLE["nb_of_independent_event"]} 'match($0, /XM:i:([0-9]+).*XO:i:([0-9]+)/, a) {if ((a[1]+a[2])<=limit) {print $0}}'
    #done < $INFILE
    #IFS=$old

    gawk -v limit="${LIMIT}" 'match($0, /XM:i:([0-9]+).*XO:i:([0-9]+)|XO:i:([0-9]+).*XM:i:([0-9]+)/, a) {if ((a[1]+a[2])<=limit) {print $0} else {print $0 > "/dev/stderr"}}' $INFILE
}

#-------------------------------------------------------------
# remove_reads_with_indels_size_greater_than_y_bases function
#
remove_reads_with_indels_size_greater_than_y_bases()
{
    LIMIT=$1
    INFILE=$2
    #while read line 
    #do
	# find a way to return all occurences of the pattern in CIGAR code
	# echo '10M5I5D8M3I7M' | gawk 'function findall(str, re) { while(match(str, re, a)) {sub("I", "", a[0]); print a[0]; str = substr(str, RSTART+RLENGTH)} } { findall($0, "([0-9]+)I"); }'
	#echo $line | 
    gawk -v limit=$LIMIT \
	    'function findAllIndelsSize(str, lim)
             {
                flag=0;
                while(match(str, /([0-9]+)[ID]/, a)) 
                {
                   sub("[ID]", "", a[0]);
                   if (a[0] <= lim)
                   {
                      flag=1;
                   } else
                   {
                      flag=0;
                      break;
                   }
                   str = substr(str, RSTART+RLENGTH)
                }
                if (flag == 1)
                {
                   print $0;
                }  else
                {
                   print $0 > "/dev/stderr"
                }
             }
             {
                if ($6 ~ "^[0-9]+M$")
                {
                   print $0;
                } else if ($6 ~ "[ID]")
                {
                   findAllIndelsSize($0, limit)
                } else # reads with soft clipping
                {
                   print $0; 
                }
             }' $INFILE
    #done < $INFILE
}