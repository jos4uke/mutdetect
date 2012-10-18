#! /bin/bash

# Function to 

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
