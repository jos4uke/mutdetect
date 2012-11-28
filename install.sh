#!/usr/bin/env bash

## Copyright (c) 2012 by Joseph Tran. This work is made available under the terms of the Creative Commons Attribution-ShareAlike 3.0 license, http://creativecommons.org/licenses/by-sa/3.0/.

################################################
##
## @script: install.sh
##
## @author: Joseph Tran (Joseph.Tran@versailles.inra.fr)
##
## @version: 0.0.1
##
## @date: 2012-11-15
##
## @description: This bash script performs the mutation detection pipeline deployment to the production server
##
###############################################


ARGV=("$@")
ARGC=$#
#echo ${ARGV[@]}
#echo $ARGC
#MUTDETECT_VERSION=$1
#echo ${ARGV[0]}
#echo ${ARGV[1]}

ERROR_LOG=/tmp/mutdetect_install_error.log
# PRINT_HELP=0
# PRINT_USAGE=0
# MUTDETECT_VERSION=""
# PREFIX="/usr/local"

# while [[ $ARGC -gt 0 ]]; do
#     case "$1" in
# 	-v|--version) MUTDETECT_VERSION=$(readlink -f "$2"); shift 2; ;;
# 	-p|--prefix) PREFIX=$(readlink -f "$2"); shift 2; ;;
# 	-h|--help) PRINT_HELP=1; shift ;;
# 	-u|--usage) PRINT_USAGE=1; shift ;;
# 	*) args+=("$1"); shift ;;
#     esac
# done

# eval set -- "${args[@]}"

# usage () {
#     echo "$0 $VERSION $PREFIX" 
#     cat <<USAGE
# $0 [-h|--help]|[-u|--usage]
# or
# $0 -v|--version <pipeline_git_tag_version> -p|--prefix <pipeline_prefix>
# USAGE
# }

# help () {
#     echo "./install.sh --version $VERSION --prefix $PREFIX" 
#     usage
#     cat <<HELP
#   DESCRIPTION:

#       This script fetches and installs mutdetect pipeline to
#       given directory prefix and for the given git tag version.

#   EXAMPLES:
#       ./install.sh v0.0.1 /usr/local

#   OPTIONS:
#       -v|--version <pipeline_git_tag_version>
#           will fetch the mutdetect pipeline from git repository
#           for the given git tag version <pipeline_git_tag_version>        

#       -p|--prefix <pipeline_prefix>
#           install pipeline under a given <pipeline_prefix>

#       -h|--help    - print this help and exit
#       -u|--usage   - print a short usage and exit
#   AUTHOR:
#       Copyright by Joseph.Tran@versailles.inra.fr
# HELP
# }

# if [ "$PRINT_HELP" = 1 ]; then help; exit; fi
# if [ "$PRINT_USAGE" = 1 ]; then usage; exit; fi

#if [[ -z $PREFIX && -z $MUTDETECT_VERSION ]]; then
if [[ $ARGC -ne 2 ]]; then
    echo -e "Incorrect number of required positional arguments. Should be equal to 2: \n- the first argument corresponding to the tag version of mutdetect pipeline in git repository, \n- and the second and last argument to the prefix used to deploy the pipeline." | tee $ERROR_LOG 2>&1
    #usage;
    exit 1
else
    echo -e "$(date '+%Y-%m-%d %H:%M:%S') Starting the mutdetect pipeline deployment" | tee $ERROR_LOG 2>&1
    echo "$(date '+%Y-%m-%d %H:%M:%S') Will deploy version tag: $MUTDETECT_VERSION" | tee -a $ERROR_LOG 2>&1
    MUTDETECT_VERSION=${ARGV[0]}
    echo "$(date '+%Y-%m-%d %H:%M:%S') Will deploy mutdetect pipeline in: $PREFIX" | tee -a $ERROR_LOG 2>&1
    PREFIX=${ARGV[1]}
fi 

# export mutdetect pipeline from git repo 
MUTDETECT_VERSION_DIR=/usr/local/archives/mutdetect/$MUTDETECT_VERSION
mkdir -p $MUTDETECT_VERSION_DIR 2>$ERROR_LOG
rtrn=$?
if [[ rtrn -ne 0 ]]; then
    echo -e "$MUTDETECT_VERSION_DIR directory creation failed with exit error code $rtrn." | tee -a $ERROR_LOG 2>&1
    echo -e "You can get more information about mutdetect pipeline installation in $ERROR_LOG." 2>&1
    exit $rtrn
else
    echo -e "$MUTDETECT_VERSION_DIR directory was created successfully." | tee -a $ERROR_LOG 2>&1
    echo -e "Will export mutdetect pipeline project from git repository." | tee -a $ERROR_LOG 2>&1
fi
git archive --prefix=mutdetect_$MUTDETECT_VERSION$PREFIX/ $MUTDETECT_VERSION | (cd $MUTDETECT_VERSION_DIR && tar xf -)
rtrn=$?
if [[ rtrn -ne 0 ]]; then
    echo -e "git archive failed to export mutdetect pipeline project with exit error code $rtrn." | tee -a $ERROR_LOG 2>&1
    echo -e "You can get more information about mutdetect pipeline installation in $ERROR_LOG." 2>&1
    exit $rtrn
else
    echo -e "$MUTDETECT_VERSION_DIR directory was created successfully." | tee -a $ERROR_LOG 2>&1
    echo -e "Will deploy mutdetect pipeline files to $PREFIX" | tee -a $ERROR_LOG 2>&1
fi

# deploy mutdetect pipeline
#cp -R * $PREFIX/. 2>&1
cd  $MUTDETECT_VERSION_DIR/mutdetect_$MUTDETECT_VERSION/ 2>$ERROR_LOG
rm -rf usr/local/test 2>$ERROR_LOG
for f in $(find . -name \* -print 2>$ERROR_LOG); do
    echo -e "cp $f in /$f" | tee -a $ERROR_LOG 2>&1
    cp --parents $f / 2>$ERROR_LOG
done

# set privileges
echo -ne "Setting executable privileges on mutdetect.sh script ..." | tee -a $ERROR_LOG 2>&1
chmod 755 $PREFIX/bin/mutdetect.sh 2>$ERROR_LOG
echo -e "OK" | tee -a $ERROR_LOG 2>&1

echo
echo "NOTE: You can get more information about mutdetect pipeline installation in $ERROR_LOG."
echo
echo
echo "************************************************************"
echo "*             INSTALLATION IS NOW COMPLETE                 *"
echo "************************************************************"
echo
echo
echo "You can now type mutdetect.sh to get usage help for the pipeline."
echo