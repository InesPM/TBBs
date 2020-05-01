#!/bin/bash

#-----------------------------------------------------------------------#
# Defining useful bash functions                                        #
#-----------------------------------------------------------------------#

title1(){
  message=$1
	i=0; x='===='
	while [[ i -lt ${#message} ]]; do x='='$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x\n"
}

title2(){
  message="$1 $observation"
	i=0; x='----'
	while [[ i -lt ${#message} ]]; do x='-'$x; ((++i)); done
  echo -e "\n\t  $message \n\t$x"
}

title3(){
  message="$1 $observation"
  echo -e "\n # $message"
}

function waitForFinish()
{
    local STRING;
    STRING=$3;

    # wait until jobs have started
    sleep 1

    # check, twice, whether all done
    for i in 1 2 ; do
	job=99
	while [ $job -gt 0 ] ; do sleep 10; top -b | head -n 40; job=`ps ax | grep ${STRING} | wc -l`; done
    done
}
