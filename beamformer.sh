#!/bin/bash

#########################################################################
#                                                                       #
# Beamforming TBB data                                                  #
#                                                                       #
#########################################################################

#-----------------------------------------------------------------------#
# Input arguments                                                       #
#-----------------------------------------------------------------------#

progname=$0

# Default parameters
TBB_DATA=/data/projects/COM_ALERT/tbb
OUTDIR=/data/projects/COM_ALERT/pipeline/analysis/marazuela/data
SCRIPTS=/home/marazuela/scripts
LOBS=L597863

CPUS=28
POL=(0 1)
HBA=(HBA0 HBA1)
reftime="false"
TEST="false"

usage()
{
  echo -e "usage: bash $progname -t <time> -s <source> [-ra RA] [-dec DEC] [-dm DM] [-tbb TBB_DATA] [-o OUTDIR] [-sc SCRIPTS] [-l LOBS] [-c CPUS] [-rt] [-test]

  -t    --time                Timestamp of the frozen data
  -s    --src                 Source name. If source in R3, B0329, Crab, 
                              the RA, Dec and DM do not need to be provided.

  -ra   --right-ascension     Right ascension of the source in radians.
  -dec  --declination         Declination of the source in radians
  -dm   --dispersion-measure  Dispersion measure of the source.

  -tbb  --tbb-data            Path to directory where TBB data is stored.
                              Default: $TBB_DATA
  -o    --outdir              Path to output directory.
                              Default: $OUTDIR
  -sc   --scripts             Path to scripts.
                              Default: $SCRIPTS
  -l    --lobs                Observation project id.
                              Default: $LOBS
  -c    --cpus                Number of nodes to run processes in parallel.
                              Default: $CPUS

  -rt   --reftime             If provided, reference time is computed.
  -test --test                If provided, only two subbands will be beamformed.
  "
}

while [[ $# -gt 0 ]]; do
case "$1" in
  -tbb|--tbb-data)           TBB_DATA=${2:-$TBB_DATA}
  shift; shift ;;
  -o|--outdir|--output-dir)  OUTDIR=${2:-$OUTDIR}
  shift; shift ;;
  -sc|--scripts)             SCRIPTS=${2:-$SCRIPTS}
  shift; shift ;;
  -l|--lobs)                 LOBS=${2:-$LOBS}
  shift; shift ;;
  -t|--time)                 TIME=${2}
  shift; shift ;;
  -s|--src|--source)         SOURCE=${2:-$SOURCE}
  shift; shift ;;
  -ra|--right-ascension)     RA=${2}  # Not needed if source among known sources
  shift; shift ;;
  -dec|--declination)        DEC=${2} # Not needed if source among known sources
  shift; shift ;;
  -dm|--dispersion-measure)  DM=${2}  # Not needed if source among known sources
  shift; shift ;;
  -c|--cpus)                 CPUS=${2:-$CPUS}
  shift; shift ;;
  -rt|--reftime)             reftime="true"
  shift;;
  -test|--test)              TEST="true"
  shift;;
  -h|--help )                usage
  exit;;
  *)                         usage
  exit 1;;
esac
done

source $SCRIPTS/bash_fcts.sh

if [[ $TEST == "false" ]]; then 
  BF_DATA=$OUTDIR/${LOBS}_${TIME:0:16}
  t=""
elif [[ $TEST == "true" ]]; then 
  BF_DATA=$OUTDIR/${LOBS}_${TIME:0:16}_test
  t="-t"
fi
if [ ! -d $BF_DATA ]; then mkdir $BF_DATA; fi

parallel_process=$BF_DATA/process_${TIME}

#-----------------------------------------------------------------------#
# Main script                                                           #
#-----------------------------------------------------------------------#

###
# Checking source
###

# List of known sources
Crab=("Crab" "B0531+21" "B0531")
B0329=("B0329+54" "B0329")
R3=("R3" "FRB180916.J0158+65" "FRB180916")

# Defining coordinates (in radian) and DM
if [[ $B0329 =~ $SOURCE ]]; then
  echo "Source = B0329"
  RA=0.92934186635; DEC=0.952579228492; DM=26.7
elif [[ $Crab =~ $SOURCE ]]; then
  echo "Source = Crab"
  RA=1.45967506759; DEC=0.384224829441; DM=56.8
elif [[ $R3 =~ $SOURCE ]]; then
  echo "Source = R3"
  RA=0.51365039886; DEC=1.146686263660; DM=349.2
fi


title1 "BEAMFORMING $SOURCE"
echo -e "\tRA                = $RA"
echo -e "\tDEC               = $DEC"
echo -e "\tDM                = $DM"
echo -e "\tTIME              = $TIME"
echo -e "\tOBSID             = $LOBS"
echo -e "\tTBB DATA FOLDER   = $TBB_DATA"
echo -e "\tOUTPUT FOLDER     = $BF_DATA"
echo -e "\tSCRIPTS           = $SCRIPTS"
echo -e "\tCPUS              = $CPUS"

echo -e "\tTEST              = $TEST"
echo -e "\tCOMPUTING REFTIME = $reftime"

###
# Sorting station list
###

title1 "Sorting station list"

filelist=($(ls $TBB_DATA/${LOBS}_*_${TIME}*tbb.h5))
stationlist=()
for f in ${filelist[@]}; do
  s=$(h5ls $f | grep STATION | awk '{print $1}' | sed 's/STATION_//g')
  stationlist+=( ${s} )
done
stations=($(echo "${stationlist[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
IFS=$'\n' stations=($(sort <<<"${stations[*]}")); unset IFS

#echo "filelist = ${filelist[@]}"
echo "stations = ${stations[@]}"
###
# Computing reference time
###

if [[ $reftime != "false" ]]; then 
  title1 "Computing reference time"
#  python $SCRIPTS/ref_time.py --dm $DM \
#      -f $BF_DATA ${filelist[@]}
else
  title2 "Skipping reference time calculation"
fi

###
# Beamforming stations
###

title1 "Beamforming stations"

# Writing commands to run in parallel
if [ -f $BF_DATA/process_${TIME} ]; then rm $BF_DATA/process_${TIME}; fi

for s in ${stations[@]}; do for p in ${POL[@]}; do  
  if [[ $s == RS* ]]; then
    echo "python $SCRIPTS/bf_data.py $TBB_DATA/${LOBS}_${s}_${TIME}*tbb.h5 -f $OUTDIR -p $p -s HBA $t -r $RA -d $DEC --dm $DM" >> $parallel_process
  elif [[ $s == CS* ]]; then for h in ${HBA[@]}; do
      echo "python $SCRIPTS/bf_data.py $TBB_DATA/${LOBS}_${s}_${TIME}*tbb.h5 -f $OUTDIR -p $p -s $h $t -r $RA -d $DEC --dm $DM" >> $parallel_process
done; fi; done; done

# Running commands in parallel
echo "Commands written to $parallel_process"
bash $SCRIPTS/parallel.sh $parallel_process $CPUS

#waitForFinish '[p]'ython

###
# Beamforming across stations
###

title1 "Beamforming across stations"

if [ ! -d $BF_DATA/bf_data ]; then mkdir $BF_DATA/bf_data; fi

for p in ${POL[@]}; do
  echo "Starting tbb_beamformer"
  #python $SCRIPTS/tbb_beamformer.py -f $BF_DATA/bf_data -b 8 \
  #  -i $BF_DATA/beam_data/${LOBS}_*_${TIME}*tbb_fft_*_pol${p}.beam
done

#########################################################################
# End of script
