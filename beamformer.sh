#!/bin/bash

#########################################################################
#                                                                       #
# Beamforming TBB data                                                  #
#                                                                       #
#########################################################################


#-----------------------------------------------------------------------#
# Input arguments                                                       #
#-----------------------------------------------------------------------#

# Default parameters
TBB_DATA=/data/projects/COM_ALERT/tbb
OUTDIR=/data/projects/COM_ALERT/pipeline/analysis/marazuela/data
SCRIPTS=/home/marazuela/scripts
LOBS=L597863

CPUS=28
POL=(0 1)
HBA=(HBA0 HBA1)
reftime="false"

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
  -t|--time|--time)          TIME=${2}
  shift; shift ;;
  -s|--src|--source)         SOURCE=${2:-$SOURCE}
  shift; shift ;;
  -r|--ra|--right-ascension) RA=${2}  # Not needed if source among known sources
  shift; shift ;;
  -d|--dec|--declination)    DEC=${2} # Not needed if source among known sources
  shift; shift ;;
  -dm|--dispersion-measure)  DM=${2}  # Not needed if source among known sources
  shift; shift ;;
  -rt|--reftime)             reftime=${2:-$reftime}
  shift; shift ;;
  -c|--cpus)                 CPUS=${2:-$CPUS}
  shift; shift ;;
esac
done

source $SCRIPTS/bash_fcts.sh

BF_DATA=$OUTDIR/${LOBS}_${TIME:0:16}
if [ ! -d $BF_DATA ]; then mkdir $BF_DATA; fi

#-----------------------------------------------------------------------#
# Main script                                                           #
#-----------------------------------------------------------------------#

###
# Checking source
###

# List of known sources
Crab=("Crab" "B0531+21" "B0531")
B0329=("B0329+54" "B0329")

# Defining coordinates (in radian) and DM
if [[ $B0329 =~ $SOURCE ]]; then
  echo "Source = B0329"
  RA=0.92934186635; DEC=0.952579228492; DM=26.7
elif [[ $Crab =~ $SOURCE ]]; then
  echo "Source = Crab"
  RA=1.45967506759; DEC=0.384224829441; DM=56.8
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

echo "filelist = ${filelist[@]}"
echo "stations = ${stations[@]}"
###
# Computing reference time
###

if [[ $reftime != "false" ]]; then 
  title1 "Computing reference time"
  python $SCRIPTS/ref_time.py --dm $DM \
      -f $BF_DATA ${filelist[@]}
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
    echo "python $SCRIPTS/bf_data.py $TBB_DATA/${LOBS}_${s}_${TIME}*tbb.h5 -f $OUTDIR -p $p -s HBA" >> $BF_DATA/process_${TIME}
  elif [[ $s == CS* ]]; then for h in ${HBA[@]}; do
      echo "python $SCRIPTS/bf_data.py $TBB_DATA/${LOBS}_${s}_${TIME}*tbb.h5 -f $OUTDIR -p $p -s $h" >> $BF_DATA/process_${TIME}
done; fi; done; done

# Running commands in parallel
echo "$SCRIPTS/parallel.sh $BF_DATA/process_${TIME} $CPUS"
bash $SCRIPTS/parallel.sh $BF_DATA/process_${TIME} $CPUS

#waitForFinish '[p]'ython

###
# Beamforming across stations
###

title1 "Beamforming across stations"

if [ ! -d $BF_DATA/bf_data ]; then mkdir $BF_DATA/bf_data; fi

for p in ${POL[@]}; do
  python $SCRIPTS/tbb_beamformer.py -f $BF_DATA/bf_data -b 8 \
    -i $BF_DATA/beam_data/${LOBS}_*_${TIME}*tbb_fft_*_pol${p}.beam
done

#########################################################################
# End of script
