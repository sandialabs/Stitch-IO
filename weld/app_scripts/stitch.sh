#!/bin/bash

# SPPARKS executable
#SPPARKS=$HOME/jaks.git/gitsync_spparks.cloned/src/spk_serial_stitch
SPPARKS=$HOME/jaks.git/gitsync_spparks.cloned/src/spk_flamer.gnu

# PYTHON
PYTHON=$HOME/local/python_3.5.2/bin/python3
STITCH_APP=stitch.weld.app_scripts.haz_box
STITCH_APP=haz_box

PREFIX=stitch
POTTS_IN=${PREFIX}_potts.in.template
WELD_IN=${PREFIX}_weld.in.template

# Initial simulation time
tn=0.0
# Initial pool position
yp_n=0.0

# weld_distance
weld_distance=1500.0
distance_traveled=0.0

#echo "       T        Distance Traveled"
#echo "     -----      -----------------"
#echo "$TN         $DISTANCE_TRAVELED"
STAGE=1
while (( $(bc<<<"$distance_traveled<$weld_distance") )); do

  STITCH_CMD="-istitch_weld.in --yp=${yp_n} --t0=${tn}";
  # First time is special (startup); need to signal stitch 
  if [ "1" -eq "${STAGE}" ]
  then
    STITCH_CMD="-istitch_weld.in --yp=${yp_n} --t0=${tn} --startup";
  fi

  # Run stitch script
  declare VARSTR=(`$PYTHON -m ${STITCH_APP} ${STITCH_CMD}`)

  SEED=$RANDOM
  POTTS_CMD="-var SEED ${SEED} -log potts.${STAGE}.log";
  SEED=$RANDOM
  WELD_CMD="-var SEED ${SEED} -log stitch_weld.${STAGE}.log";
  echo "======================================================================" >> log.bash
  echo "STAGE N = ${STAGE}" >> log.bash
  echo "STITCH command: ${STITCH_CMD}" >> log.bash
  echo "SPPARKS command: ${POTTS_CMD}" >> log.bash
  echo "SPPARKS command: ${WELD_CMD}" >> log.bash
  echo "======================================================================" >> log.bash
  echo "** stitch_potts.parameters **" >> log.bash
  cat stitch_potts.params.txt >> log.bash
  echo "** stitch_weld.parameters **" >> log.bash
  cat stitch_weld.params.txt >> log.bash
  echo "======================================================================" >> log.bash
  echo "" >> log.bash

  mpiexec -np 16 $SPPARKS $POTTS_CMD < $POTTS_IN
  #$SPPARKS $POTTS_CMD < $POTTS_IN

  mpiexec -np 16 $SPPARKS $WELD_CMD < $WELD_IN
  #$SPPARKS $WELD_CMD < $WELD_IN

  # Setup for next cv
  tn=${VARSTR[0]}
  yp_np1=${VARSTR[1]}

  # Some arithmetic to determine distance traveled
  delta=$(bc<<<"$yp_np1-$yp_n")
  distance_traveled=$(bc<<<"$distance_traveled+$delta")
  #echo "tn=${tn}; yp_np1=${yp_np1}; distance_traveled=${distance_traveled}"
  yp_n=$yp_np1
  let STAGE+=1
done

