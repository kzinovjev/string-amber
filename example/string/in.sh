#!/bin/bash -f

if [ -f string.groupfile ];
then
  rm string.groupfile
fi

nrep=$1
echo $nrep
count=0
for i in `seq 1 $nrep`
do
  REP=$i
  sed "s/XXXXX/$RANDOM/g" in > $REP.in
  if [ $i -le $(($nrep/2)) ];
  then
    echo "-O -rem 0 -i $REP.in -o $REP.out -c ../react.rst -r $REP.rst -x $REP.mdcrd -inf $REP.mdinfo -p ../ch3cl.prmtop" >> string.groupfile
  else
    echo "-O -rem 0 -i $REP.in -o $REP.out -c ../prod.rst -r $REP.rst -x $REP.mdcrd -inf $REP.mdinfo -p ../ch3cl.prmtop" >> string.groupfile
  fi
done

echo "N REPLICAS  = $nrep"
echo " Done."


