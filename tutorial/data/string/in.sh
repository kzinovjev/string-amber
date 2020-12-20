#!/bin/bash

# If groupfile exists, remove it
if [ -f string.groupfile ];
then
  rm string.groupfile
fi

PARM=../system/ch3cl.parm7  # path to topology file
REACT=../system/react.ncrst   # path to reactants structure
PROD=../system/prod.ncrst     # path to products structure
SEED=12345                  # Some random number 

# Number of string nodes is provided as command line argument
nodes=$1
for i in `seq 1 $nodes`; do
  # generate sander input for node $i 
  sed "s/__SEED__/$((SEED + i))/g" in > $i.in
  
  # use reactants structure for the first half of the nodes and products for the rest
  if [ $i -le $(($nodes/2)) ];
  then
    crd=$REACT
  else
    crd=$PROD
  fi
  
  # write out the sander command for node $i to the groupfile
  echo "-O -rem 0 -i $i.in -o $i.out -c $crd -r $i.ncrst -x $i.nc -inf $i.mdinfo -p $PARM" >> string.groupfile
done
