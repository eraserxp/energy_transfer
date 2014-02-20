#!/bin/bash
#ls p_?.dat | xargs --max-args=1 awk 'NR==2 {max=$3} {if ($3>max) max=$3} END{print max}' > max.txt

#ls p_??.dat | xargs --max-args=1 awk 'NR==2 {max=$3} {if ($3>max) max=$3} END{print max}' >> max.txt

#ls p_???.dat | xargs --max-args=1 awk 'NR==2 {max=$3} {if ($3>max) max=$3} END{print max}' >> max.txt

ls p_{0..400}.dat | xargs --max-args=1 sed -n '/  *0  *0 /p' | awk '{ print $ 3}' > max.txt

awk 'NR==1 {max=$1; s0=$1} {if ($1>max) max=$1} \
END{print "The maximum probability is enhanced by "max/s0}' max.txt

pmax=`awk 'NR==1 {max=$1} {if ($1>max) max=$1} END{print max}' max.txt`
echo "The maximum probability is: " $pmax

ntime=`grep -n $pmax max.txt | awk 'BEGIN {FS=":"} {print $1}'`
echo 'which appear at time = ' $(($ntime-1))

xmgrace max.txt
