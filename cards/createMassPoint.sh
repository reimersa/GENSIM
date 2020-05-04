#! /bin/bash
# Author: Izaak Neutelings (November, 2019)
# Description: Quickly create new cards with sed.

[[ $# -lt 2 ]] && M=$@ || M=500

for f in `ls *M1000*`; do
  cp $f `echo $f | sed "s/1000/$M/"`;
done

sed -i "s/M1000/M$M/" *_M${M}_proc_card.dat
sed -i "s/1000/$M/" *_M${M}_customizecards.dat

