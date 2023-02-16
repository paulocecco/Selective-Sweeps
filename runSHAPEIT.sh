#!/bin/bash
for bedFile in $(ls *.bed)
do 
  echo ${bedFile%%.*} 
  ./shapeit2.sh ${bedFile%%.*} 
done
