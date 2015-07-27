#!/bin/bash

cd polybench
for x in *.c
do
   if [ $x != "polybench.c" ]
   then
      echo "generating ${x%.*}.gra..."
      ./gen_graphs.sh ${x%.*}
   fi
done
