#/bin/bash

./part $1 > /dev/null
./evalord $1 > orig
./evalord $1 order.txt > reord
gnuplot -e 'plot "orig", "reord" ; pause -1'
