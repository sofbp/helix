#!/bin/bash

end=`cat input | wc -l`

for file in $(seq 1 1 $end);do
	seq=`head -$file input | tail -1`
	gnuplot<<EOF
set term png
set output "$seq"
plot "$seq\_GB" u 1:2 w l, "$seq\_GH" u 1:2 w l 
EOF
done
