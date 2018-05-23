#!/bin/bash

for k in {0..5}
do
	block=$(( 32*2**k ))
	(nvprof --metrics gld_efficiency,gld_throughput,achieved_occupancy ../bin/reduction $block) &> tmpfile.txt
	sed -n 1p tmpfile.txt >> ../measurements/efficiency.txt
	sed -n 10,19p tmpfile.txt >> ../measurements/efficiency.txt
	echo >> ../measurements/efficiency.txt
done
