#!/bin/bash

(nvprof ../bin/reduction 32) &> tmpfile.txt
sed -n 6p tmpfile.txt >> ../measurements/times.txt
sed -n 9p tmpfile.txt >> ../measurements/times.txt	
sed -n 11p tmpfile.txt >> ../measurements/times.txt	
echo >> ../measurements/times.txt

for k in {1..5}
do
	block=$(( 32*2**k ))
	(nvprof ../bin/reduction $block) &> tmpfile.txt
	sed -n 6p tmpfile.txt >> ../measurements/times.txt
	sed -n 10,11p tmpfile.txt >> ../measurements/times.txt	
	echo >> ../measurements/times.txt
done


