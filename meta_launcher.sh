#!/bin/bash

generations=$1
max_iter=$2

total_iter=$generations

if [ -f "curr_gen.GAG" ]; then
	echo "Renaming 'curr_gen.GAG' file to 'OLD_curr_gen.GAG'"
	mv curr_gen.GAG OLD_curr_gen.GAG
fi

while [ "$generations" -gt 0 ]
do
	if [ "$generations" -gt "$max_iter" ]; then
		n_iter=$max_iter
		generations=$(($generations - $max_iter))
	else
		n_iter=$generations
		generations=0
	fi
	Rscript r_function/gen_algo.R $n_iter
done


echo "Total iterations"
echo $total_iter
echo "Iterations per script launch"
echo $n_iter
