#!/bin/bash

config_file=$1 # Path to config file, absolute or relative to launching dir.
max_iter=$2 # Maximum number of generations per GA script launch.
total_iter=$3 # Total number of generations to run.

remaining_gen=$total_iter

if [ -f "curr_gen.GAG" ]; then
	echo "Renaming 'curr_gen.GAG' file to 'OLD_curr_gen.GAG'"
	mv curr_gen.GAG OLD_curr_gen.GAG
	sleep 2
fi

while [ "$remaining_gen" -gt 0 ]
do
	if [ "$remaining_gen" -gt "$max_iter" ]; then
		n_iter=$max_iter
		remaining_gen=$(($remaining_gen - $max_iter))
	else
		n_iter=$remaining_gen
		remaining_gen=0
	fi
	count_iter=$(($total_iter - $remaining_gen - $n_iter))
	Rscript r_function/gen_algo.R $config_file $n_iter $count_iter
done


echo "Total iterations:"
echo $total_iter
echo "Maximum iterations per script launch:"
echo $max_iter
