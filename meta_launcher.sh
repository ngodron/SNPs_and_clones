#!/bin/bash

config_file=$1 # Path to config file, absolute or relative to launching dir.
max_iter=$2 # Maximum number of generations per GA script launch.
total_iter=$3 # Total number of generations to run.
save=$4 # 0,1 or 2, determines what outputs are saved afterwards
verbose=$5 # 0, 1 or 2, determines how much information is printed in console


remaining_gen=$total_iter

if [ -f "./output/curr_gen.GAG" ]; then
	echo "Renaming 'curr_gen.GAG' file to 'OLD_curr_gen.GAG'"
	mv ./output/curr_gen.GAG ./output/OLD_curr_gen.GAG
	sleep .5
fi

if [ -f "./output/all_generations.tsv" ]; then
	echo "Renaming 'all_generations.tsv' file to 'OLD_all_generations.tsv'"
	mv ./output/all_generations.tsv ./output/OLD_all_generations.tsv
	sleep .5
fi

if [ -f "./output/lastgen_models_list.GAG" ]; then
	echo "Renaming 'lastgen_models_list.GAG' file to 'OLD_lastgen_models_list.GAG'"
	mv ./output/lastgen_models_list.GAG ./output/OLD_lastgen_models_list.GAG
	sleep .5
fi

if [ -f "./output/all_gen_list.GAG" ]; then
	echo "Renaming 'all_gen_list.GAG' file to 'OLD_all_gen_list.GAG'"
	mv ./output/all_gen_list.GAG ./output/OLD_all_gen_list.GAG
	sleep .5
fi



printf "gen\tscor\tn_snps\n" > ./output/all_generations.tsv

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
	Rscript r_function/gen_algo.R $config_file $n_iter $count_iter $remaining_gen $save $verbose
	cat ./output/all_gen_temp >> ./output/all_generations.tsv
	rm ./output/all_gen_temp
done

echo "Total iterations:"
echo $total_iter
echo "Maximum iterations per script:"
echo $max_iter
