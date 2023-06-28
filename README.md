# SNPs_and_clones

![Genetic Algorithm implementation visual explanation](GAGexplained.PNG)

# Goal

Implement a Genetic Algorithm to study SNP haplotypes.

# Grand Master Plan

- Better population diversity metrics 
--> (RGB colors with softmax-ed values of first three PCA axis?)
- Implement (again) min_mut in evolve.R?

# Scripts

### meta_launcher.sh
**Status:** Functional
Script to launch numerous generations, acts as a wrapper for

To launch, write the following command:
./meta_launcher.sh $1 $2 $3 $4 $5

Where:
config_file=$1 # Path to config file, absolute or relative to launching dir.
max_iter=$2 # Maximum number of generations per GA script launch.
total_iter=$3 # Total number of generations to run.
save=$4 # 0,1 or 2, determines what outputs are saved afterwards
verbose=$5 # 0, 1 or 2, determines how much information is printed in console

### gen_algo.R
**Status:** Functional
(*de facto* main script)

### load_inputs.R
**Status:** Missing weights

To Do:
- Implement weights
- Specify input files characteristics

### quality_control.R
**Status:** To be done

To Do:
- Implement covariates QC

### generate_G0.R
**Status:** Functional

This function is used to generate a starting generation, and if n_nov > 0 to insert random genomes at each generation before crossing-overs.

### calc_score.R
**Status:** Functional with rpart tree

GLM mock function for testing purposes.
Rpart tree accuracy score functional.

### evolve.R
**Status:** Functional

To Do:
- Implement unique elites to ensure diversity of elites.

#### Mutation mechanisms cull feature inflation.

Indeed, mutations occurr symmetrically (similar counts of feature gain / loss) until a threshold is met.
This unsures run time is constant and predictable between Gen 0 and Gen n_iter.

Also, results are ranked in order of minimal score first, features count second. This ensures that in the case where a plateau is reached by the optimal solution, unused features (i.e. not in the decision tree) are slowly purged from the elites by random loss mutations.

#### Crossing-over mechanisms induce the introduction of diversity in genomes.

Indeed, a rare event where a genome is shuffled fully was implemented to lower the probability of the whole generation getting stuck in the neighbourhood of a local optimum.

Right now, this event happens with a probability of 0.01 times the crossing-over rate (cr).


## Known issues
Memory fills up with unattainable R objects:
- Implemented quickfix: Dump and source genomes, outputs & parameters (meta_launcher.sh)
