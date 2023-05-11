# SNPs_and_clones

# Goal

Implement a Genetic Algorithm to study SNP haplotypes.

# Grand Master Plan

- Implement covariates, features to be one-hot encoded (e.g. MTB lineages).
- Better population diversity metrics.
- Implement (again) min_mut in evolve.R?

# Scripts

### load_inputs.R
**Status:** Prototype

To Do:
- Implement covariates loading
- Specify input files characteristics

### quality_control.R
**Status:** Halfway there

To Do:
- Implement covariates QC

### generate_G0.R
**Status:** Functional

This function is used to generate a starting generation, and if n_nov > 0 to insert random genomes before a crossing-over occurs.

### calc_score.R
**Status:** Functional

GLM mock function for testing purposes.
Rpart tree accuracy score functional.

### evolve.R
**Status:** Functional

#### Mutation mechanisms cull feature inflation.

Indeed, mutations occurr symmetrically (similar counts of feature gain / loss) until a threshold is met.
This unsures run time is constant and predictable between Gen 0 and Gen n_iter.

Also, results are ranked in order of minimal score first, features count second. This ensures that in the case where a plateau is reached by the optimal solution, unused features (i.e. not in the decision tree) are slowly purged from the elites by random loss mutations.

#### Crossing-over mechanisms induce the introduction of diversity in genomes.

Indeed, a rare event where a genome is shuffled fully was implemented to lower the probability of the whole generation getting stuck in the neighbourhood of a local optimum.

Right now, this event happens with a probability of 0.01 times the crossing-over rate (cr).

### gen_algo.R
**Status:** Reorganize temp.R to be launchable from terminal
(*de facto* main script)

## Known issues
Memory fills up with unattainable R objects (probably due to rbind):
- Potential quickfix: Dump and source genomes & parameters.
- Potential slowfix: Use only lists and vectors, no dataframe.
