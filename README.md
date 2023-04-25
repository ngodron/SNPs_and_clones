# SNPs_and_clones

# Goal

Implement a Genetic Algorithm to study SNP haplotypes.

# Grand Master Plan

- Implement covariates, features always set at 1 (e.g. MTB lineages).
- Objective function with classification trees.
- Better diversity metrics.
- Non-fixed break point for crossing-over (picked at random in uniform distribution).
- Implement (again) min_mut in evolve.R

# Scripts

### load_inputs.R
**Status:** Prototype

### quality_control.R
**Status:** Halfway there

### generate_G0.R
**Status:** Functional

### calc_score.R
**Status:** Works with GLM

### evolve.R
**Status:** Functional

#### Mutation mechanisms cull feature inflation.

Indeed, mutations occurr symmetrically (similar counts of feature gain / loss) until a threshold is met.
This unsures run time is constant and predictable between Gen 0 and Gen n_iter.


#### Crossing-over mechanisms induce the introduction of diversity in genomes.

Indeed, a rare event where a genome is shuffled fully was implemented to lower the probability of the whole generation getting stuck in the neighbourhood of a local optimum.

Right now, this event happens with a probability of 0.01 times the crossing-over rate (cr).

### gen_algo.R
**Status:** Reorganize temp.R (*de facto* main script)
