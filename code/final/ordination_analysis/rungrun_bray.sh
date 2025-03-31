#!/usr/bin/bash

# schedule Bray-Curtis matrix computation with grun.py

grun.py -q highmem.q -n brayc -c "bash code/final/grun_util/run_rscript_pathfix.sh \
code/final/ordination_analysis/compute_braycurtis_species.R"


