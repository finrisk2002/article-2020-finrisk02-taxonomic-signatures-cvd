# readme

This directery contains ordination analysis (mainly Bray-Curtis-PCoA) scripts.

There is currently some duplicate code for computing PCAs. Some functions used in survival can be found in get_pca and add_pca_components_to_meta_data, in code/final/R_functions.R. Scripts that do essentially same thing in code/final/PCA/.

NB.

Unlike most other analyses, on FIMM Atlas server, the BC matrix and PCoA analysis is best computed as a grunjob on HPC nodes.

For this purpose, there is a separate bash script, rungrun_bray.sh, which schedules the grunjob with grun.py (with various not-so-pretty hacks).

