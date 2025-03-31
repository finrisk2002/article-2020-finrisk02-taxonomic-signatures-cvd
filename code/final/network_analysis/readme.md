## Readme

installation instructions for SpiecEasi can be found at 
https://github.com/zdk123/SpiecEasi

compute_network_all.R 				# sources all of the below

set_up.R    						# load all libraries and set up output directories
compute_spieceasi_analyses.R 		# computes SpiecEasi analyses, saves results to datadir (set_up.R)
compute_sparcc_analysis.R 			# computes SparCC analyses, saves results to datadir (set_up.R)
compute_spieceasi_cox_figures.R 	# compute SpiecEasi subnet Cox analyses + create figures.
compute_sparcc_cox_figures.R 		# TODO

 # network analysis specific functions. TODO: move under aaro_functions/ for clarity?
pseq_spieceasi.R 					# run spieceasi on pseq
main_subnet_analysis_figures.R 		# main analysis and figures function
get_spiec_adj_matrix.R 				# get adjacency matrix from SpiecEasi result
*everything else*

aaro_functions/
 # helper functions used by Aaro, sometimes forks of code of Ville. -> sync with survival_analysis functions
create_cox_data.R