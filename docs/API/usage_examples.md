# Short example-drived guide to uorf4u API.  

uorf4u has a simple API allowing it programmatic usage from within a python program. Below we descrive several Python snippets that mimic results of command-line calls.



```python

import uorf4u

# loading data, parameters initialization
parameters = uorf4u.manager.Parameters()
parameters.load_config()

# Creating RefseqProtein class' object
refseq_protein = uorf4u.data_processing.RefSeqProtein(accession_number="#accession number", 
											          parameters=parameters)

# Searching for protein's homologous with blastp
homologous_list = refseq_protein.blastp_searching_for_homologous()

# Creating a Homologous class' object
homologous = uorf4u.data_processing.Homologous(homologous_list, parameters)

# Getting upstream sequences and annotating ORFs
homologous.get_upstream_sequences()
homologous.save_upstream_sequences()
homologous.annotate_orfs()
homologous.filter_orfs_by_sd_annotation()
homologous.save_annotated_orfs()

# Searching for conserved ORFs and saving the results
homologous.conserved_orf_searching()
homologous.filter_out_similar_paths()
homologous.run_msa()
homologous.save_msa()
homologous.plot_ggmsa_figs()
homologous.plot_logo_figs()
homologous.save_results_summary_table()

```
