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

# Searching for protein's homologues with blastp
homologues_list = refseq_protein.blastp_searching_for_homologues()

# Creating a Homologues class' object
   homologues = uorf4u.data_processing.Homologues(homologues_list, parameters)

# Getting upstream sequences and annotating ORFs
homologues.get_upstream_sequences()
homologues.save_upstream_sequences()
homologues.annotate_orfs()
homologues.filter_orfs_by_sd_annotation()  
homologues.save_annotated_orfs()

# Searching for conserved ORFs and saving the results
homologues.conserved_orf_searching()
homologues.filter_out_similar_paths()
homologues.run_msa()
homologues.save_orfs_sequences()
homologues.save_msa()
homologues.save_results_summary_table()
homologues.plot_annotation()
homologues.plot_logo_figs()
homologues.plot_ggmsa_figs()

```
