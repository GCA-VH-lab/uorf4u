# Configuration file parameters


uorf4u configuration file allows detailed customization of the tool's parameters. Below, comments (*placed after ; in each line*) are allowed in configuration files. Here they are used to provide short parameter descriptions.


;[General]  
ncbi_genetic_code_name = Bacterial  
alternative_start_codons = False  
main_start_codon = ATG  
min_orf_length = 9  
sd_energy_cutoff = -2  
sd_window_length = 20  
blastp_evalue_cutoff = 1e-5  
upstream_region_length = 600  
minimal_upstream_region_length = 50  
orf_length_group_range = 15  
orfs_presence_cutoff = 0.3  
paths_identity_cutoff = 0.5  
max_number_of_assemblies = 30  
check_assembly_annotation = 1  

;[Output]  
sequences_to_write = nt, aa  

;[Paths]  
ref_energies = {internal}/energyRef-CCTCCT.json  
muscle_binary = {internal}/bin/muscle5.1.macos_arm64  
plot_msa_R_script = {internal}/msa_plot.R  
palette_nt = {internal}/palette_nt.txt  
palette_aa = {internal}/palette_aa.txt  

;[Output_dir_names]  
output_dir = uorf4u_{current_date}  

;[Pairwise alignment]  
alignment_type = aa  
global_match_score = 2  
global_mismatch_score = -1  
global_open_gap_score = -1  
global_extend_gap_score = -1  
global_target_end_gap_score = -1  
global_query_end_gap_score = -1  
alignment_score_cutoff = 0. 

;[Multiple sequence Alignment]  
consensus_threshold = 0.7  