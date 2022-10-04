# Configuration file parameters


uorf4u configuration file allows detailed customization of the tool's parameters. Below, comments (*placed after ; in each line*) are allowed in configuration files. Here they are used to provide short parameter descriptions.  

***Note:***   
uorf4u has two premade configuration files: *uorf4u_eukaryotes.cfg* and *uorf4u_prokaryotes* located in: ./uorftu/uorf4u_data/ folder (internal). By default, 'prokaryotes' config file is used if no path or name of premade file is specified by a cmd parameter: `-c  prokaryotes|eukaryotes|<file.cfg>`.


***;[General]***  
**ncbi_genetic_code_name = Bacterial**; *the ncbi genetic code name ('Standard' for eukaryotes' config)*   
**upstream_region_length = 1000**; *[int or 'all'] Length of upstream region to retrieve. 'all' value is set for eukaryotes config file since by default it uses only mRNAs sequences.*  
**minimal_upstream_region_length = 500**; *[int] minimal upstream region length for sequences to retrieve. If available sequence length to retrieve is shorter then this record won't be taken in the analysis*.    
**downstream_region_length = 100**; *[int] downstream region (overlapped with CDS) length to retrieve.*  
**filter_refseq_sequences_by_regex = True**; *[bool] use or not regex parameter (below) for filtering the NCBI RefSeq sequences to retrieve.*    
**refseq_sequences_regex = ^[N]._.***; *[regex] that will be used to filter the NCBI RefSeq sequnces. For eukaryotes set as '^[NX]M_.*' that means that only sequences that start with NM_ or XM_ (mRNAs) will be taken in the analysis.* 
**max_number_of_assemblies = 5**; *[int] max number of assemblies to take into analysis for each protein. If there are more sequences in the identical protein database then random sampling will be used.


**;[blastp homologous searching]**
**blastp_evalue_cutoff = 1e-5**; *[float] blastp e-value cutoff during the searching for homologs against the RefSeq database.*  
**blastp_hit_list_size = 500**; *[int] max number of blastp hits to take in the analysis.*  
**blastp_max_number_of_alignments = 1000**; *[int] max number of alignments during the searching (there could be several alignments for 1 hit, see blastp documentation)*    
**blastp_pident_to_query_length_cutoff = 0.5**; *[float: 0-1] cutoff for hit's identity to your query protein.*  

***;[ORF annotation]***  
**alternative_start_codons = False**; *[bool] use or not set of alternative start codons.*   
**main_start_codon = ATG**;  *[str]*        
**min_orf_length = 9**; *[int] cutoff for ORFs length during annotation*    
**filter_by_sd = True**; *[bool] filter annotated ORFs by Shine-Dalgarno sequence prersence. Filtering based on calculation of binding energy between aSD sequence (UCCUCC) and putative SD sequence in an upstream to uAUG window. Energy calculation performed as described here: [Yang et.al, 2016](10.1534/g3.116.032227)*  
**sd_energy_cutoff = -3**; *[float] cutoff for aSD-SD binding energy.*    
sd_window_length = 20  
**check_assembly_annotation = False**; [bool] retrieve or not the NCBI sequences annotation to be sure that annotated uORFs are not overlapped with known CDSs (can be overriden by '-annot' cmd parameter).

**;[conserved ORFs searching]**  
fast_searching = 0  
fast_searching_fraction_of_initial_genomes = 0.3  
orf_length_group_range = 0.25  
orfs_presence_cutoff = 0.5  
paths_identity_cutoff = 0.5  
max_num_of_initial_genome_iteration = 200  

**;[Pairwise alignment]**  
alignment_type = aa  
global_match_score = 2  
global_mismatch_score = -1  
global_open_gap_score = -1  
global_extend_gap_score = -1  
global_target_end_gap_score = -1  
global_query_end_gap_score = -1  
alignment_score_cutoff = 0  

**;[Multiple Sequence Alignment]**
msa_tool = maft  
consensus_threshold = 0.7  

**;[Paths]**  
ref_energies = {internal}/energyRef-CCTCCT.json  
muscle_binary = {internal}/bin/muscle5.1.macos_arm64  
maft_binary = {internal}/bin/mafft-mac/mafft.bat  
plot_msa_R_script = {internal}/msa_plot.R  
palette_nt = {internal}/palette_nt.txt  
palette_aa = {internal}/palette_aa.txt  

**;[Output]**
sequences_to_write = nt, aa, sd  
logo_type = probability  
output_dir = uorf4u_{current_date}  

;------------------------  
**;Annotation visualisation**  
;------------------------  
**;[General figure parameters]**  
margin = 0.1  
gap = 0.03  
label_gap = 0.07  
orf_height = 0.15  
annotation_width = auto  
mm_per_nt = 0.04  
font_regular = {internal}/fonts/Lato-Regular.ttf  
font_bold = {internal}/fonts/Lato-Bold.ttf  

**;[Sequence labels]**  
label_color = #3D3D3D  
label_color_alpha = 1  
label_height_to_orf_height = 0.65  

**;[Axis tics]**  
axis_tics_font_size = auto  
axis_tics_line_width = 0.3  

**;[Loci annotations]**  
upstream_seq_line_color = #CECECE  
upstream_seq_line_color_alpha = 1  
upstream_seq_line_width = 0.5  
cds_seq_stroke_color = #489143  
cds_seq_stroke_color_alpha = 0.8  
cds_seq_fill_color = #9ee19b  
cds_seq_fill_color_alpha = 0.03  
orf_line_width = 0.5  
conserved_uorfs_stroke_color = #4e4e4e  
conserved_uorfs_stroke_color_alpha = 1  
conserved_uorfs_fill_color = #ee8fb1  
conserved_uorfs_fill_color_alpha = 0.6  
other_uorfs_stroke_color = #CECECE  
other_uorfs_stroke_color_alpha = 1  
annotated_orf_stroke_color = #3d6f8e  
annotated_orf_stroke_color_alpha = 1  
