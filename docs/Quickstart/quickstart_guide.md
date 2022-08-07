# Quickstart guide

Here we present several examples of uorf4u usage and the respective command-line parameters.  
This chapter based on the considering of well-known uORFs that can be cis-acting translation modulators (see a review article [Koreaki Ito et.al. 2013](https://www.annualreviews.org/doi/10.1146/annurev-biochem-080211-105026)).

Before start, the necessary sample data as well as adjustable tool' configuration files are priveded by uorf4u at the post-install step:    
`uorf4u --data`   


## ErmCL 

One of the well-known bacterial upstream ORFs that regulates its downstream frame is ErmCL. Inducible expression of the downstream Erm resistance gene relies on ribosome stalling on the uORF (ErmCL). Molecular mechanism of ribosome stalling based on presence so-called ribosome arrest peptide (RAP) in ErmCL amino acid sequence. RAPs act in nascent states through interaction with ribosome that leads to translation arrest.  
RAP-mediated regulations based on presence a particular amino acid sequence (also known as arrest-essential amino acids). For ErmCL (19 codons length) this sequence is **IFVI** (see the [review](https://www.annualreviews.org/doi/10.1146/annurev-biochem-080211-105026) for more detailed introduction). 

To test whether uorf4u will be able to find this ORF we can use only accession number of ErmC protein as input! The searching result for "ErmC" in ncbi protein database gives us the RefSeq id: *WP_001003263.1*, which can be directly used for our searching:
`uorf4u -an WP_001003263.1 -verbose -o ErmC `  
*Note:* `-verbose` and `-o` parameters are optional. `-verbose` used to show all progress messages, `-o` - to specify output folder name (by default it's uorf4u_{current_data} e.g. uorf4u_2022_08_09-15_00).  

The results will be saved to the ErmC folder with the following structure:

<img  src="/img/output.svg" width="700"/>

Searching with default parameters returns us only one set of conserved ORFs. Let's have a look at the respective amino acid sequence logo. Fortunately, we can see that the most conserved region of the sequence is expected IFVI arrest-essential amino acids. 🥳

<img  src="/img/ermcl_logo.svg" width="700"/>



---

