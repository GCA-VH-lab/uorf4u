import data_processing
import manager
import methods


#table = Bio.Data.CodonTable.unambiguous_dna_by_name[parameters.config['general']['ncbi_genetic_code_name']]



#print(table.forward_table)
#print(table.back_table)
#print(table.start_codons)
#print(table.stop_codons)



'''


refseq_protein = data_processing.RefSeqProtein('WP_112844288.1', parameters)

homologous_list = refseq_protein.blastp_searching_for_homologous()
homologous = data_processing.Homologous(homologous_list, parameters)
homologous.get_upstream_sequences()

homologous.annotate_orfs()
homologous.filter_orfs_by_sd_annotation()
homologous.save_annotated_orfs()
homologous.conserved_orf_searching()
homologous.get_msa_of_conseved_orfs()
'''



'''`

input_upstream_seq = data_processing.Upstream_seqs(parameters)

input_upstream_seq.annotate_orfs()
input_upstream_seq.filter_orfs_by_sd_annotation()
input_upstream_seq.save_annotated_orfs()
paths = input_upstream_seq.conserved_orf_searching()



for length in paths.keys():
    print(length)
    for path in paths[length]:
        print(path.align())

'''
