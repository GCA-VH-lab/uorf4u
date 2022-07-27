import shutil
from xml.etree import ElementTree
import Bio.Seq
import Bio.Align.AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio.Data import IUPACData
import Bio.Data.CodonTable
from Bio import Entrez
from Bio import Seq

import subprocess
import tempfile
import random
import json
import sys
import os

Entrez.email = "anonymous@mail.se"


class RefSeqProtein:
    """A RefSeqProtein object holds a RefSeq protein and information about it.

    Attributes:
        accession_number (str): RefSeq accession number.
        parameters (manager.Parameters): Parameters' class object.
        record (Bio.SeqRecord.SeqRecord): SeqRecord of the ncbi protein db. Can be obtained by the get_record() method.
        taxid (str): Taxid of the protein. Can be obtained with get_assemblies_coordinates() method.
        kingdom_taxid (str): Kingdom taxid of the protein. Can be obtained with get_assemblies_coordinates() method.
        organism (str): Organism name of the protein. Can be obtained with get_assemblies_coordinates() method.
        assemblies_coordinates (list): List of dictionaries with information about assemblies' coordinates of
        the protein obtained from ipg ncbi database.

    """

    def __init__(self, accession_number: str, parameters: manager.Parameters):
        """Create a RefSeqProtein object.

        Arguments:
            accession_number (str): RefSeq accession number.
            parameters (manager.Parameters): Parameters' class object.

        """
        self.accession_number = accession_number
        self.parameters = parameters
        self.record = None
        self.taxid = None
        self.kingdom_taxid = None
        self.organism = None
        self.assemblies_coordinates = None

    def get_record(self) -> Bio.SeqRecord.SeqRecord:
        """Get a SeqRecord object of the protein from the ncbi protein database.

        Note:
            This method returns a record and updates the record attribute.

        Returns:
            Bio.SeqRecord.SeqRecordRecord: Record of the protein.

        """
        handle = Entrez.efetch(db="protein", id=self.accession_number, rettype="gbwithparts", retmode="text")
        self.record = SeqIO.read(handle, "gb")
        return self.record

    def get_assemblies_coordinates(self) -> list:  # rename
        """Get assemblies coordinates of the protein.

        Note:
            This method returns a list of assemblies coordinates and updates the self.assemblies_coordinates attribute.

        Returns:
            list: List of dictionaries with information about assemblies' coordinates of the protein obtained
                from the ipg ncbi database.

        """
        handle = Entrez.efetch(db="protein", rettype="ipg", retmode="xml", id=self.accession_number)
        xml = (handle.read())  # .decode('utf-8')
        root = ElementTree.fromstring(xml)
        assemblies_coordinates = []
        for protein in root.iter("Protein"):
            if protein.attrib["source"] == "RefSeq":
                self.taxid = protein.attrib["taxid"]
                self.kingdom_taxid = protein.attrib["kingdom_taxid"]
                self.organism = protein.attrib["org"]
                for cds in protein.iter("CDS"):
                    if "assembly" not in cds.attrib.keys():
                        cds.attrib["assembly"] = "NA"
                    try:
                        assemblies_coordinates.append(dict(locus_id=cds.attrib["accver"],
                                                           start=(int(cds.attrib["start"]) - 1),
                                                           stop=int(cds.attrib["stop"]), strand=cds.attrib['strand'],
                                                           length=int(cds.attrib["stop"]) - (
                                                                   int(cds.attrib["start"]) - 1),
                                                           strain=cds.attrib["strain"], assembly=cds.attrib["assembly"],
                                                           org=cds.attrib["org"], taxid=cds.attrib["taxid"]))
                    except:
                        print(f"Attention: {cds.attrib} record is not completed and cannot be processed",
                              file=sys.stderr)

        if len(assemblies_coordinates) > 1:
            print(f"Warning message: {len(assemblies_coordinates)} assemblies were found for the protein "
                  f"{self.accession_number}. All assemblies will be included in the analysis by default.",
                  file=sys.stderr)
        if len(assemblies_coordinates) == 0:
            print(f"Warning message: {len(assemblies_coordinates)} assemblies were found for the protein "
                  f"{self.accession_number}. This protein record can be suppressed by ncbi.",
                  file=sys.stderr)

        self.assemblies_coordinates = assemblies_coordinates
        return assemblies_coordinates

    def blastp_searching_for_homologous(self) -> list:
        """Search for protein's homologous with blastp against refseq_protein database.

        Note:
            This function does not create a new object's attribute; It only returns a list of accession numbers.

        Returns:
            list: List of proteins' accession numbers obtained with blastp searching. This list contains the query
                protein's accession number.

        """
        handle = NCBIWWW.qblast("blastp", "refseq_protein", self.accession_number,
                                expect=self.parameters.arguments["blastp_evalue_cutoff"])
        xml = handle.read()
        hits_list = []
        hits_an_list = [self.accession_number]
        # xml = open('xml_blastp_output.xml')
        # xml = xml.read()
        root = ElementTree.fromstring(xml)
        query_length = int(root.find("BlastOutput_query-len").text)
        for hit in root.iter("Hit"):
            hit_id = hit.find("Hit_id").text.strip("ref").strip("|")
            if hit_id != self.accession_number:
                hit_description = hit.find("Hit_def").text
                subject_length = int(hit.find("Hit_len").text)
                hsp_identity_sum = 0
                hsp_positive_sum = 0
                hsp_align_length = 0
                for hsp in hit.iter("Hsp"):
                    hsp_identity_sum += int(hsp.find("Hsp_identity").text)
                    hsp_positive_sum += int(hsp.find("Hsp_positive").text)
                    hsp_align_length += int(hsp.find("Hsp_align-len").text)
                pident_to_query_length = hsp_identity_sum / query_length
                pident_to_seq_length = hsp_identity_sum / subject_length
                pident_to_alignment_length = hsp_identity_sum / hsp_align_length
                # ! pident values could be used for additional filters
                if hit_id not in hits_an_list:
                    hits_an_list.append(hit_id)

        return hits_an_list


class Homologous:
    """A Homologous object holds list of proteins homologous and information about them.

    Attributes:
        accession_numbers (str): RefSeq accession number.
        parameters (Parameters): Parameters' class object.
        records (list): list of RefSeqProtein objects of the proteins.
        upstream_sequences (list): List of SeqRecords objects of the proteins' genes' upstream sequences.
        codon_table (Bio.Data.CodonTable.CodonTable): Codon table (genetic code).
        orfs (dict): Dict with keys as upstream sequences' IDs and values as corresponding lists of ORF's objects.
        conserved_paths (dict): Dict with keys as lengths of ORFs and values as corresponding lists of Path's objects.
                (Path class holds list of ORFs from different upstream sequences and information about them).

    """

    def __init__(self, accession_numbers: list, parameters: manager.Parameters):
        """Create a Homologous object.

        Note:
            With initialisation it also creates the record attribute - a list of RefSeqProtein objects of the proteins
                based on accession numbers list.

        Arguments:
            accession_numbers (list): List of RefSeq accession numbers.
            parameters (manager.Parameters): Parameters' class object.

        """
        self.accession_numbers = accession_numbers
        self.parameters = parameters
        self.records = [RefSeqProtein(i, parameters) for i in accession_numbers]
        self.upstream_sequences = None
        self.codon_table = Bio.Data.CodonTable.unambiguous_dna_by_name[
            parameters.arguments["ncbi_genetic_code_name"]]
        self.orfs = None
        self.conserved_paths = None

    def get_upstream_sequences(self) -> list:
        """Get upstream sequences of the proteins' genes.

        Note:
            A protein may be found in several assemblies (for example in different strains). In such case
                __all__ assemblies will be used for getting upstream sequence.

        Returns:
            list: List of SeqRecords objects of the proteins' genes' upstream sequences.

        """
        upstream_sequences = []
        for record in self.records:
            record.get_assemblies_coordinates()
            record_upstream_sequences = []
            for assembly in record.assemblies_coordinates:
                handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="txt", id=assembly["locus_id"])
                locus_record = SeqIO.read(handle, "fasta")
                if assembly["strand"] == "+":
                    useq_start = max(0, assembly["start"] - self.parameters.arguments["upstream_region_length"])
                    useq_stop = assembly["start"]
                elif assembly["strand"] == "-":
                    useq_start = assembly["stop"]
                    useq_stop = min(len(locus_record.seq),
                                    assembly["stop"] + self.parameters.arguments["upstream_region_length"])
                useq_length = useq_stop - useq_start  # Add additional filtering by length!
                useq = locus_record.seq[useq_start:useq_stop]
                if assembly["strand"] == "-":
                    useq = useq.reverse_complement()
                if assembly["strain"] in assembly["org"]:
                    useq_name = assembly["org"]
                else:
                    useq_name = f"{assembly['org']} {assembly['strain']}"
                useq_record = SeqRecord(useq,
                                        id=f"{assembly['locus_id']}:{useq_start}:{useq_stop}:{assembly['strand']}",
                                        name=useq_name,
                                        description=f"{record.accession_number}, {assembly['org']}, "
                                                    f"strain: {assembly['strain']}, "
                                                    f"assembly: {assembly['assembly']}, length: {useq_length}")
                record_upstream_sequences.append(useq_record)
            upstream_sequences += record_upstream_sequences
            if len(record_upstream_sequences) == 0:
                print(f"Warning message: upstream sequences of {record.accession_number} cannot be annotated",
                      file=sys.stderr)
        self.upstream_sequences = upstream_sequences
        return self.upstream_sequences

    def annotate_orfs(self) -> None:
        """Annotate ORFs of the upstream sequences.

        Note:
            This function updates 'upstream_sequences' attribute.

        Returns:
            None

        """
        if self.upstream_sequences is None:
            raise manager.Ant4suorfError(f"Error: 'annotate_orfs()' method can't be called."
                                         f" The result of 'get_upstream_sequences()' method not found.")
        if self.parameters.arguments["alternative_start_codons"]:
            start_codons_list = self.codon_table.start_codons
        else:
            start_codons_list = [self.parameters.arguments["main_start_codon"]]
        orfs = dict()
        for useq in self.upstream_sequences:
            orfs[useq.id] = []
            for first_position in range((len(useq.seq) - 3) + 1):
                first_codon = useq.seq[first_position:first_position + 3]
                if first_codon.upper() in start_codons_list:
                    start_codon_position = first_position
                    for second_position in range(start_codon_position + 3, (len(useq.seq) - 3) + 1, 3):
                        second_codon = useq.seq[second_position:second_position + 3]
                        if second_codon.upper() in self.codon_table.stop_codons:
                            stop_codon_position = second_position
                            length = stop_codon_position - start_codon_position
                            id = f"{useq.name}:{len(useq.seq) - (start_codon_position + 1)}"
                            # id: organism:distance_from_the_start_codon_to_the_main_orf
                            sd_window_start = max(
                                [0, (start_codon_position - self.parameters.arguments["sd_window_length"])])
                            current_orf = ORF(parameters=self.parameters, id=id, start=start_codon_position,
                                              stop=stop_codon_position,
                                              nt_sequence=useq.seq[start_codon_position:stop_codon_position],
                                              sd_window_seq=useq.seq[sd_window_start:start_codon_position])
                            if current_orf.length >= self.parameters.arguments["min_orf_length"]:
                                orfs[useq.id].append(current_orf)
                                for annotated_orfs in orfs[useq.id]:
                                    if current_orf.stop == annotated_orfs.stop and \
                                            current_orf.id != annotated_orfs.id:
                                        current_orf.extended_orfs.append(annotated_orfs.id)
                            break
        self.orfs = orfs
        return None

    def filter_orfs_by_sd_annotation(self) -> None:
        """Filter annotated ORFs by presence Shine-Dalgarno sequence.

        Returns:
            None

        """
        for useq, orf_list in self.orfs.items():
            filtered_orf_list = []
            for orf in orf_list:
                orf.calculate_energies()
                if orf.min_energy < self.parameters.arguments["sd_energy_cutoff"]:
                    filtered_orf_list.append(orf)
            self.orfs[useq] = filtered_orf_list
        return None

    def save_annotated_orfs(self) -> None:
        """Save information about annotated ORFs as a set of tsv files.

        Note:
            tsv files will be saved to the subdir called 'annotated_ORFs' located in 'output_dir'.

        Returns:
            None

        """
        colnames = "\t".join(["id", "length", "nt_sequence", "aa_sequence", "sd_sequence_window", "extended_orfs"])
        if not os.path.exists(self.parameters.arguments["output_dir"]):
            os.mkdir(self.parameters.arguments["output_dir"])
        output_dir_path = os.path.join(self.parameters.arguments["output_dir"], "annotated_ORFs")
        if not os.path.exists(output_dir_path):
            os.mkdir(output_dir_path)
        for useq_id, orf_list in self.orfs.items():
            useq_name = [i.name for i in self.upstream_sequences if i.id == useq_id][0].replace(" ", "_")
            lines = [colnames]
            for orf in orf_list:
                lines.append("\t".join(
                    [orf.id, str(orf.length), str(orf.nt_sequence), str(orf.aa_sequence), str(orf.sd_window_seq_str),
                     ';'.join(orf.extended_orfs)]))
            with open(os.path.join(output_dir_path, useq_name + '.tsv'), 'w') as output:
                output.write("\n".join(lines))
        return None

    def conserved_orf_searching(self) -> dict:
        """Search for conserved orf in the upstream sequences.

        Note:
            It returns a dict with conserved ORFs and updates the self.conserved_paths attribute.

        Returns:
            dict: Dict with keys as lengths of ORFs and values as corresponding lists Path's objects.
                (Path class holds list of ORFs from different upstream sequences and information about them).

        """
        lengths = []
        for useq, orfs in self.orfs.items():
            for orf in orfs:
                lengths.append(orf.length)
        lengths = sorted(list(set(lengths)))

        global_aligner = Align.PairwiseAligner()
        global_aligner.mode = "global"
        global_aligner.match_score = self.parameters.arguments["global_match_score"]
        global_aligner.mismatch_score = self.parameters.arguments["global_mismatch_score"]
        global_aligner.open_gap_score = self.parameters.arguments["global_open_gap_score"]
        global_aligner.extend_gap_score = self.parameters.arguments["global_extend_gap_score"]
        global_aligner.target_end_gap_score = self.parameters.arguments["global_target_end_gap_score"]
        global_aligner.query_end_gap_score = self.parameters.arguments["global_query_end_gap_score"]
        length_variance = self.parameters.arguments["orf_length_group_range"]

        useqs = self.orfs.keys()
        conserved_paths = dict()

        for length in lengths:
            useqs_with_filtered_orfs = []
            filtered_orfs = dict()
            for useq in self.orfs.keys():
                filtered_orfs[useq] = []
                for orf in self.orfs[useq]:
                    if abs(length - orf.length) <= length_variance:
                        filtered_orfs[useq].append(orf)
                orfs_ids = [i.id for i in filtered_orfs[useq]]
                for orf in filtered_orfs[useq]:
                    if any(i in orf.extended_orfs for i in orfs_ids):
                        filtered_orfs[useq].remove(orf)
                if len(filtered_orfs[useq]) > 0:
                    useqs_with_filtered_orfs.append(useq)

            # print('NUMBER OF ASSEMBLIES:', len(useqs))
            # print(length, len(useqs_with_filtered_orfs) / len(useqs))
            if len(useqs_with_filtered_orfs) / len(useqs) > 0.3:  # add this cutoff to config
                conserved_paths[length] = []
                for initial_useq in filtered_orfs.keys():
                    for initial_orf in filtered_orfs[initial_useq]:
                        conserved_path = Path(self.parameters)
                        conserved_path.update(initial_orf)
                        for useq in random.sample(filtered_orfs.keys(), len(filtered_orfs.keys())):
                            if useq != initial_useq and filtered_orfs[useq] != []:
                                score_sums = []
                                for orf in filtered_orfs[useq]:
                                    score_sum = 0
                                    for path_orf in conserved_path.path:
                                        if self.parameters.arguments['type_of_alignment'] == 'nt':
                                            current_alignment = global_aligner.align(orf.nt_sequence,
                                                                                     path_orf.nt_sequence)
                                        elif self.parameters.arguments['type_of_alignment'] == 'aa':
                                            current_alignment = global_aligner.align(orf.aa_sequence,
                                                                                     path_orf.aa_sequence)
                                        score_sum += current_alignment.score
                                    score_sums.append(score_sum)
                                max_score = max(score_sums)
                                if max_score > self.parameters.arguments['alignment_score_cutoff']:
                                    if score_sums.count(max_score) == 1:
                                        selected_orf = filtered_orfs[useq][score_sums.index(max_score)]
                                    else:
                                        num_of_candidates = len(filtered_orfs[useq])
                                        highest_score_orfs = [filtered_orfs[useq][k] for k in
                                                              range(num_of_candidates)
                                                              if score_sums[k] == max_score]
                                        highest_score_orfs_length_dists = [orf_it.length - length for orf_it in
                                                                           highest_score_orfs]
                                        min_length_dist = min(highest_score_orfs_length_dists)
                                        if highest_score_orfs_length_dists.count(min_length_dist) == 1:
                                            selected_orf = highest_score_orfs[
                                                highest_score_orfs_length_dists.index(min_length_dist)]
                                        else:
                                            num_of_candidates = len(highest_score_orfs)
                                            the_closest_by_length_orfs = [highest_score_orfs[k] for k in
                                                                          range(num_of_candidates) if
                                                                          highest_score_orfs_length_dists[
                                                                              k] == min_length_dist]
                                            the_closest_by_length_orfs_lengths = [orf_it.length for orf_it in
                                                                                  the_closest_by_length_orfs]
                                            max_length = max(the_closest_by_length_orfs_lengths)
                                            selected_orf = the_closest_by_length_orfs[
                                                the_closest_by_length_orfs_lengths.index(max_length)]
                                    conserved_path.update(selected_orf, max_score)

                        if len(conserved_path) / len(filtered_orfs) >= 0.3:  # cutoff!
                            to_save_this_path = 1
                            for old_path in conserved_paths[length]:
                                fraction_of_identity = conserved_path.calculate_similarity(old_path)
                                if fraction_of_identity >= self.parameters.arguments["paths_identity_cutoff"]:
                                    if conserved_path.score > old_path.score:
                                        conserved_paths[length].remove(old_path)
                                    elif conserved_path.score <= old_path.score:
                                        to_save_this_path = 0
                            if to_save_this_path == 1:
                                conserved_paths[length].append(conserved_path)
        self.conserved_paths = conserved_paths
        return conserved_paths

    def get_msa_of_conseved_orfs(self) -> None:  # To finish
        """Run a multiple sequence alignment tool for each path object (set of conserved ORFs).

        Returns:
            None

        """
        for length, paths in self.conserved_paths.items():
            for path in paths:
                path.muscle_msa()

        return None

    def save_msa(self) -> None:
        """Save MSA of conserved ORFs as fasta files.

        Note:
            Fasta files will be saved to the subdirs: ['nucleotide_msa' - for MSA of nucleotide sequences of ORFs,
                'amino_acid_msa' - MSA of amino acid sequences of ORFs, and 'sd_msa' - MSA of SD sequence regions
                of ORFS). All of them located in 'output_dir'.

        Returns:
             None

        """
        if not os.path.exists(self.parameters.arguments["output_dir"]):
            os.mkdir(self.parameters.arguments["output_dir"])
        output_dirs = dict(nt=os.path.join(self.parameters.arguments["output_dir"], "nucleotide_msa"),
                           aa=os.path.join(self.parameters.arguments["output_dir"], "amino_acid_msa"),
                           sd=os.path.join(self.parameters.arguments["output_dir"], "sd_msa"))
        for key in output_dirs:
            if not (os.path.exists(output_dirs[key])):
                os.mkdir(output_dirs[key])
        for length, paths in self.conserved_paths.items():
            for i in range(len(paths)):
                path = paths[i]
                types = ['nt', 'aa', 'sd']
                for seq_type in types:
                    if seq_type == 'nt':
                        id = f"length-{length}.score–{round(path.score)}.{len(path.nt_msa)}_seqs.index_{i}.fa"
                        msa = path.nt_msa
                    elif seq_type == 'aa':
                        id = f"length-{round(length / 3)}.score–{round(path.score)}.{len(path.aa_msa)}-seqs.index_{i}.fa"
                        msa = path.aa_msa
                    elif seq_type == 'sd':
                        id = f"length-{round(length)}.score–{round(path.score)}.{len(path.aa_msa)}-seqs.index_{i}.fa"
                        msa = path.sd_msa
                    output = os.path.join(output_dirs[seq_type], id)
                    AlignIO.write(msa, output, "fasta")
        return None

    def plot_msa(self) -> None:
        """Plot MSA of conserved ORFs saved as fasta files.

        Note:
            R script based on ggmsa package [yulab-smu.top/ggmsa] used to produce MSA plots. R script (msa_plot.R)
                can be found in output_dir. This methods uses subprocess to run this R script in the following way:
                `Rscript output_dir/msa_plot.R --aa_msa path_to_aa_msa_ --nt_msa path_to_nt_msa
                --sd_msa path_to_sd_msa`. Since during each run of uorf4u a local copy of this script is created
                in your output_dir, you can change it without any consequences for next uorf4u runs.

        Returns:
            None

        """
        r_script_path = self.parameters.arguments['plot_msa_R_script']
        r_script_local = os.path.join(self.parameters.arguments["output_dir"], os.path.basename(r_script_path))
        shutil.copy(r_script_path, r_script_local)
        aa_msa_path = os.path.join(self.parameters.arguments["output_dir"], "amino_acid_msa")
        nt_msa_path = os.path.join(self.parameters.arguments["output_dir"], "nucleotide_msa")
        sd_msa_path = os.path.join(self.parameters.arguments["output_dir"], "sd_msa")
        subprocess.run(
            ["Rscript", r_script_local, "--aa_msa", aa_msa_path, "--nt_msa", nt_msa_path, "--sd_msa", sd_msa_path])


class ORF:
    """An ORF object holds information about an annotated ORF.

    Note:
        It's supposed that the ORFs class' objects will not be used directly by API users since
            it's only needed for other classes' methods.

    Attributes:
        parameters (manager.Parameters): Parameters' class object.
        id (str): identifier of the ORF. Format: organism (strain):distance_from_the_start_codon_to_the_proteins_orf
        sequence_id (str): identifier of the ORF's sequence (locus id from the ncbi database).
        start (int): start position of the ORF on the locus (0-based).
        stop (int): stop position of the ORF on the locus (0-based).
        length (int): ORF's nucleotide sequence length.
        nt_sequence (Bio.Seq.Seq): a Seq object of nucleotide sequence of the ORF.
        aa_sequence (Bio.Seq.Seq): a Seq object of amino acid sequence of the ORF.
        sd_window_seq (Bio.Seq.Seq): a Seq object of upstream sequence to the start codon of the ORF.
        min_energy (float): minimal value of thermodynamic interaction between aSD and putative SD sequences within the
            upstream sequences to the start codon.
        putative_sd_sequence (Bio.Seq.Seq): a Seq object of the putative SD sequence with the minimal energy value.
        extended_orfs (list): a list of ORFs with that are in frame with the ORF, but have upstream start codon.

    """

    def __init__(self, parameters: manager.Parameters, id: str, nt_sequence: Bio.Seq.Seq, sd_window_seq: Bio.Seq.Seq,
                 start: int, stop: int):
        """Create an ORF object.

        Arguments:
            parameters (manager.Parameters): Parameters' class object.
            id (str): identifier of the ORF. Format: locus_id:distance_from_the_start_codon_to_the_proteins_orf:length.
            nt_sequence (Bio.Seq.Seq): a Seq object of nucleotide sequence of the ORF.
            sd_window_seq (Bio.Seq.Seq): a Seq object of upstream sequence to the start codon of the ORF.
            start (int): start position of the ORF on the locus (0-based).
            stop (int): stop position of the ORF on the locus (0-based).

        """
        self.parameters = parameters
        codon_table = Bio.Data.CodonTable.unambiguous_dna_by_name[  # ambiguous can be needed!
            parameters.arguments['ncbi_genetic_code_name']]

        self.id = id
        self.sequence_id = id.split(':')[0]
        self.start = start
        self.stop = stop
        self.length = len(nt_sequence)
        self.nt_sequence = nt_sequence
        self.aa_sequence = self.nt_sequence.translate(table=codon_table)
        self.sd_window_seq = sd_window_seq
        self.extended_orfs = []
        self.min_energy = 0
        self.putative_sd_sequence = 'NA'

    def calculate_energies(self) -> None:
        """Calculate energies of putative SD sequences of the upstream sequence.

        Returns:
            None

        """
        # Loading reference energies json file
        with open(self.parameters.arguments['ref_energies']) as ref_energy_file:
            ref_energy = json.load(ref_energy_file)
        sd_seq_length = min([len(i) for i in ref_energy.keys()])
        # Energies calculations
        if len(self.sd_window_seq) >= min(ref_energy.values()):
            energies = []
            for position in range((len(self.sd_window_seq) - sd_seq_length) + 1):
                energies.append(
                    ref_energy[self.sd_window_seq[position:position + sd_seq_length]])
            if energies:
                self.min_energy = min(energies)
                if self.min_energy < self.parameters.arguments['sd_energy_cutoff']:
                    sd_start_position = energies.index(self.min_energy)  # Be careful, it could be more than one!
                    self.putative_sd_sequence = self.sd_window_seq[sd_start_position:sd_start_position + sd_seq_length]
                    self.sd_window_seq_str = (f"{self.sd_window_seq[0:sd_start_position].lower()}"
                                              f"{self.putative_sd_sequence.upper()}"
                                              f"{self.sd_window_seq[sd_start_position:sd_start_position + sd_seq_length:].lower()}")

        return None


class Path:
    """A Path object holds information about a list of conserved ORFs.

    Note:
        It's supposed that the Path class' objects will not be used directly by API users since
            it's only needed for other classes' methods.

    Attributes:
        parameters (manager.Parameters): Parameters' class object.
        path (list): List of the ORF class objects.
        score (float): Score of the Path (calculated as sum of pairwise alignments scores of ORFs).
        aa_msa (Bio.Align.MultipleSeqAlignment): Multiple sequence alignment (MSA) for amino acid sequences.
        aa_msa (Bio.Align.MultipleSeqAlignment): Multiple sequence alignment (MSA) for nucleotide sequences.
        sd_msa (Bio.Align.MultipleSeqAlignment): Multiple sequence alignment (MSA) for SD sequences (nt).
        aa_msa_consensus (Bio.Seq.Seq): Amino acid consensus sequence of the MSA.
        nt_msa_consensus (Bio.Seq.Seq): Nucleotide consensus sequence of the MSA.
        sd_msa_consensus (Bio.Seq.Seq): SD (nt) consensus sequence of the MSA.

    """

    def __init__(self, parameters: manager.Parameters):
        """Create a Path object.

        Arguments:
            parameters (manager.Parameters): Parameters' class object.

        """
        self.parameters = parameters
        self.path = []
        self.score = 0
        self.aa_msa = None
        self.nt_msa = None
        self.sd_msa = None
        self.nt_msa_consensus = None
        self.aa_msa_consensus = None
        self.sd_msa_consensus = None

    def update(self, orf: ORF, score=0):
        """Update a Path with a new ORF.

        Arguments:
            orf (ORF): an ORF class' object.
            score (float): a sum of pairwise alignment scores of the ORF against all ORFs in the Path.

        Returns:
            None

        """
        self.path.append(orf)
        self.score += score

    def __len__(self):
        """__len__ magic method for a Path object.

        Returns:
            int: length of the path attribute - a number of ORFs in a Path.

        """
        return len(self.path)

    def calculate_similarity(self, other) -> float:
        """Calculate fraction of identical ORFs between two Path object.

        __Note:__ If two objects have different length, the fraction will be calculated as a number of identical ORFs
                divided by length of the shortest Path.

        Returns:
            float: fraction of identical ORFs.

        """
        num_of_identical_elements = len(set(self.path) & set(other.path))
        fraction_of_identical_orfs = num_of_identical_elements / min(len(self), len(other))
        return fraction_of_identical_orfs

    def muscle_msa(self) -> None:
        """Run a multiple sequence alignment tool (muscle) for the ORFs nucleotide and amino acid sequences.

        Note:
            This method updates nt_msa and aa_msa attributes.

        Returns:
            None

        """
        types = ['nt', 'aa', 'sd']
        for seq_type in types:
            records = []
            for orf in self.path:
                if seq_type == 'nt':
                    record = SeqRecord(orf.nt_sequence, orf.id, f"length: {orf.length}", "")
                elif seq_type == 'aa':
                    record = SeqRecord(orf.aa_sequence, orf.id, f"length: {int(orf.length / 3)}", "")
                elif seq_type == 'sd':
                    record = SeqRecord(orf.sd_window_seq, orf.id, f"length: {len(orf.sd_window_seq)}", "")
                records.append(record)
            temp_input = tempfile.NamedTemporaryFile()
            SeqIO.write(records, temp_input.name, "fasta")
            temp_output = tempfile.NamedTemporaryFile()
            muscle = self.parameters.arguments["muscle_binary"]
            subprocess.run([muscle, "-align", temp_input.name, "-output", temp_output.name],
                           stderr=subprocess.DEVNULL)
            temp_input.close()
            msa = AlignIO.read(temp_output.name, "fasta")
            msa_info = Align.AlignInfo.SummaryInfo(msa)
            msa_consensus = msa_info.gap_consensus(threshold=self.parameters.arguments['consensus_threshold'])
            print(msa_consensus)
            temp_output.close()
            if seq_type == 'nt':
                self.nt_msa, self.nt_msa_consensus = msa, msa_consensus
            elif seq_type == 'aa':
                self.aa_msa, self.aa_msa_consensus = msa, msa_consensus
            elif seq_type == 'sd':
                self.sd_msa, self.sd_msa_consensus = msa, msa_consensus

        return None
