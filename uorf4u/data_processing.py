"""
This module provides data processing including uORFs annotation and conserved subset searching.
"""
import shutil
import xml.etree.ElementTree
import xml.etree.cElementTree
import Bio.Seq
import Bio.Align.AlignInfo
import Bio.Blast.NCBIWWW
import Bio.SeqRecord
import Bio.Entrez
import Bio.SeqIO
import Bio.Align
import Bio.AlignIO
import Bio.Data.IUPACData
import Bio.Data.CodonTable
import logomaker
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('agg')
import pandas
import subprocess
import tempfile
import statistics
import random
import math
import json
import sys
import re
import os

import uorf4u.manager
import uorf4u.drawing_annotation
import uorf4u.methods
import uorf4u.drawing_msa


class RefSeqProtein:
    """A RefSeqProtein object holds a RefSeq protein and information about it.

    Attributes:
        accession_number (str): RefSeq accession number.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        record (Bio.SeqRecord.SeqRecord): SeqRecord of the ncbi protein db. Can be obtained by the get_record() method.
        taxid (str): Taxid of the protein. Can be obtained with get_assemblies() method.
        kingdom_taxid (str): Kingdom taxid of a protein. Can be obtained with get_assemblies() method.
        organism (str): Organism name of a protein. Can be obtained with get_assemblies() method.
        name (str): Protein's product name from the ncbi (if available).
        assemblies_coordinates (list): List of dictionaries with information about assemblies' coordinates of
            the protein obtained from ipg ncbi database.
        loci (dict): Dict with keys as locus_ids and values as Locus class' objects.

    """

    def __init__(self, accession_number: str, parameters: uorf4u.manager.Parameters):
        """Create a RefSeqProtein object.

        Arguments:
            accession_number (str): RefSeq accession number.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        self.accession_number = accession_number
        self.name = "NA"
        self.parameters = parameters
        self.record = None
        self.taxid = None
        self.kingdom_taxid = None
        self.organism = None
        self.assemblies_coordinates = None
        self.loci = None

    def get_record(self) -> Bio.SeqRecord.SeqRecord:
        """Get a SeqRecord object of a protein from the ncbi protein database.

        Note:
            This method returns a record and updates the record attribute.

        Returns:
            Bio.SeqRecord.SeqRecordRecord: Record of the protein.

        """
        try:
            handle = Bio.Entrez.efetch(db="protein", id=self.accession_number, rettype="fasta", retmode="text")
            self.record = Bio.SeqIO.read(handle, "fasta")
            return self.record
        except Exception as error:
            raise uorf4u.manager.uORF4uError(
                "Unable to get a SeqRecord of the protein from the ncbi protein database.") from error

    def get_assemblies(self, xml_output=None) -> list:
        """Get assemblies (loci) coordinates of a protein.

        Note:
            This method returns a list of assemblies coordinates and updates the self.assemblies_coordinates attribute.

        Returns:
            list: List of dictionaries with information about assemblies' coordinates of a protein obtained
                from the ipg ncbi database.

        """
        try:
            if not xml_output:
                handle = Bio.Entrez.efetch(db="protein", rettype="ipg", retmode="xml", id=self.accession_number)
                xml_output = handle.read().decode('utf-8')
            root = xml.etree.cElementTree.fromstring(xml_output)
            list_of_kingdom_taxid = []
            assemblies_coordinates = []
            for report in root.iter("IPGReport"):
                product = report.find("Product")
                if "product_acc" in report.attrib.keys():
                    report_accession_number = report.attrib["product_acc"]
                elif "accver" in product.attrib.keys():
                    report_accession_number = product.attrib["accver"]
                else:
                    report_accession_number = ""
                if report_accession_number == self.accession_number:  # be careful
                    for protein in report.iter("Protein"):
                        if protein.attrib["source"] == "RefSeq":
                            if "name" in protein.attrib.keys():
                                self.name = protein.attrib["name"]
                            self.taxid = protein.attrib["taxid"]
                            self.kingdom_taxid = protein.attrib["kingdom_taxid"]
                            self.organism = protein.attrib["org"]
                            list_of_kingdom_taxid.append(self.kingdom_taxid)
                            for cds in protein.iter("CDS"):
                                to_add = 1
                                if self.parameters.arguments["filter_refseq_sequences_by_regex"]:
                                    if not re.search(rf"{self.parameters.arguments['refseq_sequences_regex']}",
                                                     cds.attrib["accver"]):
                                        to_add = 0
                                if to_add == 1:
                                    if "assembly" not in cds.attrib.keys():
                                        cds.attrib["assembly"] = "NA"
                                    if "strain" not in cds.attrib.keys():
                                        cds.attrib["strain"] = "NA"
                                    try:
                                        assemblies_coordinates.append(dict(locus_id=cds.attrib["accver"],
                                                                           start=(int(cds.attrib["start"]) - 1),
                                                                           stop=int(cds.attrib["stop"]),
                                                                           strand=cds.attrib['strand'],
                                                                           length=int(cds.attrib["stop"]) - (
                                                                                   int(cds.attrib["start"]) - 1),
                                                                           assembly=cds.attrib["assembly"],
                                                                           strain=cds.attrib["strain"],
                                                                           org=cds.attrib["org"],
                                                                           taxid=cds.attrib["taxid"]))
                                    except:
                                        print(f"❕Attention: {cds.attrib['accver']} record is not completed and"
                                              f" cannot be processed", file=sys.stderr)
            '''
            if len(assemblies_coordinates) == 0:
                print(f"❗Warning message:\n\tNo assembly was found for the protein "
                      f"'{self.accession_number}'.\n\tThis protein record can be suppressed by the ncbi\n\t"
                      f"or it has no sequence record that satisfies refseq_sequnces_regex config parameter.",
                      file=sys.stderr)
            '''
            self.assemblies_coordinates = assemblies_coordinates
            return assemblies_coordinates
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to get assemblies coordinates of a protein.") from error

    '''
    def get_loci(self, start=-float("inf"), end=float("inf"), strand="NA") -> dict:
        """Get Locus class objects for each sequence from the ncbi nt database on which the protein is annotated.

        Returns:
            dict: Dict with keys as locus_ids and values as Locus class' objects.

        """
        self.loci = dict()
        for assembly in self.assemblies_coordinates:
            locus_id = assembly["locus_id"]
            self.loci[locus_id] = Locus(locus_id, start_b=start, end_b=end, strand=strand)
        return self.loci
    '''

    def blastp_searching_for_homologues(self) -> list:
        """Search for a protein's homologues with blastp against the 'refseq_protein' database.

        Note:
            This function does not create a new object's attribute; It only returns a list of accession numbers.

        Returns:
            list: List of proteins' accession numbers obtained with blastp searching. This list also contains the query
                protein's accession number.

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(
                    f"👀 Searching for homologues of {self.accession_number} with blastp against the RefSeq database...",
                    file=sys.stdout)
            handle = Bio.Blast.NCBIWWW.qblast("blastp", "refseq_protein", self.accession_number,
                                              expect=self.parameters.arguments["blastp_evalue_cutoff"],
                                              hitlist_size=self.parameters.arguments["blastp_hit_list_size"],
                                              alignments=self.parameters.arguments["blastp_max_number_of_alignments"])
            xml_output = handle.read()
            hits_an_list = [self.accession_number]
            blastp_stat_dict = dict()
            blastp_stat_dict[self.accession_number] = dict(pident_to_query_length="the query",
                                                           pident_to_sequence_length="the query",
                                                           pident_to_alignment_length="the query", evalue="the query")
            root = xml.etree.ElementTree.fromstring(xml_output)
            query_length = int(root.find("BlastOutput_query-len").text)
            for hit in root.iter("Hit"):
                hit_id = hit.find("Hit_id").text.strip("ref").strip("|")
                if hit_id != self.accession_number:
                    hit_description = hit.find("Hit_def").text
                    subject_length = int(hit.find("Hit_len").text)
                    hsp_identity_sum, hsp_positive_sum, hsp_align_length = 0, 0, 0
                    evalue = []
                    for hsp in hit.iter("Hsp"):
                        hsp_identity_sum += int(hsp.find("Hsp_identity").text)
                        hsp_positive_sum += int(hsp.find("Hsp_positive").text)
                        hsp_align_length += int(hsp.find("Hsp_align-len").text)
                        evalue.append(hsp.find("Hsp_evalue").text)
                    pident_to_query_length = hsp_identity_sum / query_length
                    pident_to_seq_length = hsp_identity_sum / subject_length
                    pident_to_alignment_length = hsp_identity_sum / hsp_align_length
                    if pident_to_query_length >= self.parameters.arguments["blastp_pident_to_query_length_cutoff"]:
                        blastp_stat_dict[hit_id] = dict(pident_to_query_length=str(round(pident_to_query_length, 4)),
                                                        pident_to_sequence_length=str(round(pident_to_seq_length, 4)),
                                                        pident_to_alignment_length=str(
                                                            round(pident_to_alignment_length, 4)),
                                                        evalue=",".join(evalue))
                        if hit_id not in hits_an_list:
                            hits_an_list.append(hit_id)
            columns = "\t".join(["accession_number", "name", "pident_to_query_length", "pident_to_sequence_length",
                                 "pident_to_alignment_length", "e-value"])
            table = [columns]
            hits_records_list = [RefSeqProtein(i, self.parameters) for i in hits_an_list]
            for i in range(0, len(hits_records_list), 200):
                records_subset = hits_records_list[i:i + 200]
                accession_numbers = [record.accession_number for record in records_subset]
                handle_fasta = Bio.Entrez.efetch(db="protein", id=accession_numbers, rettype="fasta", retmode="text")
                fasta_records = Bio.SeqIO.parse(handle_fasta, "fasta")
                for f_record in fasta_records:
                    record_index = accession_numbers.index(f_record.id)
                    records_subset[record_index].name = f_record.description.replace(f_record.id, "").strip()
            for rec in hits_records_list:
                table.append("\t".join([rec.accession_number, rec.name,
                                        blastp_stat_dict[rec.accession_number]["pident_to_query_length"],
                                        blastp_stat_dict[rec.accession_number]["pident_to_sequence_length"],
                                        blastp_stat_dict[rec.accession_number]["pident_to_alignment_length"],
                                        blastp_stat_dict[rec.accession_number]["evalue"]]))
            if not os.path.exists(self.parameters.arguments["output_dir"]):
                os.mkdir(self.parameters.arguments["output_dir"])
            output_filename = os.path.join(self.parameters.arguments["output_dir"], "found_homologues.tsv")
            f = open(output_filename, "w")
            f.write("\n".join(table))
            if self.parameters.arguments["verbose"]:
                print(f"✅ {len(hits_records_list) - 1} homologues were found.\n"
                      f"💌 Summary table was saved to: {os.path.basename(output_filename)}", file=sys.stdout)
            return hits_an_list
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to perform searching for homologues with blastp.") from error


class Locus:
    """
    A Locus object holds sequence and annotation of the corresponding ncbi Reference Sequence.

    Attributes:
        locus_id (str): a NCBI locus id from the Nucleotide database.
        locus_record (Bio.SeqRecord.SeqRecord): a biopython record object of the sequence.
        CDSs (list): list of dicts with information about annotated CDS in the locus' sequence.
        start_b (int): start of region within annotation should be retrieved.
        stop_b (int): stop of region within annotation should be retrieved.

    """

    def __init__(self, locus_id: str, start_b: int = 0, stop_b: int = None, target_strand: str = "NA",
                 locus_record=None, xml_output=None):
        """Create a Locus object.

        Note:
            0-based format is used for sequence indexing.

        Arguments:
            locus_id (str): locus id from the ncbi nucleotide database.
            start_b (int): start of region within annotation should be retrieved (optional).
            stop_b (int): stop of region within annotation should be retrieved (optional).
            target_strand (str): strand of the target object (optional).

        """
        try:
            self.locus_id = locus_id
            if not locus_record:
                handle = Bio.Entrez.efetch(db="nucleotide", rettype="fasta", retmode="txt", id=locus_id)
                self.locus_record = Bio.SeqIO.read(handle, "fasta")
            else:
                self.locus_record = locus_record
            if stop_b is None:
                stop_b = len(self.locus_record.seq)
            if not xml_output:
                handle = Bio.Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="xml", id=locus_id)
                xml_output = (handle.read()).decode("utf-8")
            root = xml.etree.ElementTree.fromstring(xml_output)
            self.CDSs = []
            for gbseq in root.iter("GBSeq"):
                if gbseq.find("GBSeq_accession-version").text == self.locus_id:
                    for gbfeature in gbseq.iter("GBFeature"):
                        if gbfeature.find("GBFeature_key").text == "CDS":
                            try:
                                starts, stops = [], []
                                for interval in gbfeature.iter("GBInterval"):
                                    try:
                                        start, stop = int(interval.find("GBInterval_from").text), int(
                                            interval.find("GBInterval_to").text)
                                        if start > stop:
                                            start, stop, strand = stop - 1, start, "-"
                                        else:
                                            start, stop, strand = start - 1, stop, "+"
                                        starts.append(start)
                                        stops.append(stop)
                                    except:
                                        pass
                                if starts:
                                    coordinates = list(sorted(zip(starts, stops), key=lambda pair: pair[0]))
                                    main_start, main_stop = coordinates[0][0], coordinates[-1][-1]
                                    if strand == "+":
                                        main_stop = main_stop - 3
                                    elif strand == "-":
                                        main_start = main_start + 3
                                    relative_start, relative_stop = main_start - start_b, main_stop - start_b
                                    if strand == target_strand:
                                        relative_strand = "+"
                                    else:
                                        relative_strand = "-"
                                        useq_length = stop_b - start_b
                                    if target_strand == "-":
                                        relative_start, relative_stop = useq_length - relative_stop, useq_length - relative_start

                                    if (start_b <= main_start < stop_b) or (start_b <= main_stop < stop_b):
                                        cds_seq = self.locus_record.seq[main_start:main_stop]
                                        if strand == '-':
                                            cds_seq = cds_seq.reverse_complement()
                                        protein_id, product_name = 'NA', 'NA'
                                        for gbqualifier in gbfeature.iter("GBQualifier"):
                                            if gbqualifier.find("GBQualifier_name").text == "protein_id":
                                                protein_id = gbqualifier.find("GBQualifier_value").text
                                            if gbqualifier.find("GBQualifier_name").text == "product":
                                                product_name = gbqualifier.find("GBQualifier_value").text
                                        if protein_id != 'NA':
                                            if product_name != 'NA':
                                                product_name = f"{protein_id} ({product_name})"
                                            else:
                                                product_name = f"{protein_id}"
                                            self.CDSs.append(dict(protein_id=protein_id, product_name=product_name,
                                                                  coordinates=coordinates, nt_seq=cds_seq,
                                                                  main_start=main_start, main_stop=main_stop,
                                                                  strand=strand,
                                                                  relative_start=relative_start,
                                                                  relative_stop=relative_stop,
                                                                  relative_strand=relative_strand))
                            except:
                                pass
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to create a Locus class' object.") from error


class Homologues:
    """A Homologues object holds list of proteins homologues and information about them.

    Attributes:
        accession_numbers (list): List of RefSeq accession numbers.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        records (list): list of RefSeqProtein objects of the proteins.

    """

    def __init__(self, accession_numbers: list, parameters: uorf4u.manager.Parameters):
        """Create a Homologues object.

        Note:
            With initialisation it also creates a 'records' attribute - a list of RefSeqProtein objects of proteins
                based on accession numbers list.

        Arguments:
            accession_numbers (list): List of RefSeq accession numbers.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        try:
            self.accession_numbers = accession_numbers
            self.parameters = parameters
            self.records = [RefSeqProtein(i, parameters) for i in accession_numbers]
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to create a Homologues class' object.") from error

    def get_upstream_sequences(self) -> list:
        """Get upstream sequences of proteins' genes.

        Note:
            A protein may be found in multiple assemblies (for example in different strains).

        Returns:
            list: List of Bio.SeqRecord.SeqRecord objects of upstream sequences.

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"📡 Retrieving upstream sequences...", file=sys.stdout)
            for i in range(0, len(self.records), 200):
                records_subset = self.records[i:i + 200]
                accession_numbers = [record.accession_number for record in records_subset]
                handle = Bio.Entrez.efetch(db="protein", id=accession_numbers, rettype="ipg", retmode="xml")
                handle_txt = handle.read().decode('utf-8')
                for record in records_subset:
                    record.get_assemblies(handle_txt)
                handle_fasta = Bio.Entrez.efetch(db="protein", id=accession_numbers, rettype="fasta", retmode="text")
                fasta_records = Bio.SeqIO.parse(handle_fasta, "fasta")
                for f_record in fasta_records:
                    record_index = accession_numbers.index(f_record.id)
                    records_subset[record_index].record = f_record

            proteins_wo_assemblies = []
            if self.parameters.arguments["assemblies_list"] == 'NA':
                assemblies_table = [f"accession_number\tlocus_id\tassembly\torganism\tstrain\ttax_id"]
                list_of_protein_with_multiple_assemblies = []
                numbers_of_assemblies = []
                for record in self.records:
                    numbers_of_assemblies.append(len(record.assemblies_coordinates))
                    if len(record.assemblies_coordinates) == 0:
                        proteins_wo_assemblies.append(record.accession_number)
                    if len(record.assemblies_coordinates) > 1:
                        list_of_protein_with_multiple_assemblies.append(record.accession_number)
                    for assembly in record.assemblies_coordinates:
                        assemblies_table.append(
                            f"{record.accession_number}\t"
                            f"{assembly['locus_id']}:{assembly['start']}:{assembly['stop']}({assembly['strand']})"
                            f"\t{assembly['assembly']}"
                            f"\t{assembly['org']}\t{assembly['strain']}\t{assembly['taxid']}")
                if not os.path.exists(self.parameters.arguments["output_dir"]):
                    os.mkdir(self.parameters.arguments["output_dir"])
                assemblies_table_path = os.path.join(self.parameters.arguments["output_dir"], "assemblies_list.tsv")
                assemblies_selected_table_path = os.path.join(self.parameters.arguments["output_dir"],
                                                              "selected_assemblies_list.tsv")
                assemblies_table_file = open(assemblies_table_path, "w")
                assemblies_table_file.write("\n".join(assemblies_table))
                assemblies_table_file.close()

                proteins_wo_assemblies_txt = "\n".join(proteins_wo_assemblies) + "\n"
                proteins_wo_assemblies_path = os.path.join(self.parameters.arguments["output_dir"],
                                                           "proteins_wo_assembly.txt")
                proteins_wo_assemblies_file = open(proteins_wo_assemblies_path, "w")
                proteins_wo_assemblies_file.write(proteins_wo_assemblies_txt)

                if numbers_of_assemblies.count(0) > 0:
                    print(f"❗️Warning message:\n\tFor {numbers_of_assemblies.count(0)} proteins "
                          f"no assembly was found.\n"
                          f"\tThese proteins' records can be suppressed by the ncbi\n\t"
                          f"or they don't have loci that satisfies refseq_sequnces_regex config parameter.\n\t"
                          f"List of these proteins was saved as: {os.path.basename(proteins_wo_assemblies_path)}",
                          file=sys.stderr)
                if len(list_of_protein_with_multiple_assemblies) > 0:
                    print(f"❗️Warning message:\n\tFor {len(list_of_protein_with_multiple_assemblies)} proteins "
                          f"multiple assemblies were found in identical protein database\n"
                          f"\twith max number of assemblies per one protein as {max(numbers_of_assemblies)} 😱.\n\t"
                          f"A table with information about the assemblies was saved as a tsv file: "
                          f"{os.path.basename(assemblies_table_path)}.\n\tYou can edit it and remove lines with assemblies "
                          f"you do not want to include in your analysis.\n"
                          f"\tAfter filtering, you can use -al cmd parameter with your table as an argument.\n"
                          f"\tIn addition, config file has 'max_number_of_assemblies' parameter "
                          f"(set as {self.parameters.arguments['max_number_of_assemblies']}).\n\tBy default ❕, it's used "
                          f"by uorf4u to limit max number of assemblies included in the analysis;\n"
                          f"\tand it works only if '-al' option is not provided. In case number of assemblies is more than "
                          f"the cutoff,\n\trandom sampling 🎲 will be used to take only subset of them.\n\t"
                          f"Selected assemblies information was savead as a tsv file: "
                          f"{os.path.basename(assemblies_selected_table_path)}"
                          f"\n\tSee documentation 📖 for details.", file=sys.stderr)
            else:
                assemblies_table = pandas.read_table(self.parameters.arguments["assemblies_list"], sep="\t")
                locus_ids = assemblies_table["locus_id"].to_list()
                locus_ids = [id.split(":")[0] for id in locus_ids]

            upstream_sequences = []
            an_with_no_annotated_useq = []
            for record in self.records:
                assemblies = record.assemblies_coordinates
                if isinstance(self.parameters.arguments["max_number_of_assemblies"], int) and \
                        self.parameters.arguments["assemblies_list"] == "NA":
                    if len(assemblies) >= self.parameters.arguments["max_number_of_assemblies"]:
                        assemblies = random.sample(assemblies, self.parameters.arguments["max_number_of_assemblies"])
                if self.parameters.arguments["assemblies_list"] != "NA":
                    assemblies_filtered = [i for i in assemblies if i["locus_id"] in locus_ids]
                    assemblies = assemblies_filtered
                record.assemblies_coordinates = assemblies

            assemblies_table = [f"accession_number\tlocus_id\tassembly\torganism\tstrain\ttax_id"]
            for record in self.records:
                for assembly in record.assemblies_coordinates:
                    assemblies_table.append(
                        f"{record.accession_number}\t"
                        f"{assembly['locus_id']}:{assembly['start']}:{assembly['stop']}({assembly['strand']})"
                        f"\t{assembly['assembly']}"
                        f"\t{assembly['org']}\t{assembly['strain']}\t{assembly['taxid']}")
            assemblies_table_file = open(assemblies_selected_table_path, "w")
            assemblies_table_file.write("\n".join(assemblies_table))
            assemblies_table_file.close()

            lists_of_assemblies = [record.assemblies_coordinates for record in self.records]
            all_assemblies = [assembly for sublist in lists_of_assemblies for assembly in sublist]
            for i in range(0, len(all_assemblies), 150):
                assemblies_subset = all_assemblies[i:i + 150]
                sequences_ids = [assembly["locus_id"] for assembly in assemblies_subset]
                handle = Bio.Entrez.efetch(db="nucleotide", rettype="fasta", retmode="txt", id=sequences_ids)
                records = Bio.SeqIO.parse(handle, "fasta")
                for record, assembly in zip(records, assemblies_subset):
                    assembly["record"] = record

            for record in self.records:
                record_upstream_sequences = []
                for assembly in record.assemblies_coordinates:
                    locus_record = assembly["record"]
                    try:
                        useq_downstream_region_length = min(self.parameters.arguments["downstream_region_length"],
                                                            len(record.record.seq) * 3)
                    except:
                        useq_downstream_region_length = self.parameters.arguments["downstream_region_length"]
                    useq_upstream_region_length = self.parameters.arguments["upstream_region_length"]
                    if assembly["strand"] == "+":
                        if self.parameters.arguments["upstream_region_length"] == "all":
                            useq_start = 0
                        else:
                            useq_start = max(0, assembly["start"] - self.parameters.arguments["upstream_region_length"])
                        if useq_start == 0:
                            useq_upstream_region_length = assembly["start"]
                        useq_stop = min(assembly["start"] + self.parameters.arguments["downstream_region_length"],
                                        len(locus_record.seq))
                        if useq_stop == len(locus_record.seq):
                            useq_downstream_region_length = len(locus_record.seq) - assembly["start"]
                    elif assembly["strand"] == "-":
                        useq_start = max(0, assembly["stop"] - self.parameters.arguments["downstream_region_length"])
                        if useq_start == 0:
                            useq_downstream_region_length = assembly["stop"]
                        if self.parameters.arguments["upstream_region_length"] == "all":
                            useq_stop = len(locus_record.seq)
                        else:
                            useq_stop = min(len(locus_record.seq),
                                            assembly["stop"] + self.parameters.arguments["upstream_region_length"])
                        if useq_stop == len(locus_record.seq):
                            useq_upstream_region_length = len(locus_record.seq) - assembly["stop"]
                    useq_length = abs(useq_stop - useq_start)
                    if self.parameters.arguments["upstream_region_length"] != 'all':
                        if self.parameters.arguments["minimal_upstream_region_length"] >= self.parameters.arguments[
                            "upstream_region_length"]:
                            self.parameters.arguments["minimal_upstream_region_length"] = self.parameters.arguments[
                                "upstream_region_length"]
                    if useq_upstream_region_length >= self.parameters.arguments["minimal_upstream_region_length"] or \
                            self.parameters.arguments["upstream_region_length"] == "all":
                        useq = locus_record.seq[useq_start:useq_stop]
                        if assembly["strand"] == "-":
                            useq = useq.reverse_complement()
                        if assembly["strain"] == "NA":
                            useq_name = assembly["org"]
                        elif assembly["strain"] in assembly["org"]:
                            useq_name = f"{assembly['org'].replace(assembly['strain'], '')}{assembly['strain']}"
                        else:
                            useq_name = f"{assembly['org']} {assembly['strain']}"
                        useq_id = f"{assembly['locus_id']}|{useq_start}-{useq_stop}({assembly['strand']})|" \
                                  f"{record.accession_number}"
                        # useq_id = f"{useq_name}|{assembly['locus_id']}|{record.accession_number}"
                        useq_label = f"{useq_name}|{assembly['locus_id']}|{record.accession_number}"
                        useq_annotations = dict(RefSeq=True, locus_record=locus_record,
                                                locus_id=assembly['locus_id'], length=useq_length,
                                                start=useq_start, stop=useq_stop, strand=assembly["strand"],
                                                accession_number=record.accession_number,
                                                organism=assembly['org'], label=useq_label,
                                                upstream_region_length=useq_upstream_region_length,
                                                downstream_region_length=useq_downstream_region_length)
                        useq_record = Bio.SeqRecord.SeqRecord(useq, id=useq_id, description=useq_name,
                                                              annotations=useq_annotations)
                        record_upstream_sequences.append(useq_record)
                upstream_sequences += record_upstream_sequences
                if len(record_upstream_sequences) == 0:
                    an_with_no_annotated_useq.append(record.accession_number)
            if an_with_no_annotated_useq:
                print(f"❗Warning message:\n\tNo upstream sequences for {len(an_with_no_annotated_useq)} protein(s)"
                      f" were annotated.\n\tCorresponding loci in the nucleotide ncbi database can be too short 📏.\n"
                      f"\tSee 'minimal_upstream_region_length' config parameter description in the documentation.",
                      file=sys.stderr)
            if self.parameters.arguments["verbose"]:
                print(f"✅ {len(upstream_sequences)} upstream sequences were obtained.", file=sys.stdout)
            return upstream_sequences
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to retrieve upstream sequences.") from error


class UpstreamSequences:
    """An UpstreamSequences object holds list of upstream sequences records and information about them.

    Attributes:
        records (list): List of Bio.SeqRecord.SeqRecord objects with upstream sequences.
            (attribute 'annotations' (dict) is used for holding additional information, (e.g. downstream protein_id)).
        codon_table (Bio.Data.CodonTable.CodonTable): Codon table (genetic code).
        conserved_paths (list): list  of Path's objects (Path class holds list of ORFs from different upstream
            sequences and information about them).

    """

    def __init__(self, records: list, parameters: uorf4u.manager.Parameters):
        """Create an UpstreamSequences object.

        Arguments:
            records (list): List of Bio.SeqRecord.SeqRecord objects with upstream sequences.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        self.records = records
        self.parameters = parameters
        self.codon_table = Bio.Data.CodonTable.unambiguous_dna_by_name[
            parameters.arguments["ncbi_genetic_code_name"]]
        self.conserved_paths = None

    def save_upstream_sequences(self) -> None:
        """Save upstream sequences as a fasta file.

        Returns:
            None

        """
        try:
            output_file = os.path.join(self.parameters.arguments["output_dir"], "upstream_sequences.fa")
            if not os.path.exists(self.parameters.arguments["output_dir"]):
                os.mkdir(self.parameters.arguments["output_dir"])
            Bio.SeqIO.write(self.records, output_file, "fasta")
            if self.parameters.arguments["verbose"]:
                print(f"💌 Fasta file with upstream sequences was saved to {os.path.basename(output_file)}.",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to save a fasta file with upstream sequences.") from error

    def annotate_orfs(self) -> None:
        """Annotate ORFs in upstream sequences.

        Note:
            This function updates 'records' attribute.

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🔎 Annotating ORFs in the upstream sequences...", file=sys.stdout)
            if self.parameters.arguments["alternative_start_codons"]:
                start_codons_list = self.codon_table.start_codons
            else:
                start_codons_list = [self.parameters.arguments["main_start_codon"]]

            if self.parameters.arguments["check_assembly_annotation"] and \
                    self.records[0].annotations["RefSeq"]:
                if self.parameters.arguments["verbose"]:
                    print(f"📡 Retrieving assemblies' annotation...", file=sys.stdout)
                for i in range(0, len(self.records), 100):
                    useq_subset = [record for record in self.records[i:i + 100] if record.annotations["RefSeq"]]
                    locus_ids = [locus.annotations["locus_id"] for locus in useq_subset]
                    handle = Bio.Entrez.efetch(db="nucleotide", id=locus_ids, rettype="gb", retmode="xml")
                    handle_txt = handle.read().decode('utf-8')
                    for useq_record in useq_subset:
                        useq_record.annotations["locus_annotation"] = Locus(useq_record.annotations["locus_id"],
                                                                            start_b=useq_record.annotations["start"],
                                                                            stop_b=useq_record.annotations["stop"],
                                                                            target_strand=useq_record.annotations[
                                                                                "strand"],
                                                                            locus_record=useq_record.annotations[
                                                                                "locus_record"],
                                                                            xml_output=handle_txt)
            for useq_index in range(len(self.records)):
                useq_record = self.records[useq_index]
                useq_record.annotations["ORFs"] = []
                for first_position in range((useq_record.annotations["length"] - self.parameters.arguments[
                    "downstream_region_length"]) + 1):
                    first_codon = useq_record.seq[first_position:first_position + 3]
                    if first_codon.upper() in start_codons_list:
                        start_codon_position = first_position
                        for second_position in range(start_codon_position + 3,
                                                     (useq_record.annotations["length"] - 3) + 1, 3):
                            second_codon = useq_record.seq[second_position:second_position + 3]
                            if second_codon.upper() in self.codon_table.stop_codons:
                                stop_codon_position = second_position
                                orf_length = stop_codon_position - start_codon_position
                                distance = (useq_record.annotations["length"] - self.parameters.arguments[
                                    "downstream_region_length"]) - stop_codon_position
                                distance_sc = (useq_record.annotations["length"] - self.parameters.arguments[
                                    "downstream_region_length"]) - start_codon_position

                                if useq_record.annotations["RefSeq"]:
                                    orf_id = f"{useq_record.annotations['locus_id']}|" \
                                             f"{useq_record.annotations['accession_number']}|" \
                                             f"{distance}"
                                    orf_name = f"{useq_record.annotations['label']}|{distance}"
                                else:
                                    distance = useq_record.annotations["length"] - stop_codon_position
                                    orf_id = f"{useq_record.id}|{distance}"
                                    if useq_record.description:
                                        orf_name = f"{useq_record.description}|{orf_id}"
                                    else:
                                        orf_name = orf_id
                                sd_window_start = max(
                                    [0, (start_codon_position - self.parameters.arguments["sd_window_length"])])
                                current_orf = ORF(parameters=self.parameters, id=orf_id, name=orf_name,
                                                  distance=distance, start=start_codon_position,
                                                  stop=stop_codon_position, useq_index=useq_index,
                                                  nt_sequence=useq_record.seq[start_codon_position:stop_codon_position],
                                                  sd_window_seq=useq_record.seq[sd_window_start:start_codon_position])
                                if current_orf.length >= self.parameters.arguments[
                                    "min_orf_length"] and distance_sc != 0:
                                    useq_record.annotations["ORFs"].append(current_orf)
                                    if self.parameters.arguments["check_assembly_annotation"] and \
                                            useq_record.annotations["RefSeq"]:
                                        for cds in useq_record.annotations["locus_annotation"].CDSs:
                                            if current_orf.stop == cds["relative_stop"] and (
                                                    (current_orf.start - cds["relative_start"]) % 3 == 0):
                                                the_same_stop = 1
                                                current_orf.annotation = cds["product_name"]
                                                if current_orf.start != cds["relative_start"]:
                                                    if current_orf.start < cds["relative_start"]:
                                                        current_orf.annotation += " (extension)"
                                                    else:
                                                        current_orf.annotation += " (truncation)"
                                    for annotated_orfs in useq_record.annotations["ORFs"]:
                                        if current_orf.stop == annotated_orfs.stop and \
                                                current_orf.id != annotated_orfs.id:
                                            current_orf.extended_orfs.append(annotated_orfs.id)
                                break
            number_of_orfs = sum(len(i.annotations["ORFs"]) for i in self.records)
            if self.parameters.arguments["fast_searching"] == "auto":
                if len(self.records) < 5:
                    self.parameters.arguments["fast_searching"] = False
                elif (len(self.records) >= 100 or number_of_orfs > 1000):
                    self.parameters.arguments["fast_searching"] = True
                else:
                    self.parameters.arguments["fast_searching"] = False
            if number_of_orfs == 0:
                print(f"⛔Termination:\n\tNo ORF was annotated in upstream sequences."
                      f"\n\tThis run will be terminated.", file=sys.stderr)
                if not os.path.exists(self.parameters.arguments["output_dir"]):
                    os.mkdir(self.parameters.arguments["output_dir"])
                with open(os.path.join(self.parameters.arguments["output_dir"], "report.txt"), "w") as report_f:
                    report_f.write("Termination:\nNo ORF was annotated in upstream sequences.")
                sys.exit()
            if self.parameters.arguments["verbose"]:
                print(f"✅ {number_of_orfs} ORFs were annotated.", file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to annotate ORFs in upstream sequences.") from error

    def filter_orfs_by_sd_annotation(self) -> None:
        """Filter annotated ORFs by presence the Shine-Dalgarno sequence.

        Returns:
            None

        """
        try:
            for useq_record in self.records:
                orf_list = useq_record.annotations["ORFs"]
                filtered_orf_list = []
                for orf in orf_list:
                    orf.calculate_energies()
                    if orf.min_energy < self.parameters.arguments["sd_energy_cutoff"]:
                        filtered_orf_list.append(orf)
                useq_record.annotations["ORFs"] = filtered_orf_list

            number_of_orfs = sum(len(i.annotations["ORFs"]) for i in self.records)
            if number_of_orfs == 0:
                print(f"⛔Termination:\n\tNo ORF left after filtering by SD annotation."
                      f"\n\tThis run will be terminated.", file=sys.stderr)
                if not os.path.exists(self.parameters.arguments["output_dir"]):
                    os.mkdir(self.parameters.arguments["output_dir"])
                with open(os.path.join(self.parameters.arguments["output_dir"], "report.txt"), "w") as report_f:
                    report_f.write("Termination:\nNo ORF left after filtering by SD annotation.")
                sys.exit()
            if self.parameters.arguments["verbose"]:
                print(f"🧹 {number_of_orfs} ORFs remained in the analysis after filtering by presence "
                      f"of the SD sequence.", file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to filter uORFs by SD sequence presence.") from error

    def save_annotated_orfs(self) -> None:
        """Save information about annotated ORFs as a set of tsv files.

        Note:
            tsv files will be saved to the subdir called 'annotated_ORFs' located in 'output_dir'.

        Returns:
            None

        """
        try:
            colnames = "\t".join(
                ["id", "name", "length", "nt_sequence", "aa_sequence", "sd_sequence_window", "SD-aSD energy",
                 "SD-aSD energies list", "extended_orfs", "annotation"])
            if not os.path.exists(self.parameters.arguments["output_dir"]):
                os.mkdir(self.parameters.arguments["output_dir"])
            output_dir_path = os.path.join(self.parameters.arguments["output_dir"], "annotated_ORFs")
            if not os.path.exists(output_dir_path):
                os.mkdir(output_dir_path)
            for useq_record in self.records:
                file_name = f"{useq_record.description}|{useq_record.id}".replace(' ', '_').replace('/', '_')
                lines = [colnames]
                for orf in useq_record.annotations["ORFs"]:
                    if not orf.extended_orfs:
                        extented_orfs_value = "NA"
                    else:
                        extented_orfs_value = ';'.join(orf.extended_orfs)
                    lines.append("\t".join(
                        [orf.id, orf.name, str(orf.length), str(orf.nt_sequence), str(orf.aa_sequence),
                         str(orf.sd_window_seq_str), str(orf.min_energy),";".join(orf.sd_window_energies),
                         extented_orfs_value, orf.annotation]))
                with open(os.path.join(output_dir_path, f"{file_name}.tsv"), "w") as output:
                    output.write("\n".join(lines))
            if self.parameters.arguments["verbose"]:
                print(f"💌 tsv files with information about annotated ORFs were saved to "
                      f"{os.path.basename(output_dir_path)} folder.",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to save annotated uORFs.") from error

    def conserved_orf_searching(self) -> None:
        """Search for sets of conserved ORFs in upstream sequences.

        Note:
            This method updates the self.conserved_paths attribute.

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🔎 Searching for conserved ORFs in upstream sequences...",
                      file=sys.stdout)
            if len(self.records) == 1:
                raise uorf4u.manager.uORF4uError("At least two sequences required to perform conservation analysis")
            lengths = []
            for record in self.records:
                for orf in record.annotations["ORFs"]:
                    lengths.append(orf.length)
            lengths = sorted(list(set(lengths)))
            global_aligner = Bio.Align.PairwiseAligner()
            global_aligner.mode = "global"
            global_aligner.match_score = self.parameters.arguments["global_match_score"]
            global_aligner.mismatch_score = self.parameters.arguments["global_mismatch_score"]
            global_aligner.open_gap_score = self.parameters.arguments["global_open_gap_score"]
            global_aligner.extend_gap_score = self.parameters.arguments["global_extend_gap_score"]
            global_aligner.target_end_gap_score = self.parameters.arguments["global_target_end_gap_score"]
            global_aligner.query_end_gap_score = self.parameters.arguments["global_query_end_gap_score"]
            length_variance = self.parameters.arguments["orf_length_group_range"]
            number_of_useqs = len(self.records)
            if self.parameters.arguments["fast_searching"]:
                filtered_orfs_dict = dict()
                for length in lengths:
                    if isinstance(self.parameters.arguments["orf_length_group_range"], float):
                        length_variance = length * self.parameters.arguments["orf_length_group_range"]
                    filtered_orfs = []
                    useq_with_filtered_orfs = []
                    for useq_index in range(number_of_useqs):
                        useq_record = self.records[useq_index]
                        for orf in useq_record.annotations["ORFs"]:
                            if abs(length - orf.length) <= length_variance:
                                filtered_orfs.append(orf)
                                if useq_index not in useq_with_filtered_orfs:
                                    useq_with_filtered_orfs.append(useq_index)
                    if len(useq_with_filtered_orfs) / number_of_useqs >= self.parameters.arguments[
                        "orfs_presence_cutoff"]:
                        to_add = 1
                        keys_to_remove = []
                        for added_length in filtered_orfs_dict.keys():
                            num_of_identical_elements = len(set(filtered_orfs) & set(filtered_orfs_dict[added_length]))
                            fraction = num_of_identical_elements / min(len(filtered_orfs),
                                                                       len(filtered_orfs_dict[added_length]))
                            if fraction > 0.8:  # to add as a config parameter
                                if len(filtered_orfs) >= len(filtered_orfs_dict[added_length]):
                                    keys_to_remove.append(added_length)
                                else:
                                    to_add = 0
                        for key_to_remove in keys_to_remove:
                            filtered_orfs_dict.pop(key_to_remove)
                        if to_add:
                            filtered_orfs_dict[length] = filtered_orfs
                lengths = list(filtered_orfs_dict.keys())

            conserved_paths = []
            for length in lengths:
                if isinstance(self.parameters.arguments["orf_length_group_range"], float):
                    length_variance = length * self.parameters.arguments["orf_length_group_range"]
                useq_indexes_with_filtered_orfs = []
                filtered_orfs = dict()
                for useq_index in range(number_of_useqs):
                    useq_record = self.records[useq_index]
                    filtered_orfs[useq_index] = []
                    for orf in useq_record.annotations["ORFs"]:
                        if abs(length - orf.length) <= length_variance:
                            filtered_orfs[useq_index].append(orf)
                    orfs_ids = [i.id for i in filtered_orfs[useq_index]]
                    for orf in filtered_orfs[useq_index]:
                        if any(i in orf.extended_orfs for i in orfs_ids):
                            filtered_orfs[useq_index].remove(orf)
                    if len(filtered_orfs[useq_index]) > 0:
                        useq_indexes_with_filtered_orfs.append(useq_index)
                if len(useq_indexes_with_filtered_orfs) / number_of_useqs >= self.parameters.arguments[
                    "orfs_presence_cutoff"]:
                    if self.parameters.arguments["fast_searching"]:
                        genome_iterator = random.sample(filtered_orfs.keys(),
                                                        max(1, min(round(self.parameters.arguments["fast_searching_"
                                                                                                   "fraction_of_initial"
                                                                                                   "_genomes"] * len(
                                                            useq_indexes_with_filtered_orfs)),
                                                                   self.parameters.arguments[
                                                                       "max_num_of_initial_genome_iteration"])))
                    elif len(filtered_orfs.keys()) > self.parameters.arguments["max_num_of_initial_genome_iteration"]:
                        genome_iterator = random.sample(filtered_orfs.keys(),
                                                        self.parameters.arguments[
                                                            "max_num_of_initial_genome_iteration"])
                    else:
                        genome_iterator = filtered_orfs.keys()
                    for initial_useq in genome_iterator:
                        for initial_orf in filtered_orfs[initial_useq]:
                            conserved_path = Path(self.parameters)
                            conserved_path.update(initial_orf)
                            for useq in random.sample(filtered_orfs.keys(), len(filtered_orfs.keys())):
                                if useq != initial_useq and filtered_orfs[useq] != []:
                                    score_sums = []
                                    for orf in filtered_orfs[useq]:
                                        score_sum = 0
                                        for path_orf in conserved_path.path:
                                            if self.parameters.arguments["alignment_type"] == "nt":
                                                current_alignment = global_aligner.align(orf.nt_sequence,
                                                                                         path_orf.nt_sequence)
                                            elif self.parameters.arguments["alignment_type"] == "aa":
                                                current_alignment = global_aligner.align(orf.aa_sequence,
                                                                                         path_orf.aa_sequence)
                                            score_sum += current_alignment.score
                                        score_sums.append(score_sum)
                                    max_score = max(score_sums)
                                    if max_score > self.parameters.arguments["alignment_score_cutoff"]:
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
                            if len(conserved_path) / number_of_useqs >= self.parameters.arguments[
                                "orfs_presence_cutoff"] and len(conserved_path) > 1:
                                to_save_this_path = 1
                                for old_path in conserved_paths:
                                    fraction_of_identity = conserved_path.calculate_similarity(old_path)
                                    if fraction_of_identity >= self.parameters.arguments["paths_identity_cutoff"]:
                                        if conserved_path.score > old_path.score:
                                            conserved_paths.remove(old_path)
                                        elif conserved_path.score <= old_path.score:
                                            to_save_this_path = 0
                                if to_save_this_path == 1:
                                    #conserved_path.sort() # NOT SORTING!
                                    conserved_paths.append(conserved_path)
            self.conserved_paths = conserved_paths
            number_of_paths = len(conserved_paths)
            if number_of_paths == 0:
                print(f"⛔Termination:\n\tNo conserved ORFs set was found."
                      f"\n\tThis run will be terminated.", file=sys.stderr)
                if not os.path.exists(self.parameters.arguments["output_dir"]):
                    os.mkdir(self.parameters.arguments["output_dir"])
                with open(os.path.join(self.parameters.arguments["output_dir"], "report.txt"), "w") as report_f:
                    report_f.write("Termination:\nNo conserved ORFs set was found.")
                sys.exit()
            if self.parameters.arguments["verbose"]:
                print(f"✅ {number_of_paths} sets of conserved ORFs were found.",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to perform searching for conserved uORFs.") from error

    def filter_out_similar_paths(self) -> None:
        """Filter out duplicates in sets of annotated conserved ORFs.

        Note:
            Two paths are considered as duplicates if they share more than half of ORFs
                (default value, see 'paths_identity_cutoff' config parameter). In case two paths are found as identical,
                only one with a higher score will be saved.

        Returns:
            None

        """
        try:
            filtered_paths = []
            for path in self.conserved_paths:
                to_add = 1
                for path_filtered in filtered_paths:
                    if path.calculate_similarity(path_filtered) > self.parameters.arguments["paths_identity_cutoff"]:
                        if path.score < path_filtered.score:
                            to_add = 0
                        elif path.score == path_filtered.score and (len(path) < len(path_filtered)):
                            to_add = 0
                        else:
                            filtered_paths.remove(path_filtered)
                if to_add == 1:
                    filtered_paths.append(path)
            self.conserved_paths = filtered_paths

            if self.parameters.arguments["verbose"]:
                num_of_paths = len(self.conserved_paths)
                print(f"🧹 {num_of_paths} set(s) of conserved ORFs remained in the analysis after filtering "
                      f"out duplicates.", file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to filter out duplicates in conserved uORFs sets.") from error

    def run_msa(self) -> None:
        """Run msa tool (muscle) for each path object (set of conserved ORFs).

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🧮 Running MSA for conserved ORFs.", file=sys.stdout)
            for path in self.conserved_paths:
                path.maft_msa()
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to get MSA of conserved uORFs.") from error

    def save_msa(self) -> None:
        """Save MSA of conserved ORFs as fasta files.

        Note:
            Fasta files will be saved to the subdirs: ['nucleotide_msa' - for MSA of nucleotide sequences of ORFs,
                'amino_acid_msa' - MSA of amino acid sequences of ORFs, and 'sd_msa' - MSA of SD sequence regions
                of ORFS). All of them located in your 'output_dir'.

        Returns:
             None

        """
        try:
            if not os.path.exists(self.parameters.arguments["output_dir"]):
                os.mkdir(self.parameters.arguments["output_dir"])
            rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
            output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                                   [os.path.join(self.parameters.arguments["output_dir"],
                                                 f"{rename_dict[i]}_msa_fasta_files") for i in
                                    self.parameters.arguments['sequences_to_write']]))
            for key in output_dirs:
                if not (os.path.exists(output_dirs[key])):
                    os.mkdir(output_dirs[key])
            for path in self.conserved_paths:
                for seq_type in self.parameters.arguments["sequences_to_write"]:
                    msa = path.msa[seq_type]
                    output = os.path.join(output_dirs[seq_type], f"{path.id}.fa")
                    Bio.AlignIO.write(msa, output, "fasta")

            if self.parameters.arguments["verbose"]:
                output_dirs_v = [os.path.basename(i) for i in output_dirs.values()]
                delimiter = ",\n\t"
                print(f"💌 MSA fasta files of conserved ORFs were saved to the folders:\n"
                      f"\t{delimiter.join(output_dirs_v)} folders.", file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to save MSA of conserved uORFs.") from error

    def save_orfs_sequences(self) -> None:
        """Save sequences of conserved ORFs as fasta files.

        Note:
            Fasta files will be saved to the subdirs: ['nucleotide_orfs' - for MSA of nucleotide sequences of ORFs,
                'amino_acid_msa' - MSA of amino acid sequences of ORFs, and 'sd_msa' - MSA of SD sequence regions
                of ORFS). All of them located in your 'output_dir'.

        Returns:
             None

        """
        try:
            if not os.path.exists(self.parameters.arguments["output_dir"]):
                os.mkdir(self.parameters.arguments["output_dir"])
            rename_dict = dict(nt="nucleotide", aa="amino_acid")
            sequence_to_write = [i for i in self.parameters.arguments["sequences_to_write"] if i != "sd"]
            output_dirs = dict(zip(sequence_to_write, [os.path.join(self.parameters.arguments["output_dir"],
                                                                    f"{rename_dict[i]}_orfs_fasta_files") for i in
                                                       sequence_to_write]))
            for key in output_dirs:
                if not (os.path.exists(output_dirs[key])):
                    os.mkdir(output_dirs[key])
            for seq_type in sequence_to_write:
                for path in self.conserved_paths:
                    records = []
                    for orf in path.path:
                        if seq_type == "nt":
                            record = Bio.SeqRecord.SeqRecord(orf.nt_sequence, orf.id, "", orf.name)
                        if seq_type == "aa":
                            record = Bio.SeqRecord.SeqRecord(orf.aa_sequence, orf.id, "", orf.name)
                        records.append(record)
                    output = os.path.join(output_dirs[seq_type], f"{path.id}.fa")
                    Bio.SeqIO.write(records, output, "fasta")
            if self.parameters.arguments["verbose"]:
                delimiter = ",\n\t"
                output_dirs_v = [os.path.basename(i) for i in output_dirs.values()]
                print(f"💌 Sequences fasta files of conserved ORFs were saved to the folders: \n"
                      f"\t{delimiter.join(output_dirs_v)}.", file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to save sequences of conserved uORFs.") from error

    def save_results_summary_table(self) -> None:
        """Save results summary table.

        Note:
            A tsv table will be saved to your output_dir.

        Returns:
            None

        """
        try:
            colnames = "\t".join(
                ["id", "length", "average_distance_to_the_ORF", "aa_alignment_length", "nt_alignment_length", "score",
                 "number_of_orfs", "number_of_orfs/number_of_sequences", "consensus(aa)", "consensus(nt)",
                 "uORFs", "uORFs_annotations"])
            rows = [colnames]
            for path in self.conserved_paths:
                annotations = sorted(set([i.annotation for i in path.path]))
                if len(annotations) > 1 and "NA" in annotations:
                    pass
                    # annotations.remove("NA") # To check then
                row = "\t".join(
                    [path.id, str(path.length), str(statistics.mean([i.distance for i in path.path])),
                     str(path.msa["aa"].get_alignment_length()), str(path.msa["nt"].get_alignment_length()),
                     str(path.score), str(len(path)), str(round(len(path) / len(self.records), 3)),
                     str(path.msa_consensus["aa"]), str(path.msa_consensus["nt"]), ', '.join([i.id for i in path.path]),
                     ', '.join(annotations)])
                rows.append(row)
            output_file_path = os.path.join(self.parameters.arguments["output_dir"], "results_summary.tsv")
            f = open(output_file_path, "w")
            f.write("\n".join(rows))
            if self.parameters.arguments["verbose"]:
                print(f"💌 Results summary tsv table saved to: {os.path.basename(output_file_path)}.", file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to save results summary table.") from error

    def plot_msa_figs(self) -> None:
        """Plot MSA plots of  conserved ORFs

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🎨 Plotting MSA figures...", file=sys.stdout)
            for path in self.conserved_paths:
                path.plot_msa()

            if self.parameters.arguments["verbose"]:
                rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
                output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                                       [os.path.join(self.parameters.arguments["output_dir"],
                                                     f"{rename_dict[i]}_msa_visualisation") for i in
                                        self.parameters.arguments['sequences_to_write']]))
                output_dirs_v = [os.path.basename(i) for i in output_dirs.values()]
                delimiter = ",\n\t"
                print(f"💌 MSA figures were saved to the folders: \n\t{delimiter.join(output_dirs_v)}",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to plot sequence logo of conserved uORFs.") from error

    def plot_ggmsa_figs(self) -> None:
        """Plot MSA plots of conserved ORFs saved as fasta files.

        Note:
            R script based on ggmsa package [yulab-smu.top/ggmsa] used to produce MSA plots. R script (msa_plot.R)
                can be found in output_dir. This method uses subprocess to run this R script in the following way:
                `Rscript {output_dir}/msa_plot.R --msa_fasta path_to_fasta --output output_path --seq_type (nt/aa)
                --width N(mm) --height M(mm)`.
                Since during each run of uorf4u a local copy of this script is created
                in your output_dir, you can change it without any consequences for next uorf4u runs.
                This method based on _plot_ggmsa_ method of Path class and simply call it for each Path object.

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🎨 Plotting MSA figures...", file=sys.stdout)
            for path in self.conserved_paths:
                path.plot_ggmsa()

            if self.parameters.arguments["verbose"]:
                rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
                output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                                       [os.path.join(self.parameters.arguments["output_dir"],
                                                     f"{rename_dict[i]}_msa_visualisation") for i in
                                        self.parameters.arguments['sequences_to_write']]))
                output_dirs_v = [os.path.basename(i) for i in output_dirs.values()]
                delimiter = ",\n\t"
                print(f"💌 MSA figures were saved to the folders:\n\t{delimiter.join(output_dirs_v)}.",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to visualise MSA of conserved uORFs.") from error

    def plot_logo_figs(self) -> None:
        """Plot sequence Logo figures of conserved ORFs saved as fasta files.

        Note:
            This method uses logomaker package to produce images.

            This method based on _plot_logo_ method of Path class and simply call it for each Path object.

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🎨 Plotting sequence logo figures ...",
                      file=sys.stdout)

            for path in self.conserved_paths:
                path.plot_logo()

            if self.parameters.arguments["verbose"]:
                rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
                output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                                       [os.path.join(self.parameters.arguments["output_dir"],
                                                     f"{rename_dict[i]}_seqlogo_visualisation") for i in
                                        self.parameters.arguments['sequences_to_write']]))
                output_dirs_v = [os.path.basename(i) for i in output_dirs.values()]
                delimiter = ",\n\t"
                print(f"💌 Sequence logo figures were saved to the folders: \n\t{delimiter.join(output_dirs_v)}",
                      file=sys.stdout)
            return None
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to plot sequence logo of conserved uORFs.") from error

    def plot_annotation(self) -> None:
        """Plot loci' annotations figures with conserved ORFs highlighting.

        Returns:
            None

        """
        try:
            if self.parameters.arguments["verbose"]:
                print(f"🎨 Plotting loci annotations figures...",
                      file=sys.stdout)
            if not os.path.exists(self.parameters.arguments["output_dir"]):
                os.mkdir(self.parameters.arguments["output_dir"])
            output_dir = os.path.join(self.parameters.arguments["output_dir"], "annotation_visualisation")
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for path in self.conserved_paths:
                output_file_name = f"{os.path.join(output_dir, path.id)}.pdf"
                annotation_plot_manager = uorf4u.drawing_annotation.AnnotationPlotManager(path, self.records,
                                                                                          self.parameters)
                annotation_plot_manager.define_x_axis_coordinate_system()
                annotation_plot_manager.create_tracks()
                annotation_plot_manager.plot(output_file_name)
            if self.parameters.arguments["verbose"]:
                print(f"💌 Annotation figures were saved to the folder: {os.path.basename(output_dir)}",
                      file=sys.stdout)
        except Exception as error:
            raise uorf4u.manager.uORF4uError("Unable to plot loci' annotations figures.") from error


'''
    def plot_ggmsa(self) -> None:
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
'''


class ORF:
    """An ORF object holds information about an annotated ORF.

    Note:
        It's supposed that the ORFs class' objects will not be used directly by API users since
            it's only needed for other classes' methods.

    Attributes:
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        id (str): identifier of the ORF. Format: locus_id|accession_number|distance_from_the_start_codon_to_the_main_orf
        name (str): name of the ORF. Format: useq_name|distance_from_the_start_codon_to_the_main_orf
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

    def __init__(self, parameters: uorf4u.manager.Parameters, id: str, name: str, nt_sequence: Bio.Seq.Seq,
                 sd_window_seq: Bio.Seq.Seq, start: int, stop: int, distance: int, useq_index: int,
                 annotation: str = "NA"):
        """Create an ORF object.

        Arguments:
            parameters (uorf4u.manager.Parameters): Parameters' class object.
            id (str): identifier of the ORF. Format: locus_id:distance_from_the_start_codon_to_the_proteins_orf:length.
            nt_sequence (Bio.Seq.Seq): a Seq object of nucleotide sequence of the ORF.
            sd_window_seq (Bio.Seq.Seq): a Seq object of upstream sequence to the start codon of the ORF.
            start (int): start position of the ORF on the locus (0-based).
            stop (int): stop position of the ORF on the locus (0-based).
            distance (int): distance to the main ORF.

        """

        self.parameters = parameters
        codon_table = Bio.Data.CodonTable.unambiguous_dna_by_name[  # ambiguous can be needed!
            parameters.arguments["ncbi_genetic_code_name"]]
        codon_table_ambiguous = Bio.Data.CodonTable.ambiguous_dna_by_name[  # ambiguous can be needed!
            parameters.arguments["ncbi_genetic_code_name"]]
        self.name = name
        self.distance = distance
        self.id = id
        self.sequence_id = id.split(":")[0]
        self.start = start
        self.stop = stop
        self.length = len(nt_sequence)
        self.nt_sequence = nt_sequence
        self.annotation = annotation
        self.useq_index = useq_index
        try:
            self.aa_sequence = self.nt_sequence.translate(table=codon_table)
        except:
            self.aa_sequence = self.nt_sequence.translate(table=codon_table_ambiguous)
        self.sd_window_seq = sd_window_seq
        self.extended_orfs = []
        self.min_energy = 0
        self.putative_sd_sequence = "NA"
        self.sd_window_seq_str = "NA"
        self.sd_window_energies = []

    def calculate_energies(self) -> None:
        """Calculate energies of putative SD sequences of the upstream sequence.

        Returns:
            None

        """
        # Loading reference energies json file
        with open(self.parameters.arguments["ref_energies"]) as ref_energy_file:
            ref_energy = json.load(ref_energy_file)
        sd_seq_length = min([len(i) for i in ref_energy.keys()])
        # Energies calculations
        if len(self.sd_window_seq) >= min(ref_energy.values()):
            energies = []
            for position in range((len(self.sd_window_seq) - sd_seq_length) + 1):
                try:
                    energies.append(
                        ref_energy[self.sd_window_seq[position:position + sd_seq_length]])
                except:
                    energies.append(0)
            if energies:
                self.min_energy = min(energies)
                self.sd_window_energies = [str(i) for i in energies]
                if self.min_energy < self.parameters.arguments["sd_energy_cutoff"]:
                    sd_start_position = energies.index(self.min_energy)  # Be careful, it could be more than one!
                    self.putative_sd_sequence = self.sd_window_seq[sd_start_position:sd_start_position + sd_seq_length]
                    self.sd_window_seq_str = (f"{self.sd_window_seq[0:sd_start_position].lower()}"
                                              f"{self.putative_sd_sequence.upper()}"
                                              f"{self.sd_window_seq[sd_start_position + sd_seq_length:].lower()}")

        return None


class Path:
    """A Path object holds information about a list of conserved ORFs.

    Note:
        It's supposed that the Path class' objects will not be used directly by API users since
            it's only needed for other classes' methods.

    Attributes:
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        path (list): List of the ORF class objects.
        score (float): Score of the Path (calculated as sum of pairwise alignments scores of ORFs).
        msa (dict): Dict with Multiple sequence alignment (MSA, Bio.Align.MultipleSeqAlignment object) as values
            for different sequences (nt, aa, sd) as keys.
        msa_consensus (dict): Dict with consensus sequence (Bio.Seq.Seq object) as values
            for different sequences (nt, aa, sd) as keys.
        length: length of the nucleotide sequence alignment.
        id (str): Path's id (format: length|score|num_of_orfs|average_distance_to_the_main_ORF

    """

    def __init__(self, parameters: uorf4u.manager.Parameters):
        """Create a Path object.

        Arguments:
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        self.parameters = parameters
        self.path = []
        self.score = 0
        self.msa = dict()
        self.msa_consensus = dict()
        self.id = None
        self.length = None

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

    def sort(self) -> None:
        """Sort list of ORFs by their names.

        Returns:
            None

        """
        sorted_path = [x for _, x in sorted(zip([i.name for i in self.path], self.path), key=lambda pair: pair[0])]
        self.path = sorted_path

        return None

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
        self.msa = dict()
        for seq_type in self.parameters.arguments["sequences_to_write"]:
            records = []
            for orf in self.path:
                # record_id = f"{orf.id}"
                # record_description = f"{(orf.name.split('|')[0])}"
                record_id = f"{orf.name}"
                record_description = ""
                if seq_type == "nt":
                    record = Bio.SeqRecord.SeqRecord(orf.nt_sequence, record_id, "", record_description)
                elif seq_type == "aa":
                    record = Bio.SeqRecord.SeqRecord(orf.aa_sequence, record_id, "", record_description)
                elif seq_type == "sd":
                    record = Bio.SeqRecord.SeqRecord(orf.sd_window_seq, record_id, "", record_description)
                records.append(record)
            temp_input = tempfile.NamedTemporaryFile()
            Bio.SeqIO.write(records, temp_input.name, "fasta")
            temp_output = tempfile.NamedTemporaryFile()
            muscle = self.parameters.arguments["muscle_binary"]
            subprocess.run([muscle, "-align", temp_input.name, "-output", temp_output.name],
                           stderr=subprocess.DEVNULL)
            temp_input.close()
            msa = Bio.AlignIO.read(temp_output.name, "fasta")
            # msa.sort(key=lambda r: r.description)
            msa_info = Bio.Align.AlignInfo.SummaryInfo(msa)
            msa_consensus = msa_info.gap_consensus(threshold=self.parameters.arguments["consensus_threshold"])
            temp_output.close()
            if seq_type == "nt":
                self.length = msa.get_alignment_length()
            self.msa[seq_type], self.msa_consensus[seq_type] = msa, msa_consensus

            avr_distance = str(round(statistics.mean([i.distance for i in self.path])))
            self.id = f"length-{self.msa['nt'].get_alignment_length()}|score–{round(self.score)}|" \
                      f"num_of_orfs-{len(self.path)}|avr_dist-{avr_distance}"
        return None

    def maft_msa(self) -> None:
        """Run a multiple sequence alignment tool (MAFT) for the ORFs nucleotide and amino acid sequences.

        Note:
            This method updates nt_msa and aa_msa attributes.

        Returns:
            None

        """
        self.msa = dict()
        for seq_type in self.parameters.arguments["sequences_to_write"]:
            records = []
            for orf in self.path:
                # record_id = f"{orf.id}"
                # record_description = f"{(orf.name.split('|')[0])}"
                record_id = f"{orf.id}"
                record_description = orf.name
                if seq_type == "nt":
                    record = Bio.SeqRecord.SeqRecord(orf.nt_sequence, record_id, "", record_description)
                elif seq_type == "aa":
                    record = Bio.SeqRecord.SeqRecord(orf.aa_sequence, record_id, "", record_description)
                elif seq_type == "sd":
                    record = Bio.SeqRecord.SeqRecord(orf.sd_window_seq, record_id, "", record_description)
                records.append(record)
            temp_input = tempfile.NamedTemporaryFile()
            Bio.SeqIO.write(records, temp_input.name, "fasta")
            temp_output = tempfile.NamedTemporaryFile()
            temp_stderr = tempfile.NamedTemporaryFile()
            maft = self.parameters.arguments["maft_binary"]
            try:
                subprocess.run([maft, "--auto", "--reorder", temp_input.name], stdout=temp_output, stderr=temp_stderr)
                msa = Bio.AlignIO.read(temp_output.name, "fasta")
                temp_stderr.close()
                temp_output.close()
            except Exception as error:
                temp_stderr.seek(0)
                temp_output.seek(0)
                print(f"🤬 MAFFT error message:\n{temp_stderr.read()}", file=sys.stderr)
                temp_stderr.close()
                temp_output.close()
                raise uorf4u.manager.uORF4uError(f"mafft error. If you work on a linux machine,"
                                                 f" run uorf4 --linux.") from error
            for record in msa:
                record.description = " ".join(record.description.split(" ")[1:])
            # msa.sort(key=lambda r: r.description) # add a parameter for order setting
            msa_info = Bio.Align.AlignInfo.SummaryInfo(msa)
            msa_consensus = msa_info.gap_consensus(threshold=self.parameters.arguments["consensus_threshold"])
            temp_output.close()
            if seq_type == "nt":
                self.length = msa.get_alignment_length()
            self.msa[seq_type], self.msa_consensus[seq_type] = msa, msa_consensus
            avr_distance = str(round(statistics.mean([i.distance for i in self.path])))
            self.id = f"length-{self.msa['nt'].get_alignment_length()}|score–{round(self.score)}|" \
                      f"num_of_orfs-{len(self.path)}|avr_dist-{avr_distance}"
        return None

    def plot_msa(self) -> None:
        """Plot MSA of conserved ORFs.

        Returns:
            None

        """

        rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
        output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                               [os.path.join(self.parameters.arguments["output_dir"],
                                             f"{rename_dict[i]}_msa_visualisation") for i in
                                self.parameters.arguments["sequences_to_write"]]))
        for o_dir in output_dirs.values():
            if not (os.path.exists(o_dir)):
                os.mkdir(o_dir)

        for s_type in self.parameters.arguments["sequences_to_write"]:
            current_msa = self.msa[s_type]
            if s_type == "nt" or s_type == "sd":
                seq_type = "nt"
            else:
                seq_type = "aa"
            msa_plot_manager = uorf4u.drawing_msa.MSAPlotManager(current_msa, self.parameters, seq_type)
            msa_plot_manager.define_x_axis_coordinate_system()
            output_file = os.path.join(output_dirs[s_type], f"{self.id}.pdf")
            msa_plot_manager.create_tracks()
            msa_plot_manager.plot(output_file)

    def plot_ggmsa(self) -> None:
        """Plot MSA of conserved ORFs saved as fasta files.

        Note:
            R script based on ggmsa package [yulab-smu.top/ggmsa] used to produce MSA plots. R script (msa_plot.R)
                can be found in output_dir. This method uses subprocess to run this R script in the following way:
                `Rscript {output_dir}/msa_plot.R --msa_fasta path_to_fasta --output output_path --seq_type (nt/aa)
                --width N(mm) --height M(mm)`.
                Since during each run of uorf4u a local copy of this script is created
                in your output_dir, you can change it without any consequences for next uorf4u runs.

        Returns:
            None

        """

        rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
        output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                               [os.path.join(self.parameters.arguments["output_dir"],
                                             f"{rename_dict[i]}_msa_visualisation") for i in
                                self.parameters.arguments["sequences_to_write"]]))
        fasta_files_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                                    [os.path.join(self.parameters.arguments["output_dir"],
                                                  f"{rename_dict[i]}_msa_fasta_files") for i in
                                     self.parameters.arguments["sequences_to_write"]]))
        for o_dir in output_dirs.values():
            if not (os.path.exists(o_dir)):
                os.mkdir(o_dir)
        r_script_path = self.parameters.arguments["plot_msa_R_script"]
        r_script_local = os.path.join(self.parameters.arguments["output_dir"], os.path.basename(r_script_path))
        if not (os.path.exists(r_script_local)):
            shutil.copy(r_script_path, r_script_local)
        for s_type in self.parameters.arguments["sequences_to_write"]:
            current_msa = self.msa[s_type]
            if s_type == "nt" or s_type == "sd":
                seq_type = "nt"
            else:
                seq_type = "aa"

            output_dir = os.path.abspath(os.path.join(output_dirs[s_type]))
            input_file = os.path.abspath(os.path.join(fasta_files_dirs[s_type], f"{self.id}.fa"))
            num_sequences = len(current_msa)
            length_of_alignment = current_msa.get_alignment_length()
            page_width = (50 + length_of_alignment) * 5
            page_height = max(17, (num_sequences + 5) * 3)
            subprocess.run(["Rscript", r_script_local, "--msa_fasta", input_file, "--output", output_dir,
                            "--seq_type", seq_type, "--width", str(page_width), "--height", str(page_height)])

    def plot_logo(self) -> None:
        """Plot sequence Logo of conserved ORFs MSA saved as fasta files.

        Note:
            This method uses logomaker package to produce images.

        Returns:
            None

        """
        rename_dict = dict(nt="nucleotide", aa="amino_acid", sd="sd")
        output_dirs = dict(zip(self.parameters.arguments["sequences_to_write"],
                               [os.path.join(self.parameters.arguments["output_dir"],
                                             f"{rename_dict[i]}_seqlogo_visualisation") for i in
                                self.parameters.arguments['sequences_to_write']]))

        for o_dir in output_dirs.values():
            if not (os.path.exists(o_dir)):
                os.mkdir(o_dir)
        ambiguous_codon_table = Bio.Data.CodonTable.ambiguous_dna_by_name[
            self.parameters.arguments["ncbi_genetic_code_name"]]
        unambiguous_codon_table = Bio.Data.CodonTable.unambiguous_dna_by_name[
            self.parameters.arguments["ncbi_genetic_code_name"]]
        alphabet = dict(nt=set(ambiguous_codon_table.nucleotide_alphabet),
                        aa=set(ambiguous_codon_table.protein_alphabet))
        unambiguous_alphabet = dict(nt=set(unambiguous_codon_table.nucleotide_alphabet),
                                    aa=set(unambiguous_codon_table.protein_alphabet))
        for s_type in self.parameters.arguments["sequences_to_write"]:
            current_msa = self.msa[s_type]
            if s_type == "nt" or s_type == "sd":
                seq_type = "nt"
            elif s_type == "aa":
                seq_type = "aa"
            output_file = os.path.abspath(
                os.path.join(output_dirs[s_type], os.path.basename(self.id)))
            msa_length = current_msa.get_alignment_length()
            num_of_sequences = len(current_msa)
            current_msa_info = Bio.Align.AlignInfo.SummaryInfo(current_msa)
            pos_specific_dict = dict()
            pos_specific_score_matrix = current_msa_info.pos_specific_score_matrix()
            for i in alphabet[seq_type]:
                pos_specific_dict[i] = [0 for j in range(msa_length)]
            for i in range(msa_length):
                for element in pos_specific_score_matrix[i].keys():
                    pos_specific_dict[element.upper()][i] = (pos_specific_score_matrix[i][element] / num_of_sequences)
            pos = [i for i in range(msa_length)]
            pos_specific_dict = {k: v for k, v in pos_specific_dict.items() if
                                 sum(v) > 0 or k in unambiguous_alphabet[seq_type]}
            matrix_fr = pandas.DataFrame(pos_specific_dict, index=pos)
            colors = self.parameters.arguments[f"colors_{seq_type}"]
            colors = {k: uorf4u.methods.color_name_to_hex(v, self.parameters.arguments) for k, v in colors.items()}
            fig_size = (min(max(10, msa_length * 1.3), ((2 ** 16) - 1) / 100),
                        min(2.5, 2.5 * 10 / (msa_length ** (1 / 5))))

            if self.parameters.arguments["logo_type"] == "probability" or \
                    self.parameters.arguments["logo_type"] == "both":
                output_file_fr = f"{output_file}_prob.pdf"
                max_value_fr = 1
                logo_fr = logomaker.Logo(matrix_fr, color_scheme=colors, figsize=fig_size,
                                         alpha=self.parameters.arguments["logo_alpha"], show_spines=False,
                                         baseline_width=0)
                logo_fr.style_spines(spines=["left"], visible=True, linewidth=0.7)
                logo_fr.ax.set_xticks([])
                logo_fr.ax.set_yticks([0, max_value_fr])
                plt.savefig(output_file_fr)
                plt.close("all")

            if self.parameters.arguments["logo_type"] == "information" or \
                    self.parameters.arguments["logo_type"] == "both":
                colors["-"] = colors["_"]
                matrix_fr["-"] = round((1 - matrix_fr.sum(axis=1)), 5)
                if matrix_fr["-"].sum() == 0:
                    del matrix_fr['-']
                matrix_info = logomaker.transform_matrix(matrix_fr, from_type="probability", to_type="information")
                max_value_info = math.log2(len(pos_specific_dict.keys()))
                output_file_info = f"{output_file}_info.pdf"
                logo_info = logomaker.Logo(matrix_info, color_scheme=colors, figsize=fig_size,
                                           alpha=self.parameters.arguments["logo_alpha"], show_spines=False,
                                           baseline_width=0)
                logo_info.style_spines(spines=["left"], visible=True, linewidth=0.7)
                logo_info.ax.set_xticks([])
                logo_info.ax.set_yticks([0, max_value_info])
                plt.savefig(output_file_info)
                plt.close("all")
        return None
