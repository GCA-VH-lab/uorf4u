"""
This module provides managing classes and methods for the tool.
"""
import configs
import traceback
import argparse
import sys
import os
import time
import uorf4u.methods


class uORF4uError(Exception):
    """A helper for exceptions parsing inherited from the Exception class.

    """
    pass


class Parameters:
    """A Parameters object holds and parse cmd's and config's arguments for the tool.

    Note:
        A Parameters object have to be created in each script since it's used by each
            class of the tool as a mandatory argument.

    """

    def __init__(self):
        self.arguments = dict(assemblies_list="NA", debug=False, verbose=False)
        self.cmd_arguments = {}

    def parse_cmd_arguments(self) -> None:
        parser = argparse.ArgumentParser(prog="uorf4u", add_help=False,
                                         usage="uorf4u [-an accession_number | -hl [ac1, ac2..] | -hlf path | -fa path]"
                                               "[optional arguments]")
        mutually_exclusive_group = parser.add_mutually_exclusive_group()
        mutually_exclusive_group.add_argument("-an", dest="accession_number", type=str, default=None)
        mutually_exclusive_group.add_argument("-hl", dest="homologues_list", nargs="*", default=None)
        mutually_exclusive_group.add_argument("-hlf", dest="homologues_list_file", type=str, default=None)
        mutually_exclusive_group.add_argument("-fa", dest="fasta", type=str, default=None)
        parser.add_argument("-data", "--data", dest="uorf4u_data", action="store_true")
        parser.add_argument("-linux", "--linux", dest="linux", action="store_true", default=None)
        parser.add_argument("-bh", dest="blastp_hit_list_size", type=int, default=None)
        parser.add_argument("-bid", dest="blastp_pident_to_query_length_cutoff", type=float, default=None)
        parser.add_argument("-mna", dest="max_number_of_assemblies", type=int, default=None)
        parser.add_argument("-al", dest="assemblies_list", type=str, default="NA")
        parser.add_argument("-annot", dest="check_assembly_annotation", action="store_true", default=None)
        parser.add_argument("-ul", dest="upstream_region_length", type=int, default=None)
        parser.add_argument("-dl", dest="downstream_region_length", type=int, default=None)
        parser.add_argument("-asc", dest="alternative_start_codons", action="store_true", default=None)
        parser.add_argument("-nsd", dest="filter_by_sd", action="store_false", default=None)
        parser.add_argument("-at", dest="alignment_type", choices=['nt', 'aa', None], type=str, default=None)
        parser.add_argument("-fast", dest="fast_searching", action="store_true", default=None)
        parser.add_argument("-o", dest="output_dir", type=str, default=None)
        parser.add_argument("-c", dest="config_file", type=str, default="prokaryotes")
        parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.6.2')
        parser.add_argument("-q", "--quiet", dest="verbose", default=True, action="store_false")
        parser.add_argument("--debug", "-debug", dest="debug", action="store_true")
        parser.add_argument("-h", "--help", dest="help", action="store_true")
        args = parser.parse_args()
        args = vars(args)

        if len(sys.argv[1:]) == 0:
            args["help"] = True

        if args["uorf4u_data"]:
            uorf4u.methods.copy_package_data()
            sys.exit()

        if args["linux"]:
            uorf4u.methods.adjust_paths_for_linux()
            sys.exit()

        if args["help"]:
            help_message_path = os.path.join(os.path.dirname(__file__), 'uorf4u_data', "help.txt")
            with open(help_message_path, "r") as help_message:
                print(help_message.read(), file=sys.stdout)
                sys.exit()

        filtered_args = {k: v for k, v in args.items() if v is not None}
        self.cmd_arguments = filtered_args

    def load_config(self, path="prokaryotes"):
        try:
            if path == "prokaryotes" or path == "eukaryotes":
                path = os.path.join(os.path.dirname(__file__), "uorf4u_data", f"uorf4u_{path}.cfg")
            config = configs.load(path)
            config = config.get_config()
            internal_dir = os.path.dirname(__file__)
            config["root"]["output_dir"] = config["root"]["output_dir"].replace("{current_date}",
                                                                                time.strftime("%Y_%m_%d-%H_%M"))
            for key in config["root"].keys():
                if type(config["root"][key]) is str and "{internal}" in config["root"][key]:
                    config["root"][key] = config["root"][key].replace("{internal}",
                                                                      os.path.join(internal_dir, "uorf4u_data"))
            self.arguments.update(config['root'])
            self.arguments.update(self.cmd_arguments)
            self.load_palette()
            self.load_color_config()
        except Exception as error:
            raise uORF4uError(
                "Unable to parse the specified config file. Please check your config file or written name.") from error

    def load_palette(self) -> None:
        palette_path = self.arguments[f"palette"]
        self.arguments[f"palette"] = configs.load(palette_path).get_config()["root"]

    def load_color_config(self) -> None:
        for seq_type in ["nt", "aa"]:
            path = self.arguments[f"colors_{seq_type}"]
            colors_pre_dict = configs.load(path).get_config()["root"]
            colors_dict = dict()
            for elements, color in colors_pre_dict.items():
                for element in elements:
                    colors_dict[element] = color
            self.arguments[f"colors_{seq_type}"] = colors_dict

    def update(self, parameters):
        self.arguments.update(parameters)
