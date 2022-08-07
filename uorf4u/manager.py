import traceback
import configs
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

    def parse_cmd_arguments(self) -> None:
        parser = argparse.ArgumentParser(prog="uorf4u", add_help=False,
                                         usage="uorf4u [-an accession_number| -hl [ac1, ac2..] | -hlf path ]"
                                               "[optional arguments]")

        mutually_exclusive_group = parser.add_mutually_exclusive_group()
        mutually_exclusive_group.add_argument("-an", dest="accession_number", type=str, default=None)
        mutually_exclusive_group.add_argument("-hl", dest="homologous_list", nargs="*", default=None)
        mutually_exclusive_group.add_argument("-hlf", dest="homologous_list_file", type=str, default=None)
        # mutually_exclusive_group.add_argument('-useq', dest='upstream_sequences', type=str, default=None)
        parser.add_argument("--data", dest="uorf4u_data", action="store_true")
        parser.add_argument("-al", dest="assemblies_list", type=str, default="NA")
        parser.add_argument("-at", dest="alignment_type", choices=['nt', 'aa', None], type=str, default=None)
        parser.add_argument("-o", dest="output_dir", type=str, default=None)
        parser.add_argument("-c", dest="config_file", type=str, default="internal")
        parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.3.0')
        parser.add_argument("--verbose", "-verbose", dest="verbose", action='store_true')
        parser.add_argument("--debug", "-debug", dest="debug", action="store_true")
        parser.add_argument("-h", "--help", dest="help", action="store_true")

        args = parser.parse_args()
        args = vars(args)

        if len(sys.argv[1:]) == 0:
            args["help"] = True

        if args["uorf4u_data"]:
            print("a")
            uorf4u.methods.copy_package_data()
            sys.exit()

        if args["help"]:
            help_message_path = os.path.join(os.path.dirname(__file__), 'uorf4u_data', "help.txt")
            with open(help_message_path, "r") as help_message:
                print(help_message.read(), file=sys.stdout)
                sys.exit()

        filtered_args = {k: v for k, v in args.items() if v is not None}
        self.arguments.update(filtered_args)

    def load_config(self, path="internal"):
        try:
            if path == "internal":
                path = os.path.join(os.path.dirname(__file__), "uorf4u_data", "uorf4u.cfg")
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
            self.load_palette()
        except Exception as error:
            raise uORF4uError(
                "Unable to parse the specified config file. Please check your config file.") from error

    def load_palette(self):
        for seq_type in ["nt", "aa"]:
            palette_path = self.arguments[f"palette_{seq_type}"]
            palette_pre_dict = configs.load(palette_path).get_config()["root"]
            palette_dict = dict()
            for elements, color in palette_pre_dict.items():
                for element in elements:
                    palette_dict[element] = color
            self.arguments[f"palette_{seq_type}"] = palette_dict

    def update(self, parameters):
        self.arguments.update(parameters)
