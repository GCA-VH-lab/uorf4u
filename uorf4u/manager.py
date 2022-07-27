import traceback
import configs
import argparse
import sys
import os
import time


class Ant4suorfError(Exception):
    """A helper for exceptions parsing inherited from the Exception class.

    """
    pass


class Parameters:
    """A Parameters object holds and parse cmd's and config's arguments for the tool.

    A Parameters object have to be created in each script since it's used by each
        class of the tool as a mandatory argument.

    """

    def __init__(self):
        self.arguments = dict()

    def parse_cmd_arguments(self) -> None:
        parser = argparse.ArgumentParser(prog="uorf4u", add_help=False,
                                         usage="uorf4u [-an accession_number| -hl [ac1, ac2..] | -hlf path ]"
                                               "[optional arguments]")

        mutually_exclusive_group = parser.add_mutually_exclusive_group()
        mutually_exclusive_group.add_argument("-an", dest="accession_number", type=str, default=None)
        mutually_exclusive_group.add_argument("-hl", dest="homologous_list", nargs="*", default=None)
        mutually_exclusive_group.add_argument("-hlf", dest="homologous_list_file", type=str, default=None)
        # mutually_exclusive_group.add_argument('-useq', dest='upstream_sequences', type=str, default=None)

        parser.add_argument("-sao", dest="save_annotated_orfs", action="store_true")
        parser.add_argument("-o", dest="output_dir", type=str, default=None)
        parser.add_argument("-c", dest="config_file", type=str, default="internal")
        parser.add_argument("-v", "--version", action='version', version='%(prog)s 0.1')
        parser.add_argument("-h", "--help", dest="help", action="store_true")
        args = parser.parse_args()
        args = vars(args)
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
        except Exception as error:
            raise Ant4suorfError(
                "Unable to parse the specified config file. Please check your config file.") from error

    def update(self, parameters):
        self.arguments.update(parameters)
