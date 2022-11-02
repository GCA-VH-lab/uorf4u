"""
This module provides some methods (e.g. colors tranformation, data copying) used by the tool.
"""
import shutil
import os
import re
import uorf4u.manager
import Bio.SeqIO
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


def adjust_paths_for_linux() -> None:
    """Change paths in the internal config files for linux.

    Returns:
        None

    """
    internal_dir = os.path.join(os.path.dirname(__file__), "uorf4u_data")
    config_files = ["uorf4u_eukaryotes.cfg", "uorf4u_bacteria.cfg"]
    for config_file in config_files:
        config_file_path = os.path.join(internal_dir, config_file)
        with open(config_file_path, "r+") as config:
            config_txt = re.sub(r"/muscle5\.1\.macos_arm64", "/muscle5.1.linux_intel64", config.read())
            config_txt = re.sub(r"/mafft-mac/mafft\.bat", "/mafft-linux64/mafft.bat", config_txt)
            config.seek(0)
            config.truncate()
            config.write(config_txt)
    return None


def copy_package_data() -> None:
    """Copy the uorf4u package data folder to your current dir.

    Returns:
        None

    """
    try:
        users_dir = os.path.join(os.getcwd(), "uorf4u_data")
        internal_dir = os.path.join(os.path.dirname(__file__), "uorf4u_data")
        shutil.copytree(internal_dir, users_dir, ignore=shutil.ignore_patterns("help*", ".*", "msa_plot_dir.R"))
        return None
    except Exception as error:
        raise uorf4u.manager.uORF4uError(f"Unable to copy uorf4u_data folder in your working dir.") from error


def hex_to_rgb(value: str) -> list:
    """Convert HEX color to RGB format.

    Arguments:
        value (str): color in HEX format.

    Returns:
        list: RGB color.

    """
    try:
        value = value.lstrip("#")
        lv = len(value)
        rgb = [i / 255 for i in tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))]
        return rgb
    except Exception as error:
        raise uorf4u.manager.uORF4uError(
            f"Unable to convert color definition from HEX to RGB. Please check the palette config file.") from error


def color_name_to_hex(name: str, parameters: dict) -> str:
    return parameters["palette"][name]


def get_color(name: str, parameters: dict) -> tuple:
    """Get color code by a name.

    Arguments:
        name (str): name of a color.
        parameters (dict): Parameters' object dict.

    Returns:
        tuple: RGB color.

    """
    rgb_color = *hex_to_rgb(parameters['palette'][parameters[name]]), parameters[f"{name}_alpha"]
    return rgb_color


def string_height_to_font_size(height: float, font_type: str, parameters: dict) -> float:
    """Transform string height to the font size.

    Arguments:
        height (float): available height of the string.
        font_type (str): font type (see config file; at this moment only regular is available).
        parameters (dict): Parameters' object dict.

    Returns:
        float: font size defined by height.

    """
    pdfmetrics.registerFont(TTFont(font_type, parameters[f"font_{font_type}"]))
    face = pdfmetrics.getFont('regular').face
    font_size = (1000 * 1.38 * height) / (face.ascent - face.descent)
    return font_size


def parse_fasta_file(path: str, parameters) -> list:
    """Parse fasta file with sequences.

    Arguments
        path: path to a fasta file.

    Returns:
        list: list of processed Bio.SeqRecord.SeqRecord objects.

    """
    try:
        processed_records = []
        record_ids = []
        with open(path) as handle:
            for record in Bio.SeqIO.parse(handle, "fasta"):
                label = record.description
                record.description = record.description.replace(record.id, "")
                length = len(record.seq)
                record_annotation = dict(RefSeq=False, length=length,
                                         upstream_region_length=length - parameters.arguments[
                                             "downstream_region_length"],
                                         downstream_region_length=parameters.arguments["downstream_region_length"],
                                         label=label)
                record.annotations = record_annotation
                processed_records.append(record)
        return processed_records
    except Exception as error:
        raise uorf4u.manager.uORF4uError(f"Unable to process the fasta file with sequences.") from error
