"""
This module provides some methods (e.g. colors tranformation, data copying) used by the tool.
"""
import shutil
import os
import uorf4u.manager
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont


def copy_package_data() -> None:
    """Copy uorf4u package data folder to your current dir.

    Returns:
        None

    """
    try:
        users_dir = os.path.join(os.getcwd(), 'uorf4u_data')
        internal_dir = os.path.join(os.path.dirname(__file__), 'uorf4u_data')
        shutil.copytree(internal_dir, users_dir, ignore=shutil.ignore_patterns("help*", ".*", "msa_plot_dir.R"))
        return None
    except Exception as error:
        raise uorf4u.manager.uORF4uError(f"Unable to copy uorf4u_data folder in your working dir.") from error


def get_color(name: str, parameters: dict) -> tuple:
    return *hex_to_rgb(parameters[name]), parameters[f"{name}_alpha"]


def hex_to_rgb(value: str) -> list:
    """Convert HEX color to RGB format.

    Returns:
        list: color in rgb format

    """
    try:
        value = value.lstrip('#')
        lv = len(value)
        rgb = [i / 255 for i in tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))]
        return rgb
    except Exception as error:
        raise uorf4u.manager.uORF4uError(
            f"Unable to convert color definition from HEX to RGB. Please check the palette config file.") from error


def string_height_to_font_size(height: float, font_type: str, parameters : dict):
    """

    Arguments:
        height (float): available height of the string.
        font_type (str): font type (see config file; at this moment only regular is available)
        parameters (uorf4u.manager.Parameters): Parameters' class object.

    Returns:
        float: font size defined by height

    """
    pdfmetrics.registerFont(TTFont(font_type, parameters[f"font_{font_type}"]))
    face = pdfmetrics.getFont('regular').face
    font_size = (1000 * 1.38 * height) / (face.ascent - face.descent)
    return font_size
