import shutil
import os


def copy_package_data() -> None:
    """Copy uorf4u package data folder to your current dir.

    Returns:
        None
    """
    users_dir = os.path.join(os.getcwd(), 'uorf4u_data')
    internal_dir = os.path.join(os.path.dirname(__file__), 'uorf4u_data')
    shutil.copytree(internal_dir, users_dir, ignore = shutil.ignore_patterns("help*", ".*", "msa_plot_dir.R"))
    return None