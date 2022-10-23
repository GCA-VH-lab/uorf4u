"""
This module provides visualisation of loci annotation.
"""
import reportlab.pdfgen
import reportlab.pdfgen.canvas
from reportlab.lib.units import cm, mm
import reportlab.rl_config
import reportlab.pdfbase.ttfonts
import reportlab.pdfbase.pdfmetrics
import Bio.Seq

reportlab.rl_config.warnOnMissingFontGlyphs = 0

import uorf4u.methods
import uorf4u.manager


class MSAPlotManager:
    """
    AnnotationPlotManager object holds needed information for annotation visualisation and controls it.

    Note:
        It's supposed that the AnnotationPlotManager' objects will not be used directly by API users since
            visualisation can be controlled by 'plot_annotation' method.

    Attributes:
        msa (FILL IN): Path class' multiple sequence alignment.
        upstream_sequences (list): list of dicts with information about upstream sequences.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        coordinate_system (dict): coordinate system of figure.
        additional_data (dict): dict with data for visualisation tracks.

    """

    def __init__(self, msa, parameters: uorf4u.manager.Parameters, type: str):
        """Create a AnnotationPlotManager object.

        Arguments:
            path (uorf4u.data_processing.Path): Path class' objects that holds list of conserved ORFs.
            upstream_sequences (list): list of dicts with information about upstream sequences.
            parameters (uorf4u.manager.Parameters): Parameters' class object.
            type (str): type of sequences (sd, nt, aa)

        """
        self.msa = msa
        self.parameters = parameters
        self.coordinate_system = dict()
        self.additional_data = dict()
        self.type = type

    def define_x_axis_coordinate_system(self) -> None:
        """Define coordinate system.

        Returns:
            None

        """

        label_height = self.parameters.arguments["label_size"] * self.parameters.arguments["tile_size"] * cm
        label_font_size = uorf4u.methods.string_height_to_font_size(label_height, "regular", self.parameters.arguments)
        self.additional_data["label_font_size"] = label_font_size
        msa_length = self.msa.get_alignment_length()
        max_label_width = max([reportlab.pdfbase.pdfmetrics.stringWidth(i.description, "regular",
                                                                        label_font_size) for i in self.msa])

        char_height = self.parameters.arguments["char_size"] * self.parameters.arguments["tile_size"] * cm
        char_font_size = uorf4u.methods.string_height_to_font_size(char_height, "mono", self.parameters.arguments)
        self.additional_data["char_font_size"] = char_font_size

        self.additional_data["number_of_sequences"] = len(self.msa)
        self.coordinate_system["x_labels_start"] = self.parameters.arguments["margin"] * cm
        self.coordinate_system["x_labels_stop"] = self.coordinate_system["x_labels_start"] + max_label_width
        self.coordinate_system["x_msa_start"] = self.coordinate_system["x_labels_stop"] + \
                                                self.parameters.arguments["label_gap"] * cm
        msa_width = self.parameters.arguments["tile_size"] * msa_length * cm
        self.coordinate_system["x_msa_stop"] = self.coordinate_system["x_msa_start"] + msa_width
        self.coordinate_system["figure_width"] = 2 * self.parameters.arguments["margin"] * cm + msa_width + \
                                                 max_label_width + self.parameters.arguments["label_gap"] * cm
        self.coordinate_system["figure_height"] = self.parameters.arguments["margin"] * cm
        self.additional_data["palette"] = self.parameters.arguments[f"colors_{self.type}"]
        self.additional_data["palette"] = {k: uorf4u.methods.color_name_to_hex(v, self.parameters.arguments) for k, v in
                                           self.additional_data["palette"].items()}

        return None

    def create_tracks(self) -> None:
        """Create visualisation tracks.

        Returns:
            None

        """
        self.tracks = []
        """
        title_loader = TitleLoader(self.parameters)
        title_loader.prepare_data(self.coordinate_system, self.additional_data)
        title_track = title_loader.create_track()
        self.tracks.append(title_track)
        self.coordinate_system["figure_height"] += title_track.needed_y_space()
        """
        for record in self.msa:
            sequence_loader = SequencesLoader(self.parameters)
            sequence_loader.prepare_data(record, self.coordinate_system, self.additional_data)
            track = sequence_loader.create_track()
            self.tracks.append(track)
            self.coordinate_system["figure_height"] += track.needed_y_space()
            # self.coordinate_system["figure_height"] += self.parameters.arguments["gap"] * cm
            # if index < self.additional_data["number_of_sequences"] - 1:
        self.coordinate_system["figure_height"] += self.parameters.arguments["margin"] * cm

    def plot(self, filename):
        image = Image(filename, self.coordinate_system["figure_width"], self.coordinate_system["figure_height"])
        current_y_top = self.coordinate_system["figure_height"] - self.parameters.arguments["margin"] * cm
        for track in self.tracks:
            track.visualisation_data["y_top"] = current_y_top
            track.draw(image.canvas)
            current_y_top -= (track.needed_space)
        image.save()
        return None


class Track:
    """Parent clas for visualisation Tracks.

    Attributes:
        visualisation_data (dict): a dictionary with data needed for visualisation.
        parameters (uorf4u.manager.Parameters): Parameters' class object.

    """

    def __init__(self, visualisation_data: dict, parameters: uorf4u.manager.Parameters):
        """Parent's constructor for creating a Track object.

        Arguments:
            visualisation_data (dict): a dictionary with data needed for visualisation.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        self.visualisation_data = visualisation_data
        self.parameters = parameters

    def needed_y_space(self) -> None:
        """Empy parent's method for calculation needed vertical space for a track.

        Returns:
            None

        """
        pass

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Empy parent's method for track visualisation.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None
        """
        pass


class SequenceVis(Track):
    """SequenceVis track draws sequences and annotation.

    Attributes:
        visualisation_data (dict): a dictionary with data needed for visualisation.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        needed_space (float): needed vertical space for a track.

    """

    def __init__(self, visualisation_data: dict, parameters: uorf4u.manager.Parameters):
        """Create a SequenceVis object.

        Arguments:
            visualisation_data (dict): a dictionary with data needed for visualisation.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        super().__init__(visualisation_data, parameters)
        self.needed_space = None

    def needed_y_space(self) -> float:
        """Calculate needed vertical space for a SequenceVis track.

        Returns:
            float: needed vertical space.

        """
        self.needed_space = self.parameters.arguments["tile_size"] * cm
        return self.needed_space

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a Sequence track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None

        """
        tile_size = self.parameters.arguments["tile_size"] * cm
        y_c = self.visualisation_data["y_top"] - (tile_size * 0.5)
        y_l = self.visualisation_data["y_top"] - tile_size

        y_gap_label = tile_size * (1 - self.parameters.arguments["label_size"]) * 0.5
        y_gap_char = tile_size * (1 - self.parameters.arguments["char_size"]) * 0.5
        # Labels
        canvas.setFillColorRGB(*uorf4u.methods.get_color("label_color", self.parameters.arguments))
        canvas.setFont("regular", self.visualisation_data["label_font_size"])
        canvas.drawRightString(self.visualisation_data["coordinate_system"]["x_labels_stop"], y_l + y_gap_label,
                               self.visualisation_data["label"])

        canvas.setLineWidth(0.05 * tile_size)
        canvas.setStrokeColorRGB(1, 1, 1)
        canvas.setFont("mono", self.visualisation_data["char_font_size"])
        x = self.visualisation_data['msa_left_border']
        for symbol in self.visualisation_data["sequence"]:
            x_c = x + tile_size * 0.5
            symbol = symbol.upper()
            try:
                color = self.visualisation_data["palette"][symbol]
            except:
                color = "#FFFFFF"
            canvas.setFillColorRGB(*uorf4u.methods.hex_to_rgb(color), self.parameters.arguments["tile_alpha"])
            canvas.rect(x, y_l, tile_size, tile_size, fill=1)
            canvas.setFillColorRGB(0, 0, 0, 0.8)  # to change
            canvas.drawCentredString(x_c, y_l + y_gap_char, symbol)
            x += tile_size

        return None


class Loader:
    """Parent class for tracks loaders.

    Attributes:
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        prepared_data (dict): dict with data needed for visualisation tracks.

    """

    def __init__(self, parameters: uorf4u.manager.Parameters):
        """Parent's constructor for creating a Loader class object.

        Arguments:
            parameters (uorf4u.manager.Parameters): Parameters' class object.
        """
        self.parameters = parameters
        self.prepared_data = None

    def prepare_data(self) -> None:
        """Empty parent's method for data preparation.

        Returns:
            None

        """
        pass

    def create_track(self) -> None:
        """Empty parent's method for initialisation of a track.

        Returns:
            None

        """
        pass


class SequencesLoader(Loader):
    """A SequencesLoader object prepares data for a Sequence track object.


    Attributes:
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        prepared_data (dict): dict with data needed for a visualisation track.

    """

    def __init__(self, parameters):
        """Create a SequenceLoader object.

        Arguments:
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        super().__init__(parameters)

    def prepare_data(self, record, coordinate_system: dict, additional_data: dict) -> dict:
        """Prepare data for a Title visualisation track.

        Attributes:
            record (FILL in): record of blablabla
            coordinate_system (dict): coordinate system of a figure page.
            additional_data (dict): data needed for a track initialisation.

        Returns:
            dict: dictionary with prepared data for visualisation.

        """
        prepared_data = dict()
        prepared_data["coordinate_system"] = coordinate_system
        prepared_data["label_font_size"] = additional_data["label_font_size"]
        prepared_data["char_font_size"] = additional_data["char_font_size"]
        prepared_data["label_right_border"] = coordinate_system["x_labels_stop"]
        prepared_data["msa_left_border"] = coordinate_system["x_msa_start"]
        prepared_data["sequence"] = record.seq
        prepared_data["label"] = record.description
        prepared_data["palette"] = additional_data["palette"]
        self.prepared_data = prepared_data

        return prepared_data

    def create_track(self) -> SequenceVis:
        """Initialise a Sequence track object.

        Returns:
            SequenceVis: visualisation track.

        """
        return SequenceVis(self.prepared_data, self.parameters)


class Image:
    """An Image object holds pdf.

    Attributes:
        canvas (reportlab.pdfgen.canvas.Canvas): pdf object of the reportlab library.

    """

    def __init__(self, filename: str, width: float, height: float):
        """Create an Image object.

        Arguments:
            filename (str): path and name of a pdf.
            width (float): width of a pdf.
            height (float): height of a pdf.

        """
        self.canvas = reportlab.pdfgen.canvas.Canvas(filename, pagesize=(width, height))

    def save(self) -> None:
        """Save a pdf file.

        Returns:
            None

        """
        self.canvas.save()
        return None
