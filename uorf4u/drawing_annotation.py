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


class AnnotationPlotManager:
    """
    AnnotationPlotManager object holds needed information for annotation visualisation and controls it.

    Note:
        It's supposed that the AnnotationPlotManager' objects will not be used directly by API users since
            visualisation can be controlled by 'plot_annotation' method.

    Attributes:
        path (uorf4u.data_processing.Path): Path class' objects that holds list of conserved ORFs.
        upstream_sequences (list): list of dicts with information about upstream sequences.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        coordinate_system (dict): coordinate system of figure.
        additional_data (dict): dict with data for visualisation tracks.

    """

    def __init__(self, path, upstream_sequences: list, parameters: uorf4u.manager.Parameters):
        """Create a AnnotationPlotManager object.

        Arguments:
            path (uorf4u.data_processing.Path): Path class' objects that holds list of conserved ORFs.
            upstream_sequences (list): list of dicts with information about upstream sequences.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        self.path = path
        self.upstream_sequences = upstream_sequences
        self.parameters = parameters
        self.coordinate_system = dict()
        self.additional_data = dict()

    def define_x_axis_coordinate_system(self) -> None:
        """Define coordinate system.

        Returns:
            None

        """

        label_height = self.parameters.arguments["label_height_to_orf_height"] * self.parameters.arguments[
            "orf_height"] * cm
        label_font_size = uorf4u.methods.string_height_to_font_size(label_height, "regular", self.parameters.arguments)
        self.additional_data["label_font_size"] = label_font_size
        self.additional_data["ordered_upstream_sequences"] = [self.upstream_sequences[i] for i in
                                                              [orf.useq_index for orf in self.path.path]]
        max_label_width = max([reportlab.pdfbase.pdfmetrics.stringWidth(i.annotations["label"], "regular",
                                                                        label_font_size)
                               for i in self.additional_data["ordered_upstream_sequences"]])
        self.additional_data["number_of_sequences"] = len(self.path)
        self.additional_data["max_upstream_sequence_length"] = max(
            i.annotations["upstream_region_length"] for i in self.additional_data["ordered_upstream_sequences"])
        self.additional_data["max_downstream_sequence_length"] = max(
            i.annotations["downstream_region_length"] for i in self.additional_data["ordered_upstream_sequences"])
        window_size_nt = self.additional_data["max_upstream_sequence_length"] + self.additional_data[
            "max_downstream_sequence_length"]
        if self.parameters.arguments["annotation_width"] == "auto":
            annotation_width = window_size_nt * self.parameters.arguments["mm_per_nt"] * mm
        else:
            annotation_width = self.parameters.arguments["annotation_width"] * cm
        self.coordinate_system["transformation_coef"] = annotation_width / window_size_nt
        self.coordinate_system["x_labels_start"] = self.parameters.arguments["margin"] * cm
        self.coordinate_system["x_labels_stop"] = self.coordinate_system["x_labels_start"] + max_label_width
        self.coordinate_system["x_annotation_start"] = self.coordinate_system["x_labels_stop"] + \
                                                       self.parameters.arguments["label_gap"] * cm
        self.coordinate_system["x_annotation_stop"] = self.coordinate_system["x_annotation_start"] + annotation_width
        self.coordinate_system["figure_width"] = 2 * self.parameters.arguments["margin"] * cm + annotation_width + \
                                                 max_label_width + self.parameters.arguments["label_gap"] * cm
        self.coordinate_system["figure_height"] = self.parameters.arguments["margin"] * cm

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
        for index in range(self.additional_data["number_of_sequences"]):
            upstream_sequence = self.additional_data["ordered_upstream_sequences"][index]
            conserved_orf = self.path.path[index]
            sequence_loader = SequencesLoader(self.parameters)
            sequence_loader.prepare_data(upstream_sequence, conserved_orf, self.coordinate_system, self.additional_data)
            track = sequence_loader.create_track()
            self.tracks.append(track)
            self.coordinate_system["figure_height"] += track.needed_y_space()
            self.coordinate_system["figure_height"] += self.parameters.arguments["gap"] * cm
            # if index < self.additional_data["number_of_sequences"] - 1:
        axis_tics_loader = AxisLoader(self.parameters)
        axis_tics_loader.prepare_data(self.coordinate_system, self.additional_data)
        axis_tics_track = axis_tics_loader.create_track()
        self.tracks.append(axis_tics_track)
        self.coordinate_system["figure_height"] += axis_tics_track.needed_y_space()
        self.coordinate_system["figure_height"] += self.parameters.arguments["margin"] * cm

    def plot(self, filename):
        image = Image(filename, self.coordinate_system["figure_width"], self.coordinate_system["figure_height"])
        current_y_top = self.coordinate_system["figure_height"] - self.parameters.arguments["margin"] * cm
        for track in self.tracks:
            track.visualisation_data["y_top"] = current_y_top
            track.draw(image.canvas)
            current_y_top -= (track.needed_space + self.parameters.arguments["gap"] * cm)
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


class TitleVis(Track):
    """Title visualisation track object draws figure's title.

    Note:
        This track currently is not supported.

    Attributes:
        visualisation_data (dict): a dictionary with data needed for visualisation.
        parameters (uorf4u.manager.Parameters): Parameters' class object.

    """

    def __init__(self, visualisation_data: dict, parameters: uorf4u.manager.Parameters):
        """Create TitleVis object.

        Arguments:
            visualisation_data (dict): a dictionary with data needed for visualisation.
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        self.visualisation_data = visualisation_data
        self.parameters = parameters

    def needed_y_space(self) -> float:
        """Calculate needed vertical space for a Title track.

        Returns:
            float: needed vertical space.

        """
        font_type = self.parameters.arguments["title_font_type"]
        reportlab.pdfbase.pdfmetrics.registerFont(
            reportlab.pdfbase.ttfonts.TTFont(font_type, self.parameters.arguments[f"font_{font_type}"]))
        face = reportlab.pdfbase.pdfmetrics.getFont(font_type).face
        if self.parameters.arguments["title_font_size"] == "auto":
            text_height = self.parameters.arguments["orf_height"] * cm
            font_size = uorf4u.methods.string_height_to_font_size(text_height, font_type, self.parameters.arguments)
            self.parameters.arguments["title_font_size"] = font_size
        else:
            text_height = (self.parameters.arguments["title_font_size"] * (face.ascent - face.descent)) / (
                    1000 * 1.38)
        self.visualisation_data["text_height"] = text_height
        self.needed_space = text_height * 1.2
        return self.needed_space

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a Title track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None

        """
        x_left_border = self.visualisation_data["coordinate_system"]["x_labels_start"]
        # x_left_border = self.visualisation_data["coordinate_system"]["x_annotation_start"]
        canvas.setFillColorRGB(*uorf4u.methods.get_color("label_color", self.parameters.arguments))
        canvas.setFont(self.parameters.arguments["title_font_type"], self.parameters.arguments["title_font_size"])
        canvas.drawString(x_left_border, self.visualisation_data["y_top"] - self.visualisation_data["text_height"],
                          self.visualisation_data["title"])


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
        self.needed_space = self.parameters.arguments["orf_height"] * cm
        return self.needed_space

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw a Sequence track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None

        """
        orf_height = self.parameters.arguments["orf_height"] * cm
        y_c = self.visualisation_data["y_top"] - 0.5 * orf_height
        x_offset = 0.5 * self.parameters.arguments["upstream_seq_line_width"]
        canvas.setStrokeColorRGB(*uorf4u.methods.get_color("upstream_seq_line_color", self.parameters.arguments))
        canvas.setLineCap(0)
        canvas.setLineWidth(self.parameters.arguments["upstream_seq_line_width"])
        canvas.line(self.visualisation_data["upstream_sequence_line_start_x"], y_c,
                    self.visualisation_data["upstream_sequence_line_stop_x"] - x_offset, y_c)
        # Cleaning the space:
        canvas.setFillColorRGB
        canvas.setStrokeColorRGB(1, 1, 1, 1)

        for orf_dict in self.visualisation_data["orfs_coordinates_dict"].values():

            canvas.setLineWidth(self.parameters.arguments["upstream_seq_line_width"]*1.5)
            canvas.line(orf_dict["x_start"], y_c, orf_dict["x_stop"], y_c)
            #canvas.rect(orf_dict["x_start"], y_c - orf_height / 2, orf_dict["x_stop"] - orf_dict["x_start"], orf_height,
            #            stroke=0, fill=1)
        if self.parameters.arguments["check_assembly_annotation"] and \
                "fasta" not in self.parameters.cmd_arguments.keys():
            for protein_id, cds_dict in self.visualisation_data["CDSs_coordinates_dict"].items():
                canvas.line(cds_dict["x_start"], y_c, cds_dict["x_stop"], y_c)
        canvas.setLineWidth(self.parameters.arguments["upstream_seq_line_width"])

        # Labels
        canvas.setFillColorRGB(*uorf4u.methods.get_color("label_color", self.parameters.arguments))
        canvas.setFont("regular", self.visualisation_data["label_font_size"])
        y_l = y_c - 0.5 * (self.parameters.arguments["label_height_to_orf_height"] * orf_height)
        canvas.drawRightString(self.visualisation_data["coordinate_system"]["x_labels_stop"], y_l,
                               self.visualisation_data["useq_label"])

        # main_CDS
        canvas.setLineWidth(self.parameters.arguments["orf_line_width"])
        canvas.setStrokeColorRGB(*uorf4u.methods.get_color("cds_seq_stroke_color", self.parameters.arguments))
        canvas.setFillColorRGB(*uorf4u.methods.get_color("cds_seq_fill_color", self.parameters.arguments))
        p = canvas.beginPath()
        p.moveTo(self.visualisation_data["main_CDS_stop_x"], y_c - orf_height / 2)
        p.lineTo(self.visualisation_data["main_CDS_start_x"], y_c - orf_height / 2)
        p.lineTo(self.visualisation_data["main_CDS_start_x"], y_c + orf_height / 2)
        p.lineTo(self.visualisation_data["main_CDS_stop_x"], y_c + orf_height / 2)
        canvas.drawPath(p, stroke=1, fill=1)
        # Other ORFs:
        for orf in self.visualisation_data["annotated_orfs"]:
            orf_dict = self.visualisation_data["orfs_coordinates_dict"][orf]
            if orf != self.visualisation_data["conserved_orf"]:
                fill_color = None
                stroke_color = uorf4u.methods.get_color("other_uorfs_stroke_color", self.parameters.arguments)

            else:
                fill_color = uorf4u.methods.get_color("conserved_uorfs_fill_color", self.parameters.arguments)
                stroke_color = uorf4u.methods.get_color("conserved_uorfs_stroke_color", self.parameters.arguments)
            self.orf_object(canvas, orf_dict["x_start"], orf_dict["x_stop"], y_c, orf_dict["strand"], orf_height,
                            orf_dict["left_out"], orf_dict["right_out"], fill_color, stroke_color)

        # Annotated in RefSeq CDSs
        if self.parameters.arguments["check_assembly_annotation"] and \
                "fasta" not in self.parameters.cmd_arguments.keys():
            fill_color = None
            stroke_color = uorf4u.methods.get_color("annotated_orf_stroke_color", self.parameters.arguments)
            for protein_id, cds_dict in self.visualisation_data["CDSs_coordinates_dict"].items():
                self.orf_object(canvas, cds_dict["x_start"], cds_dict["x_stop"], y_c, cds_dict["strand"], orf_height,
                                cds_dict["left_out"], cds_dict["right_out"], fill_color, stroke_color)
        return None

    def orf_object(self, canvas: reportlab.pdfgen.canvas.Canvas, x_start: float, x_stop: float, y_c: float, strand: str,
                   height: float, left_out: bool, right_out: bool, fill_color: str, stroke_color: str) -> None:
        """Method for drawing an ORF's polygon.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.
            x_start (float): ORF's start coordinate (already transformed to pdf's)
            x_stop (float): ORF's stop coordinate (already transformed to pdf's)
            y_c: (float): centred y coordinate of a current track.
            strand (str): strand of an ORF.
            height (float): height of a polygon.
            left_out (bool): whether an ORF is out of range on the left.
            right_out (bool): whether an ORF is out of range on the right.
            fill_color (str): fill color of a polygon.
            stroke_color (str): stroke color of a polygon.

        Returns:
            None
        """
        fill, stroke = 0, 0
        if stroke_color:
            canvas.setStrokeColorRGB(*stroke_color)
            stroke = 1
        if fill_color:
            canvas.setFillColorRGB(*fill_color)
            fill = 1
        arrow_length = min(height, (x_stop - x_start))
        p = canvas.beginPath()
        if strand == "+" and not left_out and not right_out:
            p.moveTo(x_start, y_c)
            p.lineTo(x_start, y_c + height / 2)
            p.lineTo(x_stop - arrow_length, y_c + height / 2)
            p.lineTo(x_stop, y_c)
            p.lineTo(x_stop - arrow_length, y_c - height / 2)
            p.lineTo(x_start, y_c - height / 2)
            p.lineTo(x_start, y_c)
        elif strand == "+" and left_out and not right_out:
            p.moveTo(x_start, y_c + height / 2)
            p.lineTo(x_stop - arrow_length, y_c + height / 2)
            p.lineTo(x_stop, y_c)
            p.lineTo(x_stop - arrow_length, y_c - height / 2)
            p.lineTo(x_start, y_c - height / 2)
        elif strand == "+" and right_out and not left_out:
            p.moveTo(x_stop, y_c + height / 2)
            p.lineTo(x_start, y_c + height / 2)
            p.lineTo(x_start, y_c - height / 2)
            p.lineTo(x_stop, y_c - height / 2)
        elif strand == "-" and not left_out and not right_out:
            p.moveTo(x_stop, y_c)
            p.lineTo(x_stop, y_c + height / 2)
            p.lineTo(x_start + arrow_length, y_c + height / 2)
            p.lineTo(x_start, y_c)
            p.lineTo(x_start + arrow_length, y_c - height / 2)
            p.lineTo(x_stop, y_c - height / 2)
            p.lineTo(x_stop, y_c)
        elif strand == "-" and right_out and not left_out:
            p.moveTo(x_stop, y_c + height / 2)
            p.lineTo(x_start + arrow_length, y_c + height / 2)
            p.lineTo(x_start, y_c)
            p.lineTo(x_start + arrow_length, y_c - height / 2)
            p.lineTo(x_stop, y_c - height / 2)
        elif strand == "-" and left_out and not right_out:
            p.moveTo(x_start, y_c + height / 2)
            p.lineTo(x_stop, y_c + height / 2)
            p.lineTo(x_stop, y_c - height / 2)
            p.lineTo(x_start, y_c - height / 2)
        canvas.drawPath(p, stroke=stroke, fill=fill)


class TicsVis(Track):
    """TicsVis track draws axis tics.

    Attributes:
        visualisation_data (dict): a dictionary with data needed for visualisation.
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        needed_space (float): needed vertical space for a track.

    """

    def __init__(self, visualisation_data: dict, parameters: uorf4u.manager.Parameters):
        """Create a TicsVis object.

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
        font_type = "regular"
        reportlab.pdfbase.pdfmetrics.registerFont(
            reportlab.pdfbase.ttfonts.TTFont(font_type, self.parameters.arguments[f"font_{font_type}"]))
        face = reportlab.pdfbase.pdfmetrics.getFont("regular").face
        if self.parameters.arguments["axis_tics_font_size"] == "auto":
            text_height = self.parameters.arguments["label_height_to_orf_height"] * self.parameters.arguments[
                "orf_height"] * cm
            font_size = uorf4u.methods.string_height_to_font_size(text_height, "regular", self.parameters.arguments)
            self.parameters.arguments["axis_tics_font_size"] = font_size
        else:
            text_height = (self.parameters.arguments["axis_tics_font_size"] * (face.ascent - face.descent)) / (
                    1000 * 1.38)
        self.visualisation_data["tics_height"] = 0.7 * text_height
        self.visualisation_data["text_space"] = 1.2 * text_height
        self.needed_space = self.visualisation_data["tics_height"] + self.visualisation_data["text_space"]
        return self.needed_space

    def draw(self, canvas: reportlab.pdfgen.canvas.Canvas) -> None:
        """Draw an AxisTics track.

        Arguments:
            canvas (reportlab.pdfgen.canvas.Canvas): a pdf object.

        Returns:
            None

        """
        y_top = self.visualisation_data["y_top"]
        canvas.setLineCap(2)
        canvas.setLineWidth(self.parameters.arguments["axis_tics_line_width"])
        canvas.setStrokeColorRGB(*uorf4u.methods.get_color("label_color", self.parameters.arguments))
        canvas.setFillColorRGB(*uorf4u.methods.get_color("label_color", self.parameters.arguments))
        canvas.setFont("regular", self.parameters.arguments["axis_tics_font_size"])
        canvas.line(self.visualisation_data["coordinate_system"]["x_annotation_start"], y_top,
                    self.visualisation_data["coordinate_system"]["x_annotation_stop"], y_top)
        for tic_label, tic_position in self.visualisation_data["tics"].items():
            canvas.line(tic_position, y_top, tic_position, y_top - self.visualisation_data["tics_height"])
            if tic_label == -self.visualisation_data["max_upstream_sequence_length"]:
                canvas.drawString(tic_position,
                                  y_top - (self.visualisation_data["tics_height"] + self.visualisation_data[
                                      "text_space"]), str(tic_label))
            elif tic_label == self.visualisation_data["max_downstream_sequence_length"]:
                canvas.drawRightString(tic_position,
                                       y_top - (self.visualisation_data["tics_height"] + self.visualisation_data[
                                           "text_space"]), str(tic_label))
            else:
                canvas.drawCentredString(tic_position,
                                         y_top - (self.visualisation_data["tics_height"] + self.visualisation_data[
                                             "text_space"]), str(tic_label))


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


class TitleLoader(Loader):
    """A TitleLoader object prepares data for a Title track object.

    Note:
        Title track currently is not available.

    Attributes:
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        prepared_data (dict): dict with data needed for a visualisation track.

    """

    def __init__(self, parameters):
        """Create a TitleLoader object.

        Arguments:
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        super().__init__(parameters)

    def prepare_data(self, coordinate_system: dict, additional_data: dict) -> dict:
        """Prepare data for Title visualisation track.

        Attributes:
            coordinate_system (dict): coordinate system of a figure page.
            additional_data (dict): data needed for a track initialisation.

        Returns:
            dict: dictionary with prepared data for visualisation.

        """
        prepared_data = dict()
        prepared_data["title"] = "Title Testing"
        prepared_data["coordinate_system"] = coordinate_system
        self.prepared_data = prepared_data
        return prepared_data

    def create_track(self) -> TitleVis:
        """Initialise a Title track object.

        Returns:
            TitleVis: visualisation track.

        """
        return TitleVis(self.prepared_data, self.parameters)


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

    def prepare_data(self, upstream_sequence: Bio.SeqRecord.SeqRecord, conserved_orf, coordinate_system: dict,
                     additional_data: dict) -> dict:
        """Prepare data for a Title visualisation track.

        Attributes:
            upstream_sequence (dict): upstream sequence' data.
            conserved_orf (uorf4u.data_processing.ORF): conserved ORF on the upstream sequence.
            coordinate_system (dict): coordinate system of a figure page.
            additional_data (dict): data needed for a track initialisation.

        Returns:
            dict: dictionary with prepared data for visualisation.

        """
        prepared_data = dict()
        max_upstream_sequence_length = additional_data["max_upstream_sequence_length"]
        prepared_data["coordinate_system"] = coordinate_system
        prepared_data["label_font_size"] = additional_data["label_font_size"]
        prepared_data["label_right_border"] = coordinate_system["x_labels_stop"]
        prepared_data["upstream_sequence_line_start_x"] = coordinate_system["x_annotation_start"] + \
                                                          ((max_upstream_sequence_length -
                                                            upstream_sequence.annotations["upstream_region_length"]) * \
                                                           coordinate_system["transformation_coef"])
        prepared_data["upstream_sequence_line_stop_x"] = coordinate_system["x_annotation_start"] + \
                                                         (max_upstream_sequence_length *
                                                          coordinate_system["transformation_coef"])
        prepared_data["main_CDS_start_x"] = coordinate_system["x_annotation_start"] + \
                                            (max_upstream_sequence_length *
                                             coordinate_system["transformation_coef"])
        prepared_data["main_CDS_stop_x"] = coordinate_system["x_annotation_start"] + \
                                           ((max_upstream_sequence_length +
                                             upstream_sequence.annotations["downstream_region_length"]) *
                                            coordinate_system["transformation_coef"])

        prepared_data["orfs_coordinates_dict"] = {k: v for k, v in
                                                  zip(upstream_sequence.annotations["ORFs"],
                                                      [self.calculate_orf_position(i.start, i.stop, "+",
                                                                                   upstream_sequence,
                                                                                   max_upstream_sequence_length,
                                                                                   coordinate_system) for i in
                                                       upstream_sequence.annotations["ORFs"]])}
        prepared_data["useq_label"] = upstream_sequence.annotations["label"]
        prepared_data["annotated_orfs"] = [orf for orf in upstream_sequence.annotations["ORFs"] if orf != conserved_orf]
        prepared_data["annotated_orfs"].append(conserved_orf)
        prepared_data["conserved_orf"] = conserved_orf
        if self.parameters.arguments["check_assembly_annotation"] and upstream_sequence.annotations["RefSeq"]:
            prepared_data["CDSs"] = [i for i in upstream_sequence.annotations["locus_annotation"].CDSs if
                                     i["relative_start"] != upstream_sequence.annotations["upstream_region_length"]]
            prepared_data["CDSs_coordinates_dict"] = {k: v for k, v in
                                                      zip([i["protein_id"] for i in
                                                           prepared_data["CDSs"]],
                                                          [self.calculate_orf_position(i["relative_start"],
                                                                                       i["relative_stop"],
                                                                                       i["relative_strand"],
                                                                                       upstream_sequence,
                                                                                       max_upstream_sequence_length,
                                                                                       coordinate_system) for i in
                                                           prepared_data["CDSs"]])}
        else:
            prepared_data["CDSs"] = None
        self.prepared_data = prepared_data

        return prepared_data

    def calculate_orf_position(self, start: int, stop: int, strand: str, useq: Bio.SeqRecord.SeqRecord,
                               max_upstream_sequence_length: int, coordinate_system: dict) -> dict:
        """Transform an ORF's nucleotide coordinates to pdf's coordinates.

        Arguments:
            start (int): start coordinate in nt.
            stop (int): stop coordinate in nt.
            strand (str): strand of an ORF.
            useq (dict): current upstream sequence.
            max_upstream_sequence_length (int): max length of upstream sequences for visualisation.
            coordinate_system (dict): coordinate system of a figure.

        Returns:
            dict: transformed orf's coordinates.
        """
        orf_coordinates = dict()
        orf_coordinates["x_start"] = coordinate_system["x_annotation_start"] + (
                max(0, start) + (max_upstream_sequence_length - useq.annotations["upstream_region_length"])) * \
                                     coordinate_system["transformation_coef"]
        orf_coordinates["x_stop"] = coordinate_system["x_annotation_start"] + (
                min(stop, useq.annotations["length"]) + (
                max_upstream_sequence_length - useq.annotations["upstream_region_length"])) * \
                                    coordinate_system["transformation_coef"]
        orf_coordinates["strand"] = strand
        orf_coordinates["left_out"] = start < 0
        orf_coordinates["right_out"] = stop > useq.annotations["length"]
        return orf_coordinates

    def create_track(self) -> SequenceVis:
        """Initialise a Sequence track object.

        Returns:
            SequenceVis: visualisation track.

        """
        return SequenceVis(self.prepared_data, self.parameters)


class AxisLoader(Loader):
    """An AxisLoader object prepares data for an Axis track object.


    Attributes:
        parameters (uorf4u.manager.Parameters): Parameters' class object.
        prepared_data (dict): dict with data needed for a visualisation track.

    """

    def __init__(self, parameters):
        """Create an AxisLoader object.

        Arguments:
            parameters (uorf4u.manager.Parameters): Parameters' class object.

        """
        super().__init__(parameters)

    def prepare_data(self, coordinate_system: dict, additional_data: dict):
        """Prepare data for an Axis visualisation track.

        Attributes:
            coordinate_system (dict): coordinate system of a figure page.
            additional_data (dict): data needed for a track initialisation.

        Returns:
            dict: dictionary with prepared data for visualisation.

        """
        prepared_data = dict()
        prepared_data["coordinate_system"] = coordinate_system
        prepared_data["max_upstream_sequence_length"] = additional_data["max_upstream_sequence_length"]
        prepared_data["max_downstream_sequence_length"] = additional_data["max_downstream_sequence_length"]
        step = int(round(additional_data["max_upstream_sequence_length"] / 2, -2))
        tics = [-additional_data["max_upstream_sequence_length"], 0, additional_data["max_downstream_sequence_length"]]
        x_tic_centred = int(round(-additional_data["max_upstream_sequence_length"] / 2, -2))
        tics.append(x_tic_centred)
        x_tic_left, x_tic_right = x_tic_centred - step, x_tic_centred + step
        while x_tic_right < 0 and x_tic_left > - additional_data["max_upstream_sequence_length"]:
            tics.append(x_tic_left)
            tics.append(x_tic_right)
            x_tic_left -= step
            x_tic_right += step
        tics.sort()
        tics_coordinates = [self.transform_relative_position_to_x_coordinate(i, coordinate_system, additional_data[
            "max_upstream_sequence_length"]) for i in tics]
        prepared_data["tics"] = {k: v for k, v in zip(tics, tics_coordinates)}
        self.prepared_data = prepared_data
        return prepared_data

    def transform_relative_position_to_x_coordinate(self, relative_position: int, coordinate_system: dict,
                                                    max_upstream_sequence_length: int) -> float:
        """Transform nucleotide x coordinate to pdf's.

        Arguments:
            relative_position (int): nucleotide position
            coordinate_system (dict): coordinate system of a figure.
            max_upstream_sequence_length (int): max length of upstream sequences for visualisation.

        Returns:
            float: transformed x coordinate.
        """
        return coordinate_system["x_annotation_start"] + (relative_position + max_upstream_sequence_length) * \
               coordinate_system["transformation_coef"]

    def create_track(self) -> TicsVis:
        """Initialise a Tics track object.

        Returns:
            TicsVis: visualisation track.

        """
        return TicsVis(self.prepared_data, self.parameters)


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
