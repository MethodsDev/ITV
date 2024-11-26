import logging
import os
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
import re
import math
import copy
from intervaltree import Interval, IntervalTree
from ipywidgets.embed import embed_minimal_html, dependency_state
import cairosvg

RESVG = None


def split_svg_header_to_matrix(svg_header, x_width, y_height):
    # Extract width and height using regex
    width_match = re.search(r'width="(\d+)"', svg_header)
    height_match = re.search(r'height="(\d+)"', svg_header)

    if not width_match or not height_match:
        raise ValueError("Invalid SVG header: width or height not found.")

    original_width = int(width_match.group(1))
    original_height = int(height_match.group(1))

    # Split into a 2D matrix of headers
    headers_matrix = []
    x_start = 0
    while x_start < original_width:
        col = []
        y_start = 0
        while y_start < original_height:
            slice_width = min(x_width, original_width - x_start)
            slice_height = min(y_height, original_height - y_start)
            view_box = f'viewBox="{x_start} {y_start} {slice_width} {slice_height}"'
            new_header = re.sub(
                r'(width=".*?")',
                f'width="{slice_width}"',
                re.sub(
                    r'(height=".*?")',
                    f'height="{slice_height}" {view_box}',
                    svg_header
                )
            )
            col.append(new_header)
            y_start += y_height
        headers_matrix.append(col)
        x_start += x_width

    return headers_matrix


def save(doc, output_path, output_format=None, requested_converter=None):
    """
    Saves document `doc` to a file at `output_path`. By default, this file 
    will be in SVG format; if it ends with .pdf or .png, or if output_format
    is specified, the document will be converted to PDF or PNG if possible.

    Conversion to PDF and PNG require rsvg-convert (provided by librsvg), 
    inkscape or webkitToPDF (PDF conversion only).

    Attributes:
        doc: the :py:class:`integrative_transcriptomics_viewer.Document` to be saved
        output_path: a string specifying the file to save to; file extensions of
            .pdf or .png will change the default output format
        output_format: override the file format; must be one of "pdf", "png", or
            (the default) "svg"
    """
    if isinstance(output_path, bytes):
        output_path = output_path.decode()

    if output_format is None:
        if output_path.lower().endswith(".pdf"):
            output_format = "pdf"
        elif output_path.lower().endswith(".png"):
            output_format = "png"
        elif output_path.lower().endswith(".svg"):
            output_format = "svg"
        else:
            if not output_path.lower().endswith(".html"):
                output_path += ".html"
            output_format = "html"

    if output_format == "html":
        with open(output_path, "w") as f:
            f.write(doc._repr_html_())

    if output_format == "svg":
        with open(output_path, "w") as outf:
            render_to_file(doc, outf)

    elif output_format == "png":
        if output_path.lower().endswith(".png"):
            output_path_prefix = output_path[:-4]
        else:
            output_path_prefix = output_path

        if requested_converter == "cairosvg" or requested_converter is None:  # add a check that cairosvg available
            requested_converter = "cairosvg"
        elif (requested_converter is None or requested_converter == "resvg") and _checkRESVGConvert():
            requested_converter = "resvg"
        elif (requested_converter is None or requested_converter == "librsvg"):
            requested_converter = "librsvg"
        else:
            raise RuntimeError("No SVG to PNG converter provided and available.")

        full_svg = doc._repr_svg__()
        svg_header, _, svg_body = full_svg.partition("\n")
        headers_matrix = split_svg_header_to_matrix(svg_header, 10000, 10000)

        for x_idx, col in enumerate(headers_matrix):
            for y_idx, header in enumerate(col):
                filename = f"_x{x_idx:02}_y{y_idx:02}.png"  # Zero-pad the number, adjust i+1 if numbering should start from 01

                if requested_converter == "resvg":
                    split_png = _convertSVG_resvg_stdio(header + "\n" + svg_body)
                    with open(output_path_prefix + filename, 'wb') as outf:
                        outf.write(split_png)
                else:  # methods that require to save to file for convertion
                    with tempfile.TemporaryDirectory() as outdir:
                        temp_svg_path = os.path.join(outdir, "temp.svg")
                        with open(temp_svg_path, "w") as outf:
                           outf.write(header + "\n" + svg_body)

                        if requested_converter == "librsvg":
                            convert_svg(temp_svg_path, output_path_prefix + filename, output_format, requested_converter="librsvg")
                        elif requested_converter == "cairosvg":
                            cairosvg.svg2png(url=temp_svg_path, write_to=output_path_prefix + filename)

    else: # pdf
        # do something
        return


def save_html_view(view, output_path, title=""):
    embed_minimal_html(output_path, views=view, state=dependency_state(view), title=title)



def render_to_file(doc, output_handle):
    """
    Renders the document as an svg to a file-like object.
    """

    for l in doc.render():
        output_handle.write(l + "\n")




#############################################################################
######################### low-level functionality ###########################
#############################################################################

class SvgSplitter:
    # section_regex = re.compile(".*clipPath$")
    # read_regex = re.compile(".*path$")
    # text_regex = re.compile(".*text$")
    # path_dot = re.compile("L [0-9\.]+ ([0-9\.]+)")
    path_y_regex = re.compile(r'(?:M|L) (?:[-0-9\.]+) ([-0-9\.]+)')
    rect_y_regex = re.compile(r'y="([-0-9\.]+)" .* height="([-0-9\.]+)')
    line_y_regex = re.compile(r'y[12]="([-0-9\.]+)')
    text_y_regex = re.compile(r'y="([-0-9\.]+)')
    
    def __init__(self, root):
        ET.register_namespace("", "http://www.w3.org/2000/svg")  # Register the SVG namespace the avoid "ns0:" or ":ns0" being inserted in tags
        self.root = root
        self.split_template = self.make_template_svg_split(root)
        self.split_svgs = []

    def make_template_svg_split(self, root):
        return ET.Element('svg', {
                'version': root.attrib['version'],
                'width': root.attrib['width'],
                'height': "0",
                'shape-rendering': "crispEdges"
                # 'height': str(min(max_height, original_height - i * max_height)) #,
                # 'xmlns': root.attrib['xmlns']
            })
    
    def get_new_svg_split(self):
        return copy.deepcopy(self.split_template)

    # def update_splits_height(self, element, current_height = 0):
    #     for subelement in element:
    #         current_height += update_splits_height(subelement, current_height)
    #     current_height += element.attrib['height']

    def update_splits_height(self, element, current_height=0):
        for subelement in element:
            current_height += self.update_splits_height(subelement, current_height)
        current_height += float(element.attrib.get('height', 0))
        return current_height
    
    def write_splits(self, prefix):
        output_files = []
        for i, split_svg in enumerate(self.split_svgs):
            split_svg.write(f'{prefix}_slice_{i + 1}.svg')
            output_files.append(f'{prefix}_slice_{i + 1}.svg')
        return(output_files)

    def get_splits(self):
        return self.split_svgs

    def is_read_block(self, element):
        if self.has_children(element):
            return False
        path_found = False
        for subelement in element:
            if self.is_path(subelement) and len(subelement.attrib['d'].split(" ")) == 16:
                return True
        return False

    def is_clipPath(self, element):
        return True if element.tag[-8:] == "clipPath" else False

    def is_path(self, element):
        return True if element.tag[-4:] == "path" else False

    def is_rect(self, element):
        return True if element.tag[-4:] == "rect" else False

    def is_text(self, element):
        return True if element.tag[-4:] == "text" else False

    def is_line(self, element):
        return True if element.tag[-4:] == "line" else False
    
    def is_defs(self, element):
        return True if element.tag[-4:] == "defs" else False

    def is_g(self, element):
        return True if element.tag[-1] == "g" else False
 
    def has_children(self, element):
        return element.find("*")
        # for el in element.iter():
        #     if el == element:
        #         continue
        #     if self.is_clipPath(el):  # or self.is_g(el):
        #         return True
        # return False
   
    def get_path_max_y(self, element):
        max_y = 0

        matches = self.path_y_regex.finditer(element.attrib['d'])
        for match in matches:
            if float(match.group(1)) > max_y:
                max_y = float(match.group(1))

        return max_y

    def get_path_y_range(self, element):
        min_y = math.inf
        max_y = -math.inf

        matches = self.path_y_regex.finditer(element.attrib['d'])
        for match in matches:
            if float(match.group(1)) > max_y:
                max_y = float(match.group(1))
            if float(match.group(1)) < min_y:
                min_y = float(match.group(1))

        return (min_y, max_y)  # could return an Interval(min_y, max_y)
        
    # def get_line_y_range(self, element):
    #     min_y = -math.inf
    #     max_y = math.inf

    #     matches = line_y_regex.finditer(element.tag)
    #     for match in matches:
    #         if float(match.group(1)) > max_y:
    #             max_y = float(match.group(1))
    #         if float(match.group(1)) < min_y:
    #             min_y = float(match.group(1))
        
    #     return (min_y, max_y)  # could return an Interval(min_y, max_y)

    def get_line_y_range(self, element):
        min_y = min(float(element.attrib['y1']), float(element.attrib['y2']))
        max_y = max(float(element.attrib['y1']), float(element.attrib['y2']))
        return (min_y, max_y)  # could return an Interval(min_y, max_y)

    def get_line_y(self, element):  # assumes only horizontal lines
        return float(element.attrib['y1'])    
    
    # def get_rect_y_range(self, element):
    #     match = rect_regex.search(element.tag)
    #     return (float(matfch.group(1)), float(match.group(1)) + float(match.group(2))) # could return an Interval(min_y, max_y)

    def get_rect_y_range(self, element):
        return (float(element.attrib['y']), float(element.attrib['y']) + float(element.attrib['height']))
    
    def get_text_y_range(self, element):
        return (float(element.attrib['y']) - float(element.attrib['font-size']), float(element.attrib['y']))

    # can probably just keep the min/max y based on is_path elements alone
    def get_element_min_y(self, element):
        if self.is_path(element):
            return self.get_path_y_range(element)[0]
        elif self.is_rect(element):
            return self.get_rect_y_range(element)[0]
        elif self.is_text(element):
            return self.get_text_y_range(element)[0]
        return math.inf
        # elif self.is_line(element):
        #     return self.get_line_y(element)
        # elif self.is_defs(element):
        #    return None

    def get_element_max_y(self, element):
        if self.is_path(element):
            return self.get_path_y_range(element)[1]
        elif self.is_rect(element):
            return self.get_rect_y_range(element)[1]
        elif self.is_text(element):
            return self.get_text_y_range(element)[1]
        return -math.inf
        # elif self.is_line(element):
        #     return self.get_line_y(element)

    def get_element_min_y_rec(self, element):
        current_min = math.inf
        for sub_element in element:
            current_min = min(current_min, self.get_element_min_y_rec(sub_element))
        return min(current_min, self.get_element_min_y(element))


    def get_element_max_y_rec(self, element):
        current_max = -math.inf
        for sub_element in element:
            current_max = max(current_max, self.get_element_max_y_rec(sub_element))
        return max(current_max, self.get_element_max_y(element))


    def find_or_create_binterval(self, y_range, bins_tree, max_height):
        # Search for overlapping bins
        overlapping_bins = bins_tree.overlap(y_range[0], y_range[1])
        if len(overlapping_bins) == 1:
            return list(overlapping_bins)[0]
        elif len(overlapping_bins) == 2:
            overlapping_bins = list(sorted(overlapping_bins))
            overlapping_bins[0].data['max_y'] = y_range[1]
            bins_tree.remove(overlapping_bins[0])
            to_return = Interval(overlapping_bins[0].begin, y_range[1], overlapping_bins[0].data)
            bins_tree.add(to_return)
            overlapping_bins[1].data['min_y'] = y_range[0] + 0.01
            bins_tree.remove(overlapping_bins[1])
            bins_tree.add(Interval(y_range[0] + 0.01, overlapping_bins[1].end, overlapping_bins[1].data))
            return to_return

        if len(bins_tree) >= 1:
            current_binterval = sorted(bins_tree)[0]
            if y_range[0] < current_binterval.begin:
                while True:
                    previous_binterval = bins_tree.at(current_binterval.begin - 1)
                    if not previous_binterval:
                        new_bin = {'elements': [], 'min_y': current_binterval.begin - max_height, 'max_y': current_binterval.begin - 0.01}
                        new_binterval = Interval(current_binterval.begin - max_height, current_binterval.begin - 0.01, new_bin)
                        bins_tree.add(new_binterval)
                        current_binterval = new_binterval
                    else:
                        current_binterval = list(previous_binterval)[0]
                    if current_binterval.begin < y_range[0] < current_binterval.end:
                        return list(bins_tree.overlap(y_range[0], y_range[1]))[0]
            else:
                while True:
                    next_binterval = bins_tree.at(current_binterval.end + 1)
                    if not next_binterval:
                        new_bin = {'elements': [], 'min_y': current_binterval.end + 0.01, 'max_y': current_binterval.end + max_height}
                        new_binterval = Interval(current_binterval.end + 0.01, current_binterval.end + max_height, new_bin)
                        bins_tree.add(new_binterval)
                        current_binterval = new_binterval
                    else:
                        current_binterval = list(next_binterval)[0]
                    if current_binterval.begin <= y_range[0] < current_binterval.end:
                        return list(bins_tree.overlap(y_range[0], y_range[1]))[0]
    
        else:
            new_bin = {'elements': [], 'min_y': y_range[0], 'max_y': y_range[0] + max_height}
            new_binterval = Interval(y_range[0], y_range[0] + max_height - 0.01, new_bin)
            bins_tree.add(new_binterval)
            return new_binterval

    def create_bins(self, elements, max_height):
        bins_tree = IntervalTree()
        defs_elements = []
        current_bin = None
        initial_intron_line = None

        for element in elements:
            if self.is_defs(element):
                defs_elements.append(element)
                current_bin = None  # Reset current bin for elements following 'defs'
                continue

            # For 'path' elements
            if self.is_path(element):
                y_range = self.get_path_y_range(element)

                current_bin = self.find_or_create_binterval(y_range, bins_tree, max_height).data

                if initial_intron_line:
                    current_bin['elements'].append(initial_intron_line)
                    initial_intron_line = None

            # Add element to the current bin
            if current_bin:
                current_bin['elements'].append(element)
            else:
                # Handle 'text' elements following 'defs' elements
                if self.is_line(element):
                    initial_intron_line = element
                elif self.is_text(element):
                    y_range = self.get_text_y_range(element)
                    bin_for_text = self.find_or_create_binterval(y_range, bins_tree, max_height).data
                    bin_for_text['elements'].append(element)

        # Convert IntervalTree to list of bins
        bins = [interval.data for interval in sorted(bins_tree)]

        # Append 'defs' elements to each bin after creation
        for bin in bins:
            bin['elements'] = defs_elements + bin['elements']

        return bins
    
    def split_svg(self, current_element, current_split = None, max_height = 32000, current_height = 0, root = True):
        last_clip_path_height = 0
        
        for child_element in current_element:
            if current_split is None:
                current_split = self.get_new_svg_split()
            
            if child_element.tag.endswith('clipPath'):
                last_clip_path_height = float(child_element.find('.//{http://www.w3.org/2000/svg}rect').attrib['height'])
                continue
            elif self.is_g(child_element):
                child_element_height = last_clip_path_height
            else:
                child_element_height = float(child_element.attrib['height'])

            if float(current_split.attrib['height']) + child_element_height < max_height:
                # keep as whole with wtv is before, and can keep going
                current_split.append(child_element)
                current_split.attrib['height'] = str(float(current_split.attrib['height']) + child_element_height + 10) # +10 because there is a 10 between each group at the same level
            elif child_element_height < max_height:
                # finish previous part and start new one with this
                if self.has_children(current_split):
                    # need to add offsetting of everything inside like is done in the loop over bins
                    y_offset = self.get_element_min_y_rec(current_split)
                    y_trim = self.get_element_max_y_rec(current_split)

                    current_split.attrib['height'] = str(y_trim - y_offset)
                    # need to subtract min_y from all height/y coordinates
                    self.offset_y_rec(current_split, y_offset)
                    self.split_svgs.append(ET.ElementTree(current_split))

                current_split = self.get_new_svg_split()
                current_split.append(child_element)
                current_split.attrib['height'] = str(child_element_height + 10) # +10 because there is a 10 between each group at the same level

                y_offset = math.inf
                for elem in current_split.iter():
                    current_min_y = self.get_element_min_y(elem)
                    if current_min_y and current_min_y < y_offset:
                        y_offset = current_min_y
                for elem in current_split.iter():
                    self.offset_y(elem, y_offset)
                
            else:  # needs to be split somehow, adding however much can be added first, then make new splits with the rest
                # add this level branch
                if self.has_children(child_element):
                    # go deeper
                    current_split = self.split_svg(child_element, current_split, max_height, root = False)
                    # update the height of this, maybe do at the end of the whole recursion instead
                # below here, can probable close the split after this element is processed before going further to make subsections easier to handle
                elif self.is_read_block(child_element):
                    if self.has_children(current_split):
                        # need to add offsetting of everything inside like is done in the loop over bins
                        y_offset = self.get_element_min_y_rec(current_split)
                        y_trim = self.get_element_max_y_rec(current_split)

                        current_split.attrib['height'] = str(y_trim - y_offset)
                        # need to subtract min_y from all height/y coordinates
                        self.offset_y_rec(current_split, y_offset)
                        self.split_svgs.append(ET.ElementTree(current_split))

                    elements = list(child_element)
                    bins = self.create_bins(elements, max_height)

                    for bin in bins:
                        new_svg = self.get_new_svg_split()
                        new_svg.attrib['height'] = str(bin['max_y'] - bin['min_y'] + 0.01)
                        #if bin == bins[0]:
                        y_offset = math.inf
                        y_trim = -math.inf
                        for elem in bin['elements']:
                            y_offset = min(y_offset, self.get_element_min_y(elem))
                            y_trim = max(y_trim, self.get_element_max_y(elem))
                        new_svg.attrib['height'] = str(y_trim - y_offset)
                        for elem in bin['elements']:
                            # need to subtract min_y from all height/y coordinates
                            self.offset_y(elem, y_offset)
                            new_svg.append(elem)
                        self.split_svgs.append(ET.ElementTree(new_svg))
                    current_split = None
        if root and current_split and self.has_children(current_split):
            self.split_svgs.append(ET.ElementTree(current_split))
        return current_split

    def offset_y(self, element, offset):
        if self.is_path(element):
            self.offset_path_y(element, offset)
        elif self.is_line(element):
            self.offset_line_y(element, offset)
        elif self.is_rect(element):
            self.offset_rect_y(element, offset)
        elif self.is_text(element):
            self.offset_text_y(element, offset)
        # self.is_defs(element) doesn't need updating

    def offset_y_rec(self, element, offset):
        for sub_element in element:
            self.offset_y_rec(sub_element, offset)
        self.offset_y(element, offset)

    def offset_path_y(self, element, offset):
        groups = element.attrib['d'].split(" ")
        for i in range(2, len(groups), 3):
            if groups[i].replace('.', '', 1).isdigit():  # Check if it's a number
                groups[i] = str(float(groups[i]) - offset)
        element.attrib['d'] = ' '.join(groups)
        
    def offset_line_y(self, element, offset):  # assumes only horizontal lines
        element.attrib['y1'] = str(float(element.attrib['y1']) - offset)
        element.attrib['y2'] = str(float(element.attrib['y2']) - offset)

    def offset_rect_y(self, element, offset):
        element.attrib['y'] = str(float(element.attrib['y']) - offset)
    
    def offset_text_y(self, element, offset):
        element.attrib['y'] = str(float(element.attrib['y']) - offset)


def convert_svg(inpath, outpath, output_format, requested_converter=None):
    converter = _getExportConverter(output_format, requested_converter=requested_converter)

    if converter == "webkittopdf":
        exportData = _convertSVG_webkitToPDF(inpath, outpath, output_format)
    elif converter == "resvg":
        exportData = _convertSVG_resvg(inpath, outpath)
    elif converter == "librsvg":
        exportData = _convertSVG_rsvg_convert(inpath, outpath, output_format)
    elif converter == "inkscape":
        exportData = _convertSVG_inkscape(inpath, outpath, output_format)

    return exportData


def _getExportConverter(exportFormat, requested_converter=None):
    if requested_converter == "webkittopdf" and exportFormat=="png":
        logging.error("webkitToPDF does not support export to PNG; use librsvg or inkscape instead, or "
            "export to PDF")
        sys.exit(1)

    if exportFormat == "png" and requested_converter in [None, "resvg"]:
        if _checkRESVGConvert():
            return "resvg"

    if exportFormat == "png" and requested_converter in [None, "rsvg-convert", "librsvg"]:
        if _checkRSVGConvert():
            return "librsvg"

    if requested_converter in [None, "webkittopdf"]:
        if _checkWebkitToPDF():
            return "webkittopdf"

    if requested_converter in [None, "inkscape"]:
        if _checkInkscape():
            return "inkscape"

    raise Exception("No converter found for conversion to {}".format(exportFormat))
    return None



def _checkWebkitToPDF():
    try:
        subprocess.check_call("webkitToPDF", stderr=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def _checkRESVGConvert():
    global RESVG
    if RESVG is not None:
        return True

    try:
        cmd = [f"{sys.exec_prefix}/bin/resvg", "-V"]
        subprocess.check_call(cmd, stdout=subprocess.PIPE)
        RESVG = [f"{sys.exec_prefix}/bin/resvg", "--use-fonts-dir", f"{sys.exec_prefix}/fonts/"]
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        try:
            cmd = ["resvg", "-V"]
            subprocess.check_call(cmd, stdout=subprocess.PIPE)
            RESVG = ["resvg"]
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

def _checkRSVGConvert():
    try:
        subprocess.check_call("rsvg-convert -v", stdout=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def _checkInkscape():
    try:
        subprocess.check_call("inkscape --version", stdout=subprocess.PIPE, shell=True)
        return True
    except subprocess.CalledProcessError:
        return False



def _convertSVG_resvg(inpath, outpath):
    try:
        # cmd = "resvg {} {}".format(inpath, outpath)
        cmd = RESVG + [inpath, outpath]
        subprocess.check_call(cmd)#, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        return None

    return open(outpath, "rb").read()


def _convertSVG_resvg_stdio(indata):
    _checkRESVGConvert()
    try:
        if isinstance(indata, str):
            indata = indata.encode('utf-8')

        cmd = RESVG +  ["--resources-dir", "./", "-", "-c"]
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outdata, errdata = process.communicate(indata)
        if process.returncode != 0:
            print(f"Error: {errdata.decode('utf-8')}")
            return None

        return outdata

    except subprocess.CalledProcessError as e:
        print(f"Subprocess failed with error {e}")
        return None


def _convertSVG_webkitToPDF(inpath, outpath, output_format):
    if output_format.lower() != "pdf":
        return None

    try:
        cmd = "webkitToPDF {} {}".format(inpath, outpath)
        subprocess.check_call(cmd, shell=True)#, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        return None

    return open(outpath, "rb").read()

def _convertSVG_inkscape(inpath, outpath, output_format):
    options = ""
    output_format = output_format.lower()
    if output_format == "png":
        options = "--export-dpi 150 --export-background white"

    try:
        subprocess.check_call("inkscape {} {} --export-{}={}".format(options, inpath, output_format, outpath), 
            shell=True)
    except subprocess.CalledProcessError as e:
        print("EXPORT ERROR:", str(e))

    return open(outpath, "rb").read()


def _convertSVG_rsvg_convert(inpath, outpath, output_format):
    options = [""]
    output_format = output_format.lower()
    if output_format == "png":
        options = ["-a", "--background-color=white"]

    try:
        cmd = ["rsvg-convert", "-f", output_format] + options + ["-o", outpath, inpath]
        # subprocess.check_call("rsvg-convert -f {} {} -o {} {}".format(output_format, options, outpath, inpath), shell=True)
        process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outdata, errdata = process.communicate()
        if process.returncode != 0:
            print(f"Error: {errdata.decode('utf-8')}")
            return None

        return outdata

    except subprocess.CalledProcessError as e:
        print("EXPORT ERROR:", str(e))

    # return open(outpath, "rb").read()
