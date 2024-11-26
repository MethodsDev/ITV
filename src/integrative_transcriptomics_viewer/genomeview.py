import collections
from cairosvg import svg2png
import xml.etree.ElementTree as ET
import ipywidgets as widgets
from ipywidgets.embed import embed_minimal_html, dependency_state

from integrative_transcriptomics_viewer.svg import Renderer, SVG
from integrative_transcriptomics_viewer.export import SvgSplitter, _convertSVG_resvg_stdio


class Document:
    doc_counter = 0
    svg_header = """<svg version="1.1" id="{doc_counter}" width="{width}" height="{height}" shape-rendering="crispEdges" xmlns="http://www.w3.org/2000/svg">"""
    svg_footer = """</svg>"""
        
    def __init__(self, width):
        self.id = "itv_svg_" + str(self.doc_counter)
        Document.doc_counter += 1

        self.elements = []
        self.width = width
        self.height = None
        self.hide = False
        self.renderer = SVG()
        
        self.margin_x = 5
        self.margin_y = 5

        self.between_views = 5

    def add_view(self, view):
        self.elements.append(view)
        
    def add_track(self, track):
        for element in self.elements:
            try:
                element.add_track(track)
                break
            except AttributeError:
                continue
        else:
            raise Exception("No GenomeView found in Document")

    def get_tracks(self, name=None):
        matching = []
        for element in self.elements:
            try:
                matching.extend(element.get_tracks(name))
            except AttributeError:
                pass
        return matching

    def layout(self):
        self.view_width = self.width - self.margin_x*2
        for element in self.elements:
            element.layout(self.view_width)
        
    def render(self):
        self.layout()

        self.height = sum(element.height+self.between_views for element in self.elements) + self.margin_y*2
        
        if self.hide:
            yield self.svg_header.format(doc_counter=self.id, height=0, width=0)
        else:
            yield self.svg_header.format(doc_counter=self.id, height=self.height, width=self.width)
        
        cury = self.margin_y
        for element in self.elements:
            renderer = Renderer(self.renderer, self.margin_x, cury, self.view_width, element.height)
            yield from renderer.render(element)
            cury += element.height + self.between_views
            
        yield self.svg_footer

    def _repr_html_(self):
        #return widgets.HTML(self._repr_svg__())
        svg_content = self._repr_svg__()
        custom_style = """
        <style>
            .custom-svg-container svg {
                max-width: none;
                height: auto;
            }
        </style>
        """
        if self.hide:
            html_content = f"""
            {svg_content}
        <div class="custom-svg-container" style="overflow-x: auto; overflow-y: hidden; width: 100%;">
            
        </div>
        """

        else:
            html_content = f"""
        <div class="custom-svg-container" style="overflow-x: auto; overflow-y: hidden; width: 100%;">
            {svg_content}
        </div>
        """

        html_content += custom_style

        return html_content
    
    def _repr_svg__(self):
        return "\n".join(self.render())

    def get_images(self, output_format="svg"):
        svg = self._repr_svg__()
        if output_format == "svg":
            return [svg]
        else:  # "png"
            root_svg = ET.fromstring(svg)
            svg_splitter = SvgSplitter(root_svg)
            svg_splitter.split_svg(root_svg, max_height = 10000)

            pngs = []
            for split_svg in svg_splitter.get_splits():
                pngs.append(_convertSVG_resvg_stdio(ET.tostring(split_svg.getroot(), encoding='utf-8')))
                # pngs.append(svg2png(bytestring=ET.tostring(split_svg.getroot(), encoding='utf8')))
            return pngs

    def get_widget(self, output_format="svg"):
        if output_format == "png":
            return widgets.VBox([widgets.Image(value=data, format=output_format) for data in self.get_images(output_format)])
        else: # "svg"
            return widgets.HTML(self._repr_svg__())

    def render_html(self, output_path, output_format="svg", title=""):
        view = self.get_widget(output_format=output_format)
        embed_minimal_html(output_path, views=view, state=dependency_state(view), title=title)


class ViewRow:
    def __init__(self, name=None):
        self.name = name
        self.views = []

        self.width = None
        self.height = None

        self.space_between = 5
    
    def add_view(self, view):
        self.views.append(view)
    
    def get_views(self, name=None):
        matching = []
        for view in self.views:
            if name is None or view.name==name:
                matching.extend(view)
        return matching

    def get_tracks(self, name=None):
        matching = []
        for view in self.views:
            try:
                matching.extend(view.get_tracks(name))
            except AttributeError:
                pass
        return matching

    def layout(self, width):
        self.width = width
        n_views = len(self.views)

        self.each_width = {}
        
        self.all_views_have_width = True
        for view in self.views:
            if view.pixel_width is None:
                self.all_views_have_width = False

        self.height = 0
        if not self.all_views_have_width:
            self.each_width = (self.width - self.space_between*(n_views-1)) / n_views

            for view in self.views:
                view.layout(self.each_width)
                self.height = max(self.height, view.height)

        else:
            for view in self.views:
                view.layout(view.pixel_width)
                self.height = max(self.height, view.height)
    


    def render(self, renderer):
        curx = 0
        if self.all_views_have_width:
            for view in self.views:    
                subrenderer = renderer.subrenderer(x=curx, width=view.pixel_width, height=view.height)
                yield from subrenderer.render(view)
                curx += view.pixel_width + self.space_between

        else:
            for view in self.views:    
                subrenderer = renderer.subrenderer(x=curx, width=self.each_width, height=view.height)
                yield from subrenderer.render(view)
                curx += self.each_width + self.space_between


class GenomeView:
    def __init__(self, chrom, start, end, strand, source=None, name=None):
        self.name = name
        self.tracks = []

        self.scale = Scale(chrom, start, end, strand, source)

        self.pixel_width = None
        # self.pixel_height = None

        self.margin_y = 2  # 10
    
    def add_track(self, track):
        self.tracks.append(track)

    def get_tracks(self, name=None):
        matching = []
        for track in self.tracks:
            if name is None or track.name == name:
                matching.append(track)
        return matching
        
    def layout(self, width):
        self.pixel_width = width
        self.scale.pixel_width = width

        self.height = 0
        for track in self.tracks:
            track.layout(self.scale)
            self.height += track.height + self.margin_y
    
    def render(self, renderer):
        cury = 0
        for track in self.tracks:
            subrenderer = renderer.subrenderer(y=cury, height=track.height)
            yield from subrenderer.render(track)
            cury += track.height + self.margin_y
        
    
    
class Scale:
    """
    Maintains information about a projection of a specific genomic
    interval into screen coordinates.

    That is, we're interested in visualizing an interval (chrom:start-end) 
    on a canvas of a specified pixel width. The scale enables converting
    genomic coordinates into the display coordinates.
    """
    def __init__(self, chrom, start, end, strand, source):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

        self.pixel_width = None
        self._param = None

        self.source = source
        self._seq = None

        if self.start >= self.end:
            raise ValueError("End coordinate must be greater than start coordinate; you specified {}:{}-{}".format(chrom, start, end))

    def __len__(self):
        return self.end - self.start + 1

    def _setup(self):
        if self._param == self.pixel_width: return
        self._param = self.pixel_width
        self._seq = None

        nt_width = self.end - self.start
        
        self.bases_per_pixel = nt_width / self.pixel_width

    def topixels(self, genomic_position):
        """
        Converts a genomic position to a pixel location in the current
        coordinate system.
        """
        self._setup()

        pos = (genomic_position - self.start) / float(self.bases_per_pixel)
        return pos

    def relpixels(self, genomic_size):
        """
        Takes a genomic length (ie, a number of basepairs) and converts
        it to a relative screen length in pixels.
        """
        self._setup()

        dist = genomic_size / float(self.bases_per_pixel)
        return dist
    
    def get_seq(self, start=None, end=None, strand=None):
        """
        Gets the nucleotide sequence of an interval. By default, returns the 
        sequence for the current genomic interval.
        """
        self._setup()

        assert self.source is not None
        if strand is not None and strand != self.strand:
            raise Exception("ack")

        if start is None:
            start = self.start
        if end is None:
            end = self.end
        
        assert start >= self.start
        assert end <= self.end

        if self._seq is None:
            self._seq = self.source.get_seq(self.chrom, self.start, self.end, "+").upper()

        cur_seq = self._seq[start-self.start:end-self.start]

        return cur_seq

