import collections
import numpy

from integrative_transcriptomics_viewer.track import Track
from integrative_transcriptomics_viewer.utilities import match_chrom_format

COLORS = ["blue", "red", "green", "black"]
BINNED_COLORS = ["#1F618D", "#b2182b", "#f4a582", "#9B59B6", "#85929E", 
                 "#1c9099", "#74add1", "#053061", "#1b7837", "#b8e186",
                 "#bebada", "#fed976", "#e7298a", "#47E3FF", "#F6222E",
                 "#771155", "orange",  "#6347FF", "#A93226", "#270e26",
                 "#b8bc53", "#5628ce", "#fa909c", "#8ff331", "#FF6347",
                 "#6347FF", "#556270", "#4ECDC4", "#C7F464", "#FF6B6B",
                 "#C44D58", "#E3FF47", "#FF4787", "#771155", "#AA4488", 
                 "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777",
                 "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77",
                 "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455",
                 "#DD7788", "#FFC312", "#C4E538", "#12CBC4", "#FDA7DF",
                 "#ED4C67", "#F79F1F", "#A3CB38", "#1289A7", "#D980FA",
                 "#B53471", "#EE5A24", "#009432", "#0652DD", "#9980FA",
                 "#833471", "#EA2027", "#006266", "#1B1464", "#5758BB",
                 "#6F1E51"]
SECONDARY_COLORS = ["#AAAAAA", "#AFFFFF", "#FFAFFF", "#FFFFAF"]
SECONDARY_COLORS = ["#000000", "#AFFFFF", "#FFAFFF", "#FFFFAF"]


class Series:
    def __init__(self, x, y, color=None, label=None):
        x, y = zip(*sorted(zip(x,y)))
        self.x = x
        self.y = y
        self.color = color
        self.label = label
        

class GraphTrack(Track):
    """
    Visualizes quantitative data as a line across coordinates within the current genomic view.

    One or more datasets can be visualized (with different colors) on the same track using the
    ``add_series()`` method.
    """
    def __init__(self, name=None, x=None, y=None):
        super().__init__(name)

        self.series = collections.OrderedDict()
        
        self.min_y = 0
        self.max_y = 0
        
        if x is not None:
            self.add_series(x, y)

        self.height = 100
        self.ymargin = 5

        self.fill_coverage = False
        self.draw_y_axis = True
        
    def add_series(self, x, y, color=None, label=None):
        """
        Add a dataset corresponding to a single line in the track (ie, a "series"). Note that
        while a single GraphTrack can visualize multiple datasets, they are all plotted on 
        the same y-axis and so should share the same units.

        Arguments:
            x: a list of genomic coordinates
            y: a list of data values; each y-value must correspond to a single x-value
            color: an SVG color for the line being plotted
            label: an optional text label for the graph being plotted (currently unused)
        """
        if label is None:
            label = "series_{}".format(len(self.series))
            
        assert label not in self.series

        x = numpy.asarray(x)
        y = numpy.asarray(y)

        if color is None:
            color = COLORS[len(self.series) % len(COLORS)]
            
        self.series[label] = Series(x, y, color, label)

        self.min_y = min(self.min_y, y[numpy.isfinite(y)].min())
        self.max_y = max(self.max_y, y[numpy.isfinite(y)].max())

    def ytopixels(self, yval):
        height = self.max_y - self.min_y
        if height == 0:  # when no reads to plot so coverage is empty
            return self.height
        return self.height - ((yval - self.min_y) / height * (self.height-2*self.ymargin) + self.ymargin)
        
    def render(self, renderer):
        for label, series in self.series.items():
            if not self.fill_coverage:
                for i in range(len(series.x)-1):
                    if any(numpy.isnan(series.x[i:i+2])) or any(numpy.isnan(series.y[i:i+2])):
                        continue
                    x1 = self.scale.topixels(series.x[i])
                    x2 = self.scale.topixels(series.x[i+1])
                    y1 = self.ytopixels(series.y[i])
                    y2 = self.ytopixels(series.y[i+1])
                    
                    yield from renderer.line(x1, y1, x2, y1, 
                        **{"stroke-width":1, "stroke":series.color, "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
                    yield from renderer.line(x2, y1, x2, y2, 
                        **{"stroke-width":1, "stroke":series.color, "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
            else:
                current_min_y = 0  # keeping track of min because higher coverage is lower y value (closer to top)

                # need to add checks that series.x[:somevalue] is not < scale.start or continue in loop
                #      but keep track of last series.y for where to start
                # need to add checks that series.x[somevalue:] is not > scale.end or break in loop
                #      and then loop path back to start

                full_path = "<path d=\"M "
                started = False
                starting_x = 0
                starting_y = 0
                prev_y = 0
                # y1 = self.ytopixels(0) + renderer.y  # in case of empty coverage, so that the "past scale end" check doesn't error

                for (x, y) in zip(series.x, series.y):
                    if x < self.scale.start:
                        prev_y = y
                        continue

                    if x >= self.scale.end - 1:
                        if started:
                            # if coverage wasn't 0 before going out of bound
                            if y1 != self.ytopixels(0) + renderer.y:
                                x1 = self.scale.topixels(self.scale.end) + renderer.x
                                full_path += f" L {x1:.2f} {y1:.2f}"
                                current_min_y = min(current_min_y, y1)

                            y1 = self.ytopixels(0) + renderer.y
                            full_path += f" L {x1:.2f} {y1:.2f}"
                        break

                    if not started:
                        # starting anchor points at (0, 0) and (0, prev_y)
                        if prev_y > 0:
                            x1 = self.scale.topixels(self.scale.start) + renderer.x
                            y1 = self.ytopixels(0) + renderer.y
                            full_path += f"{x1:.2f} {y1:.2f}"
                            starting_x = x1
                            starting_y = y1

                            y1 = self.ytopixels(prev_y) + renderer.y
                            full_path += f" L {x1:.2f} {y1:.2f}"

                        # starting anchor point at (x, 0)
                        else:
                            y1 = self.ytopixels(0) + renderer.y
                            starting_y = y1

                        current_min_y = min(current_min_y, y1)
                        started = True

                    # else:
                    x1 = self.scale.topixels(x) + renderer.x
                    if len(full_path) > 11:
                        full_path += f" L {x1:.2f} {y1:.2f}"
                    else:
                        full_path += f"{x1:.2f} {y1:.2f}"
                        starting_x = x1

                    y1 = self.ytopixels(y) + renderer.y
                    full_path += f" L {x1:.2f} {y1:.2f}"
                    current_min_y = min(current_min_y, y1)

                if len(full_path) > 11:
                    full_path += f" L {starting_x:.2f} {starting_y:.2f}\" xcenter=\"{((self.ytopixels(0) + current_min_y)/2):.2f}\" stroke=\"{series.color}\" fill=\"{series.color}\" stroke_width=\"1\"></path>"
                else:
                    full_path = ""

                yield full_path

        if self.draw_y_axis:
            # since the labels are drawn at the top of the ticks, let's make sure the top tick/label is 
            # more than 12 pixels from the top of the track so it doesn't get clipped
            # TODO: this ignores the margin, as of now
            axis_max_y = self.min_y + (self.max_y - self.min_y) * (1-7/self.height)

            # ticks = get_ticks(self.min_y, axis_max_y, 4)
            ticks = numpy.linspace(self.min_y, axis_max_y, 4)

            yield from renderer.line(1, self.ytopixels(ticks[0]), 1, self.ytopixels(ticks[-1]), 
                                     **{"stroke-width":2, "stroke":"gray", "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
            for tick in ticks:
                if self.max_y > 1_000:
                    label = "{:.1g}".format(tick)
                elif self.max_y < 1:
                    label = "{:.1f}".format(tick)
                else:
                    label = "{:,.0f}".format(tick)

                y = self.ytopixels(tick)
                yield from renderer.line(1, y, 10, y, 
                                         **{"stroke-width":2, "stroke":"gray", "stroke-linecap":"square", "shape-rendering":"geometricPrecision"})
                yield from renderer.text(14, y, label, anchor="start", fill="gray")
                
        for x in self.render_label(renderer):
            yield x


class BigWigTrack(GraphTrack):
    """
    Visualizes continuous-valued data from a bigwig file. Requires pyBigWig
    module to be installed. Supports using web URLs as well as file paths.
    """
    def __init__(self, path, nbins=1000, name=None):
        super().__init__(name)
        
        import pyBigWig
        self.bigwig = pyBigWig.open(path)
        self.nbins = 1000

    def layout(self, scale):
        super().layout(scale)
        x = []
        y = []
        binsize = max(1, int((scale.end-scale.start) / self.nbins))
        
        chrom = match_chrom_format(scale.chrom, self.bigwig.chroms().keys())

        for i in range(scale.start, scale.end, binsize):
            value = self.bigwig.stats(chrom, i, i+binsize)[0]
            value = 0.0 if value==None else value
            x.append(i+binsize/2)
            y.append(value)
        
        self.series = {"vals":Series(x, y, color="black")}
        
        self.min_y = min(y)
        self.max_y = max(y)
