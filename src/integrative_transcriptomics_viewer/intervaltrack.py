from integrative_transcriptomics_viewer.track import Track
# import random


class Interval:
    def __init__(self, id_, chrom, start, end, strand=True, label=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.id = id_
        self.label = label
        
        #if isinstance(strand, bool):
        #    strand = {True:"+", False:"-"}[strand]
        self.strand = strand

    def overlaps(self, other, ignore_strand=True):
        if self.chrom != other.chrom:
            return False
        if other.start > self.end or self.start > other.end:
            return False
        if not ignore_strand and self.strand != other.strand:
            return False
        return True

    def __repr__(self):
        return "{}:{:,}-{:,}{}".format(self.chrom, self.start, self.end, self.strand)
        
# def color_by_strand(interval):
#     # brightness = 0.2 + (cur_reads[0].mapq/40.0*0.8)
#     color = "#E89E9D"
#     if interval.strand == "-":
#         color = "#8C8FCE"
#     return color

def color_by_strand(interval):
    # brightness = 0.2 + (cur_reads[0].mapq/40.0*0.8)
    if interval.strand:
        return "#E89E9D"
    return "#8C8FCE"


class IntervalTrack(Track):
    def __init__(self, intervals, name=None):
        super().__init__(name)
        self.rows = []
        self.intervals_to_rows = {}
        
        self.row_height = 8
        self.margin_x = 15
        self.margin_y = 2
        self.label_distance = 3
        
        self.intervals = intervals

        self.vertical_layout = False
        self.strand_specific = False
        self.color_fn = color_by_strand

    def layout_interval(self, interval):
        row = None
        interval_start = self.scale.topixels(interval.start)

        if self.strand_specific and interval.strand != self.scale.strand:
            return

        if self.vertical_layout:
            row = len(self.rows)
            # if not self.max_depth or (self.max_depth and row <= self.max_depth):
            self.rows.append(None)
            # else:
            #     return
        else:
            # if haven't reached max number of reads to display, we can try to fit it on an existing row, max_depth doesn't need to be checked here because the populated rows already are within that limit
            # if not self.max_reads or len(self.intervals_to_rows) < self.max_reads:  
            for rowi, row_end in enumerate(self.rows):
                if interval_start > row_end:  # could keep track of row_start as well, in case of random sorted
                    row = rowi
                    break
            # if row is None:
            #     if (not self.max_reads and not self.max_depth) or (self.max_depth and len(self.rows) < self.max_depth) or (self.max_reads and len(self.intervals_to_rows) < self.max_reads):
            row = len(self.rows)
            self.rows.append(None)
            #     else:
            #         return

            new_end = self.scale.topixels(interval.end) + self.margin_x
            if interval.label is not None:
                new_end += len(interval.label) * self.row_height * 0.75
            self.rows[row] = new_end

        assert not interval.id in self.intervals_to_rows
        self.intervals_to_rows[interval.id] = row


    def layout(self, scale):
        super().layout(scale)

        self.rows = []
        self.intervals_to_rows = {}

        # if self.max_depth:
        #     intervals = [_ for _ in self.intervals]
        #     random.shuffle(intervals)
        #     for interval in intervals:
        #         self.layout_interval(interval) #, max_rows = self.max_depth)
        #     if len(self.rows) > self.max_depth:
        #         self.rows = self.rows[:self.max_depth]
        # elif self.max_reads and len(self.intervals) > self.max_reads:  # max reads and it's more than the number of reads
        #     # implement resevoir sample 
        #     #intervals = self.intervals[:self.max_reads]
        #     #for inter in self.intervals[self.max_reads:]:
        #     pass
        # else:
        for interval in self.intervals:
            self.layout_interval(interval)
        
        self.height = max(1, len(self.rows)) * (self.row_height + self.margin_y)
    

    def draw_interval(self, renderer, interval, extra_args={}):
        start = self.scale.topixels(interval.start)
        end = self.scale.topixels(interval.end)
        
        # interval might not be mapped due to sampling  
        if interval.id not in self.intervals_to_rows:
            return
        row = self.intervals_to_rows[interval.id]
        top = row*(self.row_height+self.margin_y)
        
        color = self.color_fn(interval)

        if "label" not in extra_args:
            if interval.label is None:
                extra_args["label"] = interval.id
            else:
                extra_args["label"] = interval.label

        if "title" not in extra_args:
            if interval.label is None:
                extra_args["title"] = interval.id
            else:
                extra_args["title"] = interval.label

        # yield from renderer.rect(start, top, end-start, self.row_height, fill=color, 
        #     **{"stroke":"none", "id":temp_label})

        if interval.strand is None:
            yield from renderer.rect(start, top, end-start, self.row_height, fill=color, 
                **{"stroke":"none", "id":temp_label})
        else:
            arrow_width = min(self.row_height / 2, self.margin_x * 0.7, self.scale.relpixels(30))
            direction = "right" if interval.strand else "left"

            yield from renderer.block_arrow(start, top, end-start, self.row_height, 
                arrow_width=arrow_width, direction=direction,
                fill=color, **{"stroke":"none", "id":extra_args["label"], "title":extra_args["title"]})

        if interval.label is not None:
            yield from renderer.text(end+self.label_distance, top+self.row_height-2, interval.label, anchor="start")
        
    def render(self, renderer):
        for interval in self.intervals:
            yield from self.draw_interval(renderer, interval)
            
        yield from self.render_label(renderer)
