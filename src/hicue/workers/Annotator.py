from .imports import *

def initialize_globals():
    pass

class AnnotatorScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, gff = None, tracks = None, overlap = "strict", save_to = "", gff_type = "gene", **kwargs):
        super(AnnotatorScheduler, self).__init__(**kwargs)

        self._gff = gff or None
        self._gff_type = gff_type or None
        self._tracks = tracks or None
        self._overlap = overlap
        self._save_to = save_to

        self._input_queue = input_queue
        self._output_queues = output_queues

        self.start()

    def run(self):

        annotators = []
        if self._gff:
            annotators.append(GffAnnotator(gff = self._gff, overlap = self._overlap, gff_type = self._gff_type, save_to = self._save_to))
        if self._tracks:
            annotators.append(TracksAnnotator(tracks = self._tracks))
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break
            index, position = val
            for annotator in annotators:
                position = annotator.annotate_position(index, position)

            if position:
                for queue in self._output_queues:
                    queue.put((index, position)) 

class GffAnnotator():
    def __init__(self, gff, overlap, gff_type, save_to):
        self._gff = gff
        self._overlap = overlap
        self._save_to = save_to
        self._gff_type = gff_type

    @staticmethod
    def format_gff_feature(chromosome, feature):
        """Formats a GFF feature into a HiCue position."""
        loc = feature.location
        position = {
            "Name": feature.qualifiers['Name'][0] if "Name" in feature.qualifiers else feature.id,
            "Chromosome": chromosome,
            "Start": loc.start,
            "End": loc.end,
            "Strand": loc.strand
        }
        return position

    @staticmethod
    def compile_gff_annotation(annotation):
        """Compiles a list of positions into a single name and returns it as a string. If provided, will save the positions annotation in append mode in save_to."""
        known_annotation = {}
        unknown_annotation = []
        strands = []

        for annotated_pos in annotation:
            # formated genes have three lowercase characters and might have an uppercase character afterward. If more character, is a code.
            name = annotated_pos["Name"]
            strands.append(annotated_pos["Strand"])
            if len(name) <= 4 and name[:3].islower():
                prefix = name[:3]
                identificator = name[3] if len(name) == 4 else ""
                if prefix in known_annotation:
                    known_annotation[prefix].append(identificator)
                else:
                    known_annotation[prefix] = [identificator]
            else:
                unknown_annotation.append(name)
        concat_known = "/".join([prefix + "".join(known_annotation[prefix]) for prefix in known_annotation])
        concat_unknown = "/" + "/".join(unknown_annotation) if len(concat_known) > 0 else "/".join(unknown_annotation)
        concat = concat_known + concat_unknown
        sum_strands = sum(strands)
        nb_annotated_pos = len(annotation)
        return concat, (-1 if sum_strands== -nb_annotated_pos else 1 if sum_strands == nb_annotated_pos else 0)
            
    def annotate_position(self, index, position):
        """Annotates a position from a GFF file."""
        chromosome = position["Chromosome"]
        if chromosome not in self._gff:
            return 
        if self._gff_type:
            limit_info = dict(gff_type=["gene"])
            position_gff = GFF.parse(self._gff[chromosome], limit_info=limit_info)
        else:
            position_gff = GFF.parse(self._gff[chromosome])
        start = position["Start"]
        stop = position["End"]
        
        annotation = []
        outfile = open(f"{self._save_to}/position_{index}_annotation.bed", 'w') if len(self._save_to) > 0 else None
        for rec in position_gff:
            for feature in rec.features:
                rec_start = feature.location.start + 1 # the GFF parser corrects the 1 based for the location interval to be used as limits for sequence.
                rec_stop = feature.location.end
                match self._overlap:
                    case "strict":
                        if start <= rec_start and stop >= rec_stop:
                            annotation.append(self.format_gff_feature(rec.id, feature))
                            if outfile:
                                outfile.write(f"{rec.id}\t{rec_start}\t{rec_stop}\t{feature.qualifiers['Name'][0] if 'Name' in feature.qualifiers else feature.id}\t{feature.location.strand}\n")
                    case "flex":
                        if start < rec_stop and stop > rec_start:
                            annotation.append(self.format_gff_feature(rec.id, feature))
                            if outfile:
                                outfile.write(f"{rec.id}\t{rec_start}\t{rec_stop}\t{feature.qualifiers['Name'][0] if 'Name' in feature.qualifiers else feature.id}\t{feature.location.strand}\n")
                if stop < rec_start:
                    break
        if outfile:
            outfile.close()
        name, strand = self.compile_gff_annotation(annotation)
        position["Name"] = name if len(name) > 0 else position["Name"]
        position["Strand"] = strand

        return position

class TracksAnnotator():
    def __init__(self, tracks):
        self._tracks = tracks
            
    def annotate_position(self, index, position, threshold = None):
        """Annotates a position from a track file. Adds a columns called Tracks.
        If the track requires selecting that can be computed at this step, the position can be discarded."""  #TODO: add tracks options
        position["Tracks"]  = self._tracks.stats(position["Chromosome"], position["Start"], position["End"])[0]
        if threshold:
            threshold_value, min_max = threshold
            return position if (min_max == "min" and position["Tracks"] >= threshold_value) or (min_max == "max" and position["Tracks"] <= threshold_value) else None
        return position
        
