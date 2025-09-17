from .imports import *

def initialize_globals():
    pass

class SeparatorScheduler(threading.Thread):
    def __init__(self, input_queue, output_queue, **separgs):
        super(SeparatorScheduler, self).__init__()
        self._input_queue = input_queue
        temp_queue = output_queue
        if type(temp_queue) != list:
            temp_queue = [temp_queue]
        self._output_queues = temp_queue
        self._separgs = separgs
        self.start()

    def run(self):
        separator = Separator(**self._separgs)
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break
            
            pos1, pos2, i, j, distance, circ1, circ2 = val
            sep_id = separator.sort(pos1, pos2, distance, is_paired = (i != j))
            if len(sep_id) > 0:
                for output_queue in self._output_queues:
                    output_values = (i, j, distance, circ1, circ2, sep_id)
                    output_queue.put(output_values)

class Separator():
    """Class implementing the separate_by option: computes the group of a pair of positions"""
    def __init__(self, separate_by, overlap = "strict", center = "start", contact_range = [], separate_regions = None, has_trans = False):
        self._separate_by = separate_by.split(",") if separate_by and len(separate_by) > 0 else []
        self._overlap = overlap
        self._center = center
        self._contact_range = contact_range if contact_range else []
        self._has_trans = has_trans
        self._range_list = np.array(list(
            range(self._contact_range[0], self._contact_range[1], self._contact_range[2])
        )) if len(self._contact_range) > 0 else []
        self._range_dict = {
            self._range_list[i]: self._range_list[i+1] for i in range(len(self._range_list) - 1)
        }
        self._regions = None
        if "regions" in self._separate_by:
            if separate_regions != None:
                self._regions = pd.read_csv(separate_regions, header=0, index_col=None)
                for col in ["Id", "Chromosome", "Start", "End"]:
                    if col not in self._regions.columns:
                        sys.stderr.write("Region separation file error: please ensure the file is in csv format with the columns Id, Chromosome, Start, End. This separation will not be treated.\n")
                        self._separate_by.remove("regions")
            else:
                sys.stderr.write("Region separation file error: please ensure the file has been provided in csv format with the columns Id, Chromosome, Start, End. This separation will not be treated.\n")
                self._separate_by.remove("regions")
        
    @staticmethod
    def is_in_region(position, region, overlap = "strict"):
        """Assesses if a position is contained (strictly or not depending on the overlap parameter) in a region."""
        if position["Chromosome"] != region["Chromosome"]:
            return False
        
        pos1 = min(position["Start"], position["End"])
        pos2 = max(position["Start"], position["End"])
        pos1_reg = min(region["Start"], region["End"])
        pos2_reg = max(region["Start"], region["End"])
        
        match overlap:
            case "flex":
                # included if the position overlaps the regions, even not completely
                if (pos1 >= pos1_reg and pos1 <= pos2_reg) or (pos2 >= pos1_reg and pos2 <= pos2_reg):
                    return True
            case "strict":
                # included only if the position is completely in the region
                if pos1 >= pos1_reg and pos2 <= pos2_reg:
                    return True
        return False
    
    @staticmethod
    def separate_direction(pos1, pos2):
        """Returns the directionality group of the pair. None if undirected or different chromosomes."""
        if pos1["Chromosome"] == pos2["Chromosome"]:
            if pos1["Strand"] == 1:
                return "convergent" if pos2["Strand"] == -1 else "forward"
            if pos1["Strand"] == -1:
                return "divergent" if pos2["Strand"] == 1 else "reverse"
        return None
    
    @staticmethod
    def separate_chroms(pos1, pos2):
        """Returns the chromosome group of the pair. Both chromsomes in alphabetical order separaed by a _, if both are not in the same chromosome."""
        return pos1["Chromosome"] if pos1["Chromosome"] == pos2["Chromosome"] else "_".join(sorted([pos1["Chromosome"], pos2["Chromosome"]]))
    
    def separate_cis_trans(self, pos1, pos2):
        """Returns the pair's contact type: cis in the same chromosome, trans in different ones."""
        return "cis" if pos1["Chromosome"] == pos2["Chromosome"] else "trans" if self._has_trans else None
    
    def separate_regions(self, pos1, pos2, is_paired):
        """Returns the region group of the pair. None if not in any region, or both positions are not in the same region."""
        for region_id in np.unique(self._regions["Id"]):
            selected = self._regions[self._regions["Id"] == region_id]
            pos1_in_region = False
            for _,region in selected[selected["Chromosome"] == pos1["Chromosome"]].iterrows():
                if self.is_in_region(pos1, region):
                    if not is_paired:
                        return region_id
                    pos1_in_region = True
                    break
            if pos1_in_region:
                for _, region in selected[selected["Chromosome"] == pos2["Chromosome"]].iterrows():
                    if self.is_in_region(pos2, region):
                        return region_id
        return None
        
    def separate_distance(self, pos1, pos2, distance):
        """Returns the distance interval in which the pair's distance is located. If not included in the contact range, returns None."""
        if not distance:
            return None
        pair_dist = abs(distance)
        if self._contact_range:
            pair_dist_bin = ((pair_dist - self._contact_range[0]) // self._contact_range[2]) * self._contact_range[2] + self._contact_range[0]
            if pair_dist_bin in self._range_list[:-1]:
                return f"{pair_dist_bin}-{self._range_dict[pair_dist_bin]}"
        return None
            
    def sort(self, pos1, pos2, distance, is_paired = False):
        """Applies all the separations to the pair. Returns its group."""
        group_list = []
        
        if not is_paired and len(self._separate_by) == 0:
            return "default"

        # separate_by
        for sep in self._separate_by:
            group = None
            match sep:
                case "direct":
                    group = self.separate_direction(pos1, pos2)
                case "regions":
                    group = self.separate_regions(pos1, pos2, is_paired)
                case "chroms":
                    group = self.separate_chroms(pos1, pos2)
            if group == None:
                return ""
            group_list.append(group)
            
        # contact_range
        if is_paired and len(self._contact_range) > 0 :
            group = self.separate_distance(pos1, pos2, distance)
            if group == None:
                return ""
            group_list.append(group)

        # cis_trans (default)
        if is_paired:
            group = self.separate_cis_trans(pos1, pos2)
            if group == None:
                return ""
            group_list.append(group)
        
        return "_".join(group_list)
