from .imports import *

# global parameters
global position_id
global hash_pos_id

positions_lock = threading.Lock()

def initialize_globals():
    with positions_lock:
        global position_id
        position_id = 0

        global hash_pos_id
        hash_pos_id = {}

class ParserScheduler(threading.Thread):
    
    def __init__(self, file_type, record_type, input_queue, output_queues, is_loop = False, **kwargs):
        super(ParserScheduler, self).__init__(**kwargs)

        self._file_type = file_type
        self._record_type = record_type
        self._is_loop = is_loop
        self._input_queue = input_queue
        self._output_queues = output_queues

        self.start()

    def run(self):
        parser = Parser(self._file_type, self._record_type, self._output_queues[0], self._output_queues[1], is_loop = self._is_loop)
        while True:
            try:
                val = self._input_queue.get(timeout = 10)
            except Empty:
                break
            if val == 'DONE':
                break
        
            parser.parse_line(val)
                        

class Parser():
    def __init__(self, file_type, record_type, position_queue, pairing_queue, is_loop = False):
        self._file_type = file_type
        self._position_queue = position_queue
        self._pairing_queue = pairing_queue
        self._record_type = record_type
        self._is_loop = is_loop

    @staticmethod
    def parse_bed2d_line(line):
        """Parses a line from a BED2D formated file"""
        fields = line.split('\t')

        record1 = {
            "Chromosome": fields[0],
            "Start": int(fields[1]) + 1, # adapting the base for cooler (1-based)
            "End": int(fields[2]) + 1, # adapting the base for cooler (1-based)
            "Name": fields[3] if len(fields) > 6 else None,
            "Strand": -1 if len(fields) > 10 and fields[5] == '-' else 1
        }
        hash1 = f"{record1['Chromosome']},{record1['Start']},{record1['End']},{record1['Strand']}".encode("ascii")

        k = len(fields) // 2 - 1 

        record2 = {
            "Chromosome": fields[k + 1],
            "Start": int(fields[k + 2]) + 1, # adapting the base for cooler (1-based)
            "End": int(fields[k + 3]) + 1, # adapting the base for cooler (1-based)
            "Name": fields[k + 4] if len(fields) > 6 else None,
            "Strand": -1 if len(fields) > 10 and fields[k + 6] == '-' else 1
        }
        hash2 = f"{record2['Chromosome']},{record2['Start']},{record2['End']},{record2['Strand']}".encode("ascii")
        
        return [record1, record2], [hash1, hash2]

    @staticmethod
    def parse_gff_line(line, record_type = "gene"):
        """Parses a line from a GFF formated file"""
        fields = line.split('\t')

        if record_type and fields[2] != record_type:
            return None

        subfields = [f.split("=") for f in fields[8].split(";")]
        qualifiers = {sf[0]:sf[1] for sf in subfields}
        record = {
            "Chromosome": fields[0],
            "Start": int(fields[3]),
            "End": int(fields[4]),
            "Name": qualifiers["Name"] if "Name" in qualifiers else qualifiers["ID"],
            "Strand": -1 if fields[6] == '-' else 1 if fields[6] == '+' else 0
        }
        hash = f"{record['Chromosome']},{record['Start']},{record['End']},{record['Strand']}".encode("ascii")

        return [record], [hash]

    @staticmethod
    def parse_bed_line(line):
        """Parses a line from a BED formated file"""
        fields = line.split('\t')

        record = {
            "Chromosome": fields[0],
            "Start": int(fields[1]) + 1,
            "End": int(fields[2]) + 1,
            "Name": fields[3] if len(fields) > 3 else None,
            "Strand": -1 if len(fields) > 5 and fields[5] == '-' else 1
        }
        hash = f"{record['Chromosome']},{record['Start']},{record['End']},{record['Strand']}".encode("ascii")


        return [record], [hash]

    def parse_line(self, line):
        # parsing line
        line_positions = []
        match self._file_type:
            case "bed":
                line_positions, line_hashes = self.parse_bed_line(line)
            case "bed2d":
                line_positions, line_hashes = self.parse_bed2d_line(line)
            case "gff":
                value = self.parse_gff_line(line, record_type=self._record_type)
                if not value:
                    return 
                line_positions, line_hashes = value

        # retrieving positions indexes
        indexes = []
        for i in range(len(line_hashes)):
            hash = line_hashes[i]
            index = -1
            put_in_queue = False
            with positions_lock:
                global hash_pos_id
                if hash in hash_pos_id:
                    index = hash_pos_id[hash]
                else:
                    global position_id
                    index = position_id
                    position_id += 1
                    hash_pos_id[hash] = index
                    put_in_queue = True
            if put_in_queue:
                if not line_positions[i]["Name"]:
                    line_positions[i]["Name"] = f"locus_{index}"
                self._position_queue.put((index, line_positions[i]))
            indexes.append(index)

        # writing the line index pair # TODO add the is_loop option (all indexes are matched with the new one)
        if len(indexes) == 1 and self._is_loop:
            for k in range(indexes[0]):
                index_pair = (k, indexes[0])
                self._pairing_queue.put(index_pair)
        else:
            index_pair = (indexes[0], indexes[0]) if len(indexes) < 2 else (indexes[0], indexes[1])
            self._pairing_queue.put(index_pair)
        


        
