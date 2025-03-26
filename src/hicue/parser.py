from .cli.imports import *

class BedRecord:
    def __init__(self, chrom, chromStart, chromEnd, name=None, score=None, strand=None):
        self.chrom = chrom
        self.chromStart = int(chromStart)
        self.chromEnd = int(chromEnd)
        self.name = name
        self.score = score
        self.strand = strand

    def __repr__(self):
        return f"BedRecord(chrom={self.chrom}, chromStart={self.chromStart}, chromEnd={self.chromEnd}, " \
               f"name={self.name}, score={self.score}, strand={self.strand})"


def parse_bed_file(bed_file_path):
    records = []
    
    with open(bed_file_path, 'r') as file:
        for line in file:
            # Skip empty lines or comments
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            
            # BED file typically has at least 3 columns, handle additional ones if present
            chrom = fields[0]
            chromStart = fields[1]
            chromEnd = fields[2]
            
            # Optional fields
            name = fields[3] if len(fields) > 3 else None
            score = fields[4] if len(fields) > 4 else None
            strand = fields[5] if len(fields) > 5 else None

            # Create a BedRecord object and add it to the list
            record = BedRecord(chrom, chromStart, chromEnd, name, score, strand)
            records.append(record)
    
    return records

def parse_gff(in_file):
    """Parses gff file into python dataframe."""
    gff = GFF.parse(in_file)
    gff_df = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand"])
    
    for rec in gff:
        for feature in rec.features:
                loc = feature.location
                gff_tmp = {
                    "Name": feature.qualifiers['Name'][0] if "Name" in feature.qualifiers else feature.id,
                    "Chromosome": rec.id,
                    "Start": loc.start,
                    "End": loc.end,
                    "Strand": loc.strand,
                }
                gff_df = gff_df._append(gff_tmp, ignore_index=True)
    return gff_df

def parse_bed(in_file, default_strand=1):
    """Parses bed file into python dataframe."""
    bed = parse_bed_file(in_file)
    bed_df = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand"])
    
    for rec in bed:
        bed_tmp = {
                    "Name": rec.name,
                    "Chromosome": rec.chrom,
                    "Start": rec.chromStart,
                    "End": rec.chromEnd,
                    "Strand": 1 if rec.strand == '+' else -1 if rec.strand == '-' else default_strand
                }
        bed_df = bed_df._append(bed_tmp, ignore_index=True)
    return bed_df

def genes_in_rec(bed_rec, gff, overlap="flex"):
    """Returns the records from gff included in the bed record bed_rec interval."""
    pos1_rec = min(bed_rec.chromStart, bed_rec.chromEnd)
    pos2_rec = max(bed_rec.chromStart, bed_rec.chromEnd)

    genes = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand"])
    for _, position in gff.iterrows():
        if position["Chromosome"] != bed_rec.chrom:
            continue

        pos1 = min(position["Start"], position["End"])
        pos2 = max(position["Start"], position["End"])
        match overlap:
            case "flex":
                # included if the position overlaps the regions, even not completely
                if (pos1 >= pos1_rec and pos1 <= pos2_rec) or (pos2 >= pos1_rec and pos2 <= pos2_rec):
                    genes = genes._append(position, ignore_index = True)
            case "strict":
                # included only if the position is completely in the region
                if pos1 >= pos1_rec and pos2 <= pos2_rec:
                    genes = genes._append(position, ignore_index = True)

    return genes

def parse_bed_annotated(bed_file, gff_file, overlap="flex"):
    """Annotates each interval in the provided bed file with the found sequences in the gff file."""
    bed = parse_bed_file(bed_file)
    gff = parse_gff(gff_file)
    genes_df = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand"])
    for rec in bed:
        rec_genes = genes_in_rec(rec, gff, overlap=overlap)
        genes_df = pd.concat([genes_df, rec_genes])
    return genes_df

def parse_tracks(tracks, threshold = None, percentage = None, gff=None, default_strand=1, binning=1000):
    """Parses a track file into a positions table. Applies threshold AND percentage if provided. If a gff file has been provided, will compute average track value per record."""
    bigwig_file = pyBigWig.open(tracks)

    percentage_threshold = None
    if percentage != None:
        app_end, quantile = percentage
        match app_end:
            case 'high':
                quantile = 1 - quantile/100
            case 'low':
                quantile = quantile/100
        all_values = np.concatenate([bigwig_file.values(chrom, 0, bigwig_file.chroms(chrom)) for chrom in bigwig_file.chroms().keys()])
        percentage_threshold = np.nanquantile(all_values, quantile)

    if threshold != None:
        min_max, threshold_value = threshold

    if gff != None:
       positions = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand", "Track"])
       for rec in GFF.parse(gff):
            for feature in rec.features:
                    loc = feature.location
                    value = bigwig_file.stats(rec.id, loc.start, loc.end)[0]

                    if threshold != None and value != None:
                        match min_max:
                            case 'min':
                                value = value if value >= threshold_value else None
                            case 'max':
                                value = value if value <= threshold_value else None

                    if percentage_threshold != None and value != None:
                        match app_end:
                            case 'high':
                                value = value if value >= percentage_threshold else None
                            case 'low':
                                value = value if value <= percentage_threshold else None

                    if value != None:
                        pos_tmp = {
                            "Name": feature.qualifiers['Name'][0] if "Name" in feature.qualifiers else feature.id,
                            "Chromosome": rec.id,
                            "Start": loc.start,
                            "End": loc.end,
                            "Strand": loc.strand,
                            "Track": value
                        }
                        positions = positions._append(pos_tmp,  ignore_index = True)
    else:
        positions = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand", "Track"])
        k = 1
        for chrom in bigwig_file.chroms().keys():
            for i in range(0, bigwig_file.chroms(chrom) - binning, binning):
                start = i
                stop = i + binning
                value = bigwig_file.stats(chrom, start, stop)[0]

                if threshold != None and value != None:
                    match min_max:
                        case 'min':
                            value = value if value >= threshold_value else None
                        case 'max':
                            value = value if value <= threshold_value else None

                if percentage_threshold != None and value != None:
                    match app_end:
                        case 'high':
                            value = value if value >= percentage_threshold else None
                        case 'low':
                            value = value if value <= percentage_threshold else None
                if value != None:
                    pos_tmp = {
                        "Name": f"Selection_{k}",
                        "Chromosome": chrom,
                        "Start": start,
                        "End": stop,
                        "Strand": default_strand,
                        "Track": value
                    }
                    positions = positions._append(pos_tmp, ignore_index = True)
                    k += 1

    return positions, bigwig_file



        

