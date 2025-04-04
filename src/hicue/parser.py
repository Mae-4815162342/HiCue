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

def parse_bed2d_file(bed_file_path):
    records = []
    
    with open(bed_file_path, 'r') as file:
        for line in file:
            # Skip empty lines or comments
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            fields = line.split('\t')
            
            # BED2D file typically has at least 3 columns, handle additional ones if present
            chrom = fields[0]
            chromStart = fields[1]
            chromEnd = fields[2]
            name = fields[3] if len(fields) > 6 else None
            score = fields[4] if len(fields) > 8 else None
            strand = fields[5] if len(fields) > 10 else None
            k = len(fields) // 2 - 1 
            record = BedRecord(chrom, chromStart, chromEnd, name, score, strand)

            chrom2 = fields[k + 1]
            chromStart2 = fields[k + 2]
            chromEnd2 = fields[k + 3]
            name2 = fields[k + 4] if len(fields) > 6 else None
            score2 = fields[k + 5] if len(fields) > 8 else None
            strand2 = fields[k + 6] if len(fields) > 10 else None
            record2 = BedRecord(chrom2, chromStart2, chromEnd2, name2, score2, strand2)
            
            
            # Create a BedRecord object and add it to the list
            records.append(tuple((record, record2)))
    
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

def parse_bed(in_file, default_strand=1, default_name="Locus"):
    """Parses bed file into python dataframe."""
    bed = parse_bed_file(in_file)
    bed_df = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand"])
    
    i = 0
    for rec in bed:
        bed_tmp = {
                    "Name": rec.name if rec.name != None else f"{default_name}_{i}",
                    "Chromosome": rec.chrom,
                    "Start": rec.chromStart,
                    "End": rec.chromEnd,
                    "Strand": 1 if rec.strand == '+' else -1 if rec.strand == '-' else default_strand
                }
        bed_df = bed_df._append(bed_tmp, ignore_index=True)
        if rec.name == None:
            i += 1
    return bed_df

def parse_bed2d(in_file, default_strand=1, default_name="Pair"):
    """Parses bed file into python dataframe."""
    bed2d = parse_bed2d_file(in_file)
    bed2d_df = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand", "Name2", "Chromosome2", "Start2", "End2", "Strand2"])
    
    i = 0
    for rec1, rec2 in bed2d:
        bed2d_tmp = {
                    "Name": rec1.name if rec1.name != None else f"{default_name}_{i}_loc1",
                    "Chromosome": rec1.chrom,
                    "Start": rec1.chromStart,
                    "End": rec1.chromEnd,
                    "Strand": 1 if rec1.strand == '+' else -1 if rec1.strand == '-' else default_strand,
                    "Name2": rec2.name if rec2.name != None else f"{default_name}_{i}_loc2",
                    "Chromosome2": rec2.chrom,
                    "Start2": rec2.chromStart,
                    "End2": rec2.chromEnd,
                    "Strand2": 1 if rec2.strand == '+' else -1 if rec2.strand == '-' else default_strand
                }
        bed2d_df = bed2d_df._append(bed2d_tmp, ignore_index=True)
        if rec1.name == None or rec2.name == None:
            i += 1
    return bed2d_df

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
        # as it requires all gff sequences annotations with tracks, percentage is applied after said annotation
        if percentage != None:
            app_end, quantile = percentage
            match app_end:
                case 'high':
                    quantile = 1 - quantile/100
                case 'low':
                    quantile = quantile/100
            percentage_threshold = np.nanquantile(list(positions['Track']), quantile)
            positions = positions[positions['Track'] >= percentage_threshold].reset_index(drop=True)

    else:
        if percentage != None:
            app_end, quantile = percentage
            match app_end:
                case 'high':
                    quantile = 1 - quantile/100
                case 'low':
                    quantile = quantile/100
            all_values = np.concatenate([bigwig_file.values(chrom, 0, bigwig_file.chroms(chrom)) for chrom in bigwig_file.chroms().keys()])
            percentage_threshold = np.nanquantile(all_values, quantile)

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



        

