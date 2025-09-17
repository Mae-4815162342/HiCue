from .imports import *

def split_gff(gff_path, outpath = None):
    """Splits a gff file into an ensemble of gff files, one for each gff id, usually chromosomes.
    Writes in a tmp folder created in the gff directory if no outpath is provided.
    Returns the dictionnary of each id associated to its gff path."""
    outdir = f"{'/'.join(gff_path.split('/')[:-1])}"
    folder_name = f"tmp_{'.'.join(gff_path.split('/')[-1].split('.')[:-1])}"
    outdir = (outdir + ("/" if len(outdir) > 0 else "") + folder_name) if not outpath else outpath + folder_name
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        shutil.rmtree(outdir)
        os.mkdir(outdir)
        
    chrom_gff_path = f"{outdir}/{gff_path.split('/')[-1][:-len('.gff')]}"
    chrom_files = {}
    chrom_files_path = {}
    with open(gff_path, 'r') as file:
        header = ""
        take_header = True
        while True:
            line = file.readline()
            if not line:
                break
            if line[0] == "#": 
                if take_header:
                    header +=line
            else:
                take_header = False
                chrom = line[:line.find("\t")]
                if chrom not in chrom_files:
                    chrom_path = f"{chrom_gff_path}_{chrom}.gff"
                    chrom_files_path[chrom] = chrom_path
                    chrom_files[chrom] = open(chrom_path, 'a')
                    chrom_files[chrom].write(header)
                chrom_files[chrom].write(line)
    for file in chrom_files.values():
        file.close()
    return chrom_files_path

def compute_distance(locus1, locus2, center = "start"): # TODO: be aware of duplications in utils_opti
    """Returns the distance in base pairs between two position. None if not in the same chromosome."""
    if locus1["Chromosome"] != locus2["Chromosome"]:
        return None
    match center:
        case "start":
            start1 = min(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else max(locus1["Start"], locus1["End"])
            start2 = min(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else max(locus2["Start"], locus2["End"])
        case "center":
            start1 = (locus1["Start"] + locus1["End"]) // 2
            start2 = (locus2["Start"] + locus2["End"]) // 2
        case "end":
            start1 = max(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else min(locus1["Start"], locus1["End"])
            start2 = max(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else min(locus2["Start"], locus2["End"])
    return start2 - start1