from hicue.classes.Reader import *

test_files = [
    "SRP_dd_genes.bed", # bed file
    "TYs.gff", #gff file
    "plants_loops.bed2d" # bed2d file
]

expected_lengths = {
    "SRP_dd_genes.bed": 226, # track file
    "TYs.gff":53, #gff file
    "plants_loops.bed2d": 4 # both file
}

position_path = "test_data/positions/"

annotations = {
    "SRP_dd_genes.bed": {'tracks':"test_data/tracks/PGA_stat.bw"}, # track file
    "TYs.gff":{'gff':"test_data/gff/W303.gff3"}, #gff file
    "plants_loops.bed2d": {'tracks':"test_data/tracks/PGA_stat.bw",'gff':"test_data/gff/D_dadantii.gff"} # both file
}

### Reader tests

# testing simple read
for file in test_files: 
    reader = Reader(position_path + file, file.split('.')[-1])
    positions, _ = reader.read_file(threads=8)
    assert(expected_lengths[file] == len(positions))

reader = Reader(position_path + "TYs.gff", "gff", record_type = "gene")
positions, _ = reader.read_file(threads=8)
assert(len(positions) == 0) # TY are not referenced as genes: none must be returned by the reader

# strict overlap annotation
for file in test_files: 
    reader = Reader(position_path + file, file.split('.')[-1], annotation_files = annotations[file], save_to = f"test_out/classes.test.reader.{file.split('.')[0]}")
    positions, _ = reader.read_file(threads=8)
    assert(expected_lengths[file] == len(positions))
    if "gff" in annotations[file]:
        assert(len(os.listdir(f"test_out/classes.test.reader.{file.split('.')[0]}")) == len(positions)) # checking the ouput of each annotation
    if "tracks" in annotations[file]:
        assert("Tracks" in positions.columns)



